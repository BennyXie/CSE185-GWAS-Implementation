import numpy
import pandas as pd
import numpy as np
import statsmodels.api as sm
import matplotlib.pyplot as plt
import datetime
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor

__MULTIPROCESSING__ = ThreadPoolExecutor
class InvalidFileFormatError(Exception):
    pass


class GenoPhenoMismatch(Exception):
    pass

class InvalidFilterError(Exception):
    pass

def read_geno(vcf_file: str):
    """
    Read genotype file in VCF format

    :param vcf_file: genotype file in VCF format
    :return: genotype data with sample ID and numeric genotype values
    """
    column_name = []
    mapping = {'0/0': numpy.int8(1),
               '1/0': numpy.int8(2),
               '0/1': numpy.int8(2),
               '1/1': numpy.int8(3),
               '0|0': numpy.int8(1),
               '1|0': numpy.int8(2),
               '0|1': numpy.int8(2),
               '1|1': numpy.int8(3)
               }
    with open(vcf_file, 'r') as file:
        for line in file:
            if line.startswith(
                    '#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT'):
                column_name = line.split('\t')
                column_name[0] = column_name[0][1:]  # remove #
                column_name[-1] = column_name[-1][:-1]  # remove \n
                break

    if len(column_name) == 0:
        raise InvalidFileFormatError("Invalid VCF file format")

    try:
        vcf: pd.DataFrame = pd.read_csv(vcf_file, comment="#", sep="\t", names=column_name)
    except pd.errors.ParserError as e:
        raise InvalidFileFormatError("Invalid VCF file format", str(e))

    if vcf.shape[0] == 0:
        raise InvalidFileFormatError("Invalid VCF file format")

    # map genotypes to numeric values 1,2,3
    for i in range(9, len(column_name)):
        vcf[column_name[i]] = vcf[column_name[i]].map(mapping)
    vcf[column_name[1]] = vcf[column_name[1]].astype(numpy.int32)
    return vcf


def read_pheno(pheno_file: str):
    """
    Read phenotype file
    :param pheno_file: phenotype file path
    :return: phenotype data with sample ID and float phenotype measurements
    """

    column_name = ["ID", "PHENO"]

    try:
        pheno: pd.DataFrame = pd.read_csv(pheno_file, sep="\t", names=column_name)
    except pd.errors.ParserError as e:
        raise InvalidFileFormatError("Invalid phenotype file format", str(e))

    if pheno.shape[0] == 0:
        InvalidFileFormatError("Invalid phenotype file format")

    try:
        pheno[column_name[1]] = pheno[column_name[1]].astype(float)
    except ValueError as e:
        raise InvalidFileFormatError("Invalid phenotype file format, "
                                     "only numerical phenotypes are accepted.", str(e))

    transpose = pheno.transpose()
    transpose.columns = transpose.iloc[0]
    return transpose


def verify_geno_pheno(geno: pd.DataFrame, pheno: pd.DataFrame):
    """
    Verify that the genotype and phenotype files match in number of samples and sample IDs
    :param geno: genotype data
    :param pheno: phenotype data
    :return: True if genotype and phenotype files match
    """
    # column mismatch
    if geno.shape[1] - 9 != pheno.shape[0]:
        raise GenoPhenoMismatch("Number of samples in genotype and phenotype files do not match")
    for i in range(9, geno.shape[1]):
        if geno.columns[i] != pheno.iloc[i - 9, 0]:
            raise GenoPhenoMismatch("Sample IDs in genotype and phenotype files do not match")
    return True


def generate_qqplot(data: pd.DataFrame, out: str = None):
    """
    Generate a QQ plot of the expected and observed -log10(p-values)
    :param out: output file path
    :param data: p-value data
    """
    if 'PVALUE' not in data.columns:
        raise ValueError("Input DataFrame should contain a 'pvalues' column.")
    data.dropna(inplace=True)

    num_pvalues = len(data)
    expected_neglogp = -np.log10(np.random.uniform(0, 1, num_pvalues))
    expected_neglogp = np.sort(expected_neglogp)
    # Calculate observed -log10(p) values
    pvalues = data['PVALUE'].values
    observed_neglogp = np.where(pvalues != 0, -np.log10(pvalues), 0)
    observed_neglogp = np.sort(observed_neglogp)

    # Create a new DataFrame with observed and expected -log10(p) values- might have to generate expected values (?)
    sorted_data = pd.DataFrame({'Expected -log10(p)': expected_neglogp,
                                'Observed -log10(p)': observed_neglogp})

    # Plot the expected vs observed -log10(p) values
    plt.figure(figsize=(8, 6))
    plt.scatter(sorted_data['Expected -log10(p)'], sorted_data['Observed -log10(p)'],
                color='b', alpha=0.5)
    plt.plot([min(sorted_data['Observed -log10(p)']), max(sorted_data['Observed -log10(p)'])],
             [min(sorted_data['Observed -log10(p)']), max(sorted_data['Observed -log10(p)'])],
             color='r', linestyle='--')
    plt.xlabel('Expected -log10(p)')
    plt.ylabel('Observed -log10(p)')
    plt.title('GWAS QQ Plot')
    if out is not None:
        plt.savefig(out)
    else:
        plt.show()


def generate_manhattan_plot(geno_with_stats: pd.DataFrame, out=None):
    """
    Generate a Manhattan plot of the expected and observed -log10(p-values)
    :param out: output file path
    :param geno_with_stats: p-value data and chromosome data (DataFrame)
    """
    if 'PVALUE' not in geno_with_stats.columns:
        raise InvalidFilterError("Input DataFrame should contain a 'PVALUE' column.")
    if 'CHROM' not in geno_with_stats.columns:
        raise InvalidFilterError("Input DataFrame should contain a 'CHROM' column.")
    p_values = geno_with_stats['PVALUE'].values.astype(float)
    chromosome_data = geno_with_stats['CHROM'].values.astype(int)
    positions = geno_with_stats['POS'].values.astype(int)
    # Calculate -log10(p) values
    p_values = -np.log10(p_values)
    sorted_indices = np.argsort(positions)
    sorted_chromosome_data = chromosome_data[sorted_indices]
    sorted_p_values = p_values[sorted_indices]
    sorted_positions = positions[sorted_indices]
    unique_chromosomes = sorted(geno_with_stats['CHROM'].unique())
    colors = {0: 'brown', 1: 'red', 2: 'orange', 3: 'yellow', 4: 'green', 5: 'blue', 6: 'purple', 7: 'pink', 8: 'gray',
              9: 'indigo'}

    assigned_colors = []
    for chrom in sorted_chromosome_data:
        assigned_colors.append(colors.get(chrom % 10, 'black'))

    # Calculate the positions of the chromosome labels
    chrom_pos = []
    for i in unique_chromosomes:
        low = np.where(sorted_chromosome_data == i)[0][0]
        high = np.where(sorted_chromosome_data == i)[-1][-1]
        low_pos = sorted_positions[low]
        high_pos = sorted_positions[high]
        for j in range(low, high+1):
            chrom_pos.append(i + (sorted_positions[j]-low_pos)/(high_pos-low_pos) - 0.5)


    fig, ax = plt.subplots()

    ax.scatter(chrom_pos, sorted_p_values, c=assigned_colors, alpha=0.5, s=10, marker='s')

    ax.axhline(-np.log10(0.05), color='red', linestyle='--')
    ax.set_xlim(0, max(chromosome_data) + 1)
    ax.set_xticks(unique_chromosomes)
    ax.set_xticklabels(unique_chromosomes)

    ax.set_xlabel('Chromosome')  # Use 'set_xlabel' instead of 'set_x_label'
    ax.set_ylabel('-log10(p-value)')  # Use 'set_ylabel' instead of 'set_y_label'
    ax.set_title('Manhattan Plot')
    if out is not None:
        fig.savefig(out)
    else:
        plt.show()


def write_stats(genotypes_with_stats: pd.DataFrame, out: str):
    """
    Write a stats.csv file containing the statistics

    :param genotypes_with_stats: dataframe containing genotypes and statistics
    :param output_folder: output file path
    :return: true if successful
    """
    output_df = pd.DataFrame({
        "CHROM": genotypes_with_stats["CHROM"],
        "ID": genotypes_with_stats["ID"],
        "SLOPE": genotypes_with_stats["SLOPE"],
        "INTERCEPT": genotypes_with_stats["INTERCEPT"],
        "RVALUE": genotypes_with_stats["RVALUE"],
        "PVALUE": genotypes_with_stats["PVALUE"],
        "STDERR": genotypes_with_stats["STDERR"]
    })
    output_df.dropna(inplace=True)
    output_df.to_csv(out, sep="\t", index=False)

def __calc_row_stats(row:np.array, y:np.array):
        x = row  # Exclude th first column
        X = sm.add_constant(x)
        model = sm.OLS(y, X)
        results = model.fit()
        if np.count_nonzero(x == x[0]) == len(x):
            return {
                'SLOPE': 0,
                'INTERCEPT': 0,
                'RVALUE': results.rsquared,
                'PVALUE': results.f_pvalue,
                'STDERR': results.bse[0]
            }
        return {
            'SLOPE': results.params[1],
            'INTERCEPT': results.params[0],
            'RVALUE': results.rsquared,
            'PVALUE': results.f_pvalue,
            'STDERR': results.bse[0]
        }

def calc_stats(genotypes: pd.DataFrame, phenotypes: pd.DataFrame, threads: int = 1):
    """
    Calculate statistics for GWAS, iDs that doesn't match are excluded in calculation

    :param phenotypes: phenotype data
    :param genotypes: genotype data
    :return: dataframe with genotypes and statistics
    """

    # create a copy of the genotypes dataframe to store the genotypes with statistics
    # add debug flag


    
    genotypes = genotypes.dropna().copy()
    phenotypes = phenotypes.dropna().copy()
    common_columns = genotypes.columns.intersection(phenotypes.columns)
    genotypes_matched = genotypes[common_columns]
    phenotypes_matched = phenotypes[common_columns]

    # Convert to numpy arrays
    genotypes_arr = genotypes_matched.values
    y = phenotypes_matched.iloc[1].to_numpy()
    if np.count_nonzero(y == y[0]) != len(y):

        slopes = []
        intercepts = []
        rvalues = []
        pvalues = []
        stderrs = []
        # calculate the statistics
        with __MULTIPROCESSING__(max_workers=threads) as executor:
            results = executor.map(__calc_row_stats, genotypes_arr, [y] * len(genotypes_arr))

            for result in results:
                    slopes.append(result['SLOPE'])
                    intercepts.append(result['INTERCEPT'])
                    rvalues.append(result['RVALUE'])
                    pvalues.append(result['PVALUE'])
                    stderrs.append(result['STDERR'])

        genotypes = genotypes.assign(
            SLOPE=pd.Series(slopes),
            INTERCEPT=pd.Series(intercepts),
            RVALUE=pd.Series(rvalues),
            PVALUE=pd.Series(pvalues),
            STDERR=pd.Series(stderrs)
        )

    return genotypes

def __filter_row(row: np.ndarray, maf, mac, total_allele_count):
        allele_count = [0, 0]
        for i in [1, 2, 3]:
            count = np.count_nonzero(row == i)
            if i == 1:
                allele_count[0] += count * 2
            elif i == 2:
                allele_count[0] += count
                allele_count[1] += count
            elif i == 3:
                allele_count[1] += count * 2

        if min(allele_count[0], allele_count[1]) / total_allele_count < maf:
            return True
        if min(allele_count[0], allele_count[1]) < mac:
            return True
        return False
def filter_maf_mac(geno: pd.DataFrame, maf: float = 0, mac: int = 0, threads: int = 1):
    if maf < 0 or maf > 1:
        raise ValueError("maf must be between 0 and 1")
    if mac < 0:
        raise ValueError("mac filter only accepts positive number of minor alleles")

    total_allele_count = (geno.shape[1] - 9) * 2
    with __MULTIPROCESSING__(max_workers=threads) as executor:
        rows = [row[9:] for row in geno.values]
        maf = [maf for i in range(len(rows))]
        mac = [mac for i in range(len(rows))]
        total_allele_count = [total_allele_count for i in range(len(rows))]
        results = list(executor.map(__filter_row, rows, maf , mac, total_allele_count))
    drop = [not result for result in results]
    geno = geno[drop]
    return geno
def run_gwas(phenotypes: pd.DataFrame, genotypes: pd.DataFrame, out: str = None, maf=None, mac=None, threads: int = 1, process = False):
    """
    Run GWAS analysis on phenotypes and genotypes
    :param phenotypes: phenotype file
    :param genotypes: genotype file
    :param out: output directory
    :param maf: minor allele frequency between 0 and 1
    :param mac: minimum sample count
    :return: vcf data and statistics of linear regression
    """
    if process == True:
        __MULTIPROCESSING__ = ProcessPoolExecutor
    print("Running GWAS")
    # filter genotypes
    print("Filtering genotypes")
    if mac is not None and maf is not None:
        filter_maf_mac(genotypes, maf, mac, threads)
    elif maf is not None:
        genotypes = filter_maf_mac(genotypes, maf=maf, threads=threads)
    elif mac is not None:
        genotypes = filter_maf_mac(genotypes, mac, threads)
    # calculate statistics
    print("Calculating statistics")
    geno_with_stats = calc_stats(genotypes, phenotypes, threads)
    # write statistics to file
    current_time = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S-")
    print("Writing statistics to file")
    if out is not None:
        write_stats(geno_with_stats, out + current_time + "stats.csv")
        # plot manhattan plot
        generate_manhattan_plot(geno_with_stats, out + current_time + "manhattan-plot.png")
        # plot qq plot
        generate_qqplot(geno_with_stats, out + current_time + "qq-plot.png")
    else:
        # plot manhattan plot
        generate_manhattan_plot(geno_with_stats)
        # plot qq plot
        generate_qqplot(geno_with_stats)
    return geno_with_stats


if __name__ == '__main__':
    print("run gwas_tools-cli from command line")
