import numpy
import pandas as pd
import numpy as np
import statsmodels.api as sm
import matplotlib.pyplot as plt
import datetime
from concurrent.futures import ThreadPoolExecutor


class InvalidFileFormatError(Exception):
    pass


class GenoPhenoMismatch(Exception):
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

    num_pvalues = len(data)
    expected_neglogp = -np.log10(np.random.uniform(0, 1, num_pvalues))
    expected_neglogp = np.sort(expected_neglogp)
    # Calculate observed -log10(p) values
    observed_neglogp = -np.log10(data['PVALUE'])
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
        raise ValueError("Input DataFrame should contain a 'PVALUE' column.")
    if 'CHROM' not in geno_with_stats.columns:
        raise ValueError("Input DataFrame should contain a '#CHROM' column.")
    p_values = geno_with_stats['PVALUE']
    chromosome_data = geno_with_stats['CHROM']

    # Calculate -log10(p) values
    p_values = -np.log10(p_values)

    sorted_indices = np.argsort(chromosome_data)
    sorted_chromosome_data = np.array(chromosome_data)[sorted_indices]
    sorted_p_values = np.array(p_values)[sorted_indices]

    # unique_chromosomes = sorted_chromosome_data.unique()
    colors = {0: 'brown', 1: 'red', 2: 'orange', 3: 'yellow', 4: 'green', 5: 'blue', 6: 'purple', 7: 'pink', 8: 'gray',
              9: 'indigo'}

    x_ticks = []
    x_tick_labels = []
    prev_chrom = None

    for i, chrom in enumerate(sorted_chromosome_data):
        if prev_chrom is None or chrom != prev_chrom:
            x_ticks.append(i)
            x_tick_labels.append(str(chrom))
            prev_chrom = chrom

        plt.plot(i, sorted_p_values[i], 'o', color=colors.get(chrom % 10, 'black'), alpha=0.5)

    plt.axhline(-np.log10(0.05), color='red', linestyle='--')
    plt.xlabel('Chromosome')
    plt.ylabel('-log10(p-value)')
    plt.title('Manhattan Plot')
    plt.xticks(x_ticks, x_tick_labels)
    if out is not None:
        plt.savefig(out)
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

    output_df.to_csv(out, sep="\t", index=False)


def calc_stats(genotypes: pd.DataFrame, phenotypes: pd.DataFrame, threads: int = 1):
    """
    Calculate statistics for GWAS, iDs that doesn't match are excluded in calculation

    :param phenotypes: phenotype data
    :param genotypes: genotype data
    :return: dataframe with genotypes and statistics
    """

    # create a copy of the genotypes dataframe to store the genotypes with statistics
    # add debug flag

    def calc_row_stats(row):
        x = row  # Exclude the first column
        X = sm.add_constant(x)
        model = sm.OLS(y, X)
        results = model.fit()

        return {
            'SLOPE': results.params[1],
            'INTERCEPT': results.params[0],
            'RVALUE': results.rsquared,
            'PVALUE': results.f_pvalue,
            'STDERR': results.bse[0]
        }
    
    genotypes = genotypes.dropna()
    phenotypes = phenotypes.dropna()
    common_columns = list(set(genotypes.columns) & set(phenotypes.columns))
    genotypes_matched = genotypes[common_columns]
    phenotypes_matched = phenotypes[common_columns]
    # Convert to numpy arrays
    genotypes_arr = genotypes_matched.values
    y = phenotypes_matched.iloc[1].to_numpy()
    slopes = []
    intercepts = []
    rvalues = []
    pvalues = []
    stderrs = []
    # calculate the statistics
    with ThreadPoolExecutor(max_workers=threads) as executor:
        results = executor.map(calc_row_stats, genotypes_arr)

        for result in results:
            slopes.append(result['SLOPE'])
            intercepts.append(result['INTERCEPT'])
            rvalues.append(result['RVALUE'])
            pvalues.append(result['PVALUE'])
            stderrs.append(result['STDERR'])

    genotypes['SLOPE'] = pd.Series(slopes)
    genotypes['INTERCEPT'] = pd.Series(intercepts)
    genotypes['RVALUE'] = pd.Series(rvalues)
    genotypes['PVALUE'] = pd.Series(pvalues)
    genotypes['STDERR'] = pd.Series(stderrs)
    return genotypes


def filter_maf_mac(geno: pd.DataFrame, maf: float = 0, mac: int = 0):
    if maf < 0 or maf > 1:
        raise ValueError("maf must be between 0 and 1")
    if mac < 0:
        raise ValueError("mac filter only accepts positive number of minor alleles")
    # create a freq df which only contains the genotypes
    total_allele_count = (geno.shape[1] - 9) * 2
    drop = []
    for index, row in geno.iterrows():
        allele_count = [0, 0]
        for i in [1, 2, 3]:
            count = row.iloc[9:].isin([i]).sum()
            if i == 1:
                allele_count[0] += count * 2
            elif i == 2:
                allele_count[0] += count
                allele_count[1] += count
            elif i == 3:
                allele_count[1] += count * 2

        if min(allele_count[0], allele_count[1]) / total_allele_count < maf:
            drop.append(index)
        if min(allele_count[0], allele_count[1]) < mac:
            drop.append(index)
    geno.drop(drop, inplace=True)
    return geno


def run_gwas(phenotypes: pd.DataFrame, genotypes: pd.DataFrame, out: str = None, maf=None, mac=None, threads: int = 1):
    """
    Run GWAS analysis on phenotypes and genotypes
    :param phenotypes: phenotype file
    :param genotypes: genotype file
    :param out: output directory
    :param maf: minor allele frequency between 0 and 1
    :param mac: minimum sample count
    :return: vcf data and statistics of linear regression
    """
    if mac is not None and maf is not None:
        filter_maf_mac(genotypes, maf, mac)
    elif maf is not None:
        genotypes = filter_maf_mac(genotypes, maf=maf)
    elif mac is not None:
        genotypes = filter_maf_mac(genotypes, mac)
    # calculate statistics
    geno_with_stats = calc_stats(genotypes, phenotypes, threads)
    # write statistics to file
    current_time = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S-")

    if out is not None:
        write_stats(geno_with_stats, out + current_time + "stats.csv")
    # plot manhattan plot
    generate_manhattan_plot(geno_with_stats, out + current_time + "manhattan-plot.png")
    # plot qq plot
    generate_qqplot(geno_with_stats, out + current_time + "qq-plot.png")
    return geno_with_stats


if __name__ == '__main__':
    print("run gwas_tools-cli from command line")
