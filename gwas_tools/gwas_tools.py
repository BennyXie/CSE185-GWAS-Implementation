import numpy
import pandas as pd
import numpy as np
import statsmodels.api as sm
import matplotlib.pyplot as plt

VERSION = "1.0.0"


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
    vcf[column_name[1]] = vcf[column_name[1]].astype(int)
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
    :param statistics: p-value data
    :param chromosome_data: chromosome data
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


def write_stats(genotypes_with_stats: pd.DataFrame, output_folder: str):
    """
    Write a stats.csv file containing the statistics

    :param genotypes_with_stats: dataframe containing genotypes and statistics
    :param output_folder: output folder
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

    output_df.to_csv(output_folder + "stats.csv", sep="\t", index=False)


def calc_stats(genotypes: pd.DataFrame, phenotypes: pd.DataFrame):
    """
    Calculate statistics for GWAS, iDs that doesn't match are excluded in calculation

    :param phenotypes: phenotype data
    :param genotypes: genotype data
    :return: dataframe with genotypes and statistics
    """
    # create a copy of the genotypes dataframe to store the genotypes with statistics
    # add debug flag
    genotypes_with_stats = genotypes.copy()
    common_columns = list(set(genotypes.columns) & set(phenotypes.columns))
    genotypes_matched = genotypes[common_columns]
    phenotypes_matched = phenotypes[common_columns]
    slopes = []
    intercepts = []
    rvalues = []
    pvalues = []
    stderrs = []
    # calculate the statistics
    for i, row in genotypes_matched.iterrows():
        y = phenotypes_matched.iloc[1].to_numpy()
        x = row.to_numpy()
        X = sm.add_constant(x)
        model = sm.OLS(y, X)
        results = model.fit()
        # add the statistics as columns to the genotypes dataframe
        slopes.append(results.params[1])
        intercepts.append(results.params[0])
        rvalues.append(results.rsquared)
        pvalues.append(results.f_pvalue)
        stderrs.append(results.bse[0])

    genotypes_with_stats['SLOPE'] = pd.Series(slopes)
    genotypes_with_stats['INTERCEPT'] = pd.Series(intercepts)
    genotypes_with_stats['RVALUE'] = pd.Series(rvalues)
    genotypes_with_stats['PVALUE'] = pd.Series(pvalues)
    genotypes_with_stats['STDERR'] = pd.Series(stderrs)
    return genotypes_with_stats


def filter_maf(geno: pd.DataFrame, maf: float):
    """
    Filter genotypes by minor allele frequency
    :param geno: genotype data
    :param maf: minor allele frequency between 0 and 1
    :return: genotype data with minor allele frequency greater than maf
    """
    if maf < 0 or maf > 1:
        raise ValueError("maf must be between 0 and 1")
    # create a freq df which only contains the genotypes
    freq = geno.drop(columns=geno.columns[0:9])
    total_allele_count = len(freq.iloc[0]) * 2

    for index, row in freq.iterrows():
        value_counts = row.value_counts().to_dict()
        minor_allele_count = 0
        # counting maf
        for type, allele_count in value_counts.items():
            minor_allele_count += (type - 1) * allele_count
        # dropping the rows that has too few maf
        if minor_allele_count / total_allele_count < maf:
            geno = geno.drop(index)

    return geno


def filter_mac(geno: pd.DataFrame, mac: int):
    """
    Filter genotypes by sample count
    :param geno: genotype data
    :param mac: count of minor alleles
    :return: genotype data with minor alleles occuring at least 'mac' times.
    """
    if mac < 0:
        raise ValueError("mac filter only accepts positive number of minor alleles")
    # create a freq df which only contains the genotypes
    freq = geno.drop(columns=geno.columns[0:9])
    
    for index, row in freq.iterrows():
        value_counts = row.value_counts().to_dict()
        minor_allele_count = 0
        # counting maf
        for type, allele_count in value_counts.items():
            minor_allele_count += (type - 1) * allele_count
        # dropping the rows that has too few maf
        if minor_allele_count < mac:
            geno = geno.drop(index)
    return geno


def run_gwas(phenotypes: pd.DataFrame, genotypes: pd.DataFrame, out: str = None, maf=None, mac=None):
    """
    Run GWAS analysis on phenotypes and genotypes
    :param phenotypes: phenotype file
    :param genotypes: genotype file
    :param out: output directory
    :param maf: minor allele frequency between 0 and 1
    :param mac: minimum sample count
    :return: vcf data and statistics of linear regression
    """

    if maf is not None:
        genotypes = filter_maf(genotypes, maf)
    if mac is not None:
        genotypes = filter_mac(genotypes, mac)
    # calculate statistics
    geno_with_stats = calc_stats(genotypes, phenotypes)
    # write statistics to file
    if out is not None:
        write_stats(geno_with_stats, out)
    # plot manhattan plot
    generate_manhattan_plot(geno_with_stats, out+"manhattan.png")
    # plot qq plot
    generate_qqplot(geno_with_stats, out+"qqplot.png")
    return geno_with_stats

if __name__ == '__main__':
    print("run gwas_tools-cli from command line")