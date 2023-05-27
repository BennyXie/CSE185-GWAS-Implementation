import numpy
import pandas as pd
import numpy as np
import argparse
import scipy.stats as stats
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
    # TODO: may need to trim data
    return vcf


def read_pheno(pheno_file: str):
    """
    Read phenotype file
    :param pheno_file: phenotype file path
    :return: phenotype data with sample ID and float phenotype measurements
    """

    column_name = ["ID", "Phenotype"]

    try:
        pheno: pd.DataFrame = pd.read_csv(pheno_file, sep="\t", names=column_name)
    except pd.errors.ParserError as e:
        raise InvalidFileFormatError("Invalid phenotype file format", str(e))

    if pheno.shape[0] == 0:
        InvalidFileFormatError("Invalid phenotype file format")

    # TODO: need error handling for non-numeric phenotype
    pheno[column_name[1]] = pheno[column_name[1]].astype(float)

    return pheno


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


def generate_qqplot(data: pd.DataFrame):
    """
    Generate a QQ plot of the expected and observed -log10(p-values)
    :param data: p-value data
    """
    if 'pvalues' not in data.columns:
        raise ValueError("Input DataFrame should contain a 'pvalues' column.")

    num_pvalues = len(data)
    expected_neglogp = -np.log10(np.random.uniform(0, 1, num_pvalues))
    expected_neglogp = np.sort(expected_neglogp)
    # Calculate observed -log10(p) values
    observed_neglogp = -np.log10(data['pvalues'])
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
    plt.show()

def generate_manhattanplot(data:pd.DataFrame, chromosome_data:pd.DataFrame):
    """
    Generate a Manhattan plot of the expected and observed -log10(p-values)
    :param data: p-value data
    :param chromosome_data: chromosome data
    """
    p_values = data['pvalue']
    chromosome_data = chromosome_data['#CHROM']

    # Calculate -log10(p) values
    p_values = -np.log10(p_values)


    sorted_indices = np.argsort(chromosome_data)
    sorted_chromosome_data = np.array(chromosome_data)[sorted_indices]
    sorted_p_values = np.array(p_values)[sorted_indices]

    #unique_chromosomes = sorted_chromosome_data.unique()
    colors={0: 'brown', 1:'red', 2:'orange', 3:'yellow', 4:'green', 5:'blue', 6:'purple', 7:'pink', 8: 'gray', 9: 'indigo'}

    x_ticks = []
    x_tick_labels = []
    prev_chrom = None

    for i, chrom in enumerate(sorted_chromosome_data):
        if prev_chrom is None or chrom != prev_chrom:
            x_ticks.append(i)
            x_tick_labels.append(str(chrom))
            prev_chrom = chrom
        
        plt.plot(i, sorted_p_values[i], 'o', color=colors.get(chrom%10, 'black'), alpha=0.5)

    
    plt.axhline(-np.log10(0.05), color='red', linestyle='--')
    plt.xlabel('Chromosome')
    plt.ylabel('-log10(p-value)')
    plt.title('Manhattan Plot')
    plt.xticks(x_ticks, x_tick_labels)
    plt.show()



def generate_stats():
    return True


def write_stats():
    return True


def calc_stats(phenotypes: np.vstack, genotypes: np.vstack):
    ## TODO: Use dataframe instead of vstack
    """
    Calculate statistics for GWAS

    :param phenotypes: phenotype data
    :param genotypes: genotype data
    :return: statistics of linear regression
    """
    slopes, intercepts, r_values, p_values, std_errs = stats.linregress(genotypes, phenotypes)
    return {"slopes": slopes, "intercepts": intercepts, "rvalues": r_values, "pvalues": p_values,
            "stderrs": std_errs}


def filter_maf(geno: pd.DataFrame, maf: float):
    """
    Filter genotypes by minor allele frequency
    :param geno: genotype data
    :param maf: minor allele frequency between 0 and 1
    :return: genotype data with minor allele frequency greater than maf
    """
    # create a freq df which only contains the genotypes
    freq = geno.drop(columns=geno.columns[0:9])
    total_allele_count = len(freq.iloc[0]) * 2
    
    # determine which row to remove
    for index, row in freq.iterrows():
        minor_allele_count = 0
        for i in row:
            minor_allele_count += i - 1
        if minor_allele_count / total_allele_count < maf:
            geno = geno.drop(index)
    return geno


def filter_count(geno: pd.DataFrame, count: int):
    """
    Filter genotypes by sample count
    :param geno: genotype data
    :param count: count of samples
    :return: genotype data with sample count greater than count
    """

    return geno


def run_gwas(phenotypes: str, genotypes: str, out: str, maf=None, count=None):
    """
    Run GWAS analysis on phenotypes and genotypes
    :param phenotypes: phenotype file
    :param genotypes: genotype file
    :param out: output directory
    :param maf: minor allele frequency between 0 and 1
    :param count: minimum sample count
    :return: vcf data and statistics of linear regression
    """
    # read genotype file
    geno = read_geno(genotypes)
    # read phenotype file
    pheno = read_pheno(phenotypes)
    ## TODO: filter and generate graphs
    stats = calc_stats(pheno, geno)

    return geno, stats

if __name__ == '__main__':
    print("run gwas-tools-cli from command line")