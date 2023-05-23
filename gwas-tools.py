import pandas as pd
import numpy as np
import csv
import stats
import argparse


class InvalidFileFormatError(Exception):
    pass

class GenoPhenoMismatch(Exception):
    pass

def read_geno(vcf_file: str):
    column_name = []
    mapping = {"0|0": 1, "1|0": 2, "0|1": 2, "1|1": 3}

    # TODO(Yifei Ding): need error handling for open file
    with open(vcf_file, 'r') as file:
        for line in file:
            if line.startswith(
                    '#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT'):
                column_name = line.split('\t')
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
    #TODO: need error handling for outside mapping genotype
    for i in range(9, len(column_name)):
        vcf[column_name[i]] = vcf[column_name[i]].map(mapping)

    return vcf


def read_pheno(pheno_file: str):
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

    # column mismatch
    if geno.shape[1] - 9 != pheno.shape[0]:
        raise GenoPhenoMismatch("Number of samples in genotype and phenotype files do not match")

    return True

def generate_plot():
    return True


def generate_stats():
    return True


def write_stats():
    return True


def calc_stats(phenotypes: np.vstack, genotypes: np.vstack):
    slopes, intercepts, r_values, p_values, std_errs = stats.linregress(genotypes, phenotypes)
    return {"slopes": slopes, "intercepts": intercepts, "rvalues": r_values, "pvalues": p_values,
            "stderrs": std_errs}


def run_gwas(phenotypes: str, genotypes: str, out: str, maf=None, count=None):
    print("run gwas")


def main():
    parser = argparse.ArgumentParser(prog='gwas-tools',
                                     description='Perform GWAS analysis on phenotypes and genotypes.')
    parser.add_argument('--pheno', '-p', type=str, help='CSV file containing phenotypes. The first column must be \
    sample ID and the second column must be numeric phenotype measurements.', required=True)
    parser.add_argument('--geno', '-g', type=str, help='VCF file containing genotypes',
                        required=True)
    parser.add_argument('--out', '-o', type=str, help='Output directory for results', required=True)
    parser.add_argument('--maf', type=float,
                        help='Minimum minor allele frequency (between 0 and 1)')
    parser.add_argument('--count', type=int, help='Number of samples to include')
    parser.add_argument('--version', action='version', version='gwas-tools 1.0.0')

    args = parser.parse_args()
    if args.pheno == '-h' or args.geno == '-h' or args.out == '-h':
        parser.print_help()
    else:
        run_gwas(args.pheno, args.geno, args.out, args.maf, args.count)


if __name__ == '__main__':
    main()