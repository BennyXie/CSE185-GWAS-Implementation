from gwas_tools import *
import argparse
import os


def main():
    parser = argparse.ArgumentParser(prog='gwas-tools-cli',
                                     description='Perform GWAS analysis on phenotypes and genotypes.')
    parser.add_argument('--pheno', '-p', type=str, help='CSV file containing phenotypes. The first column must be \
        sample ID and the second column must be numeric phenotype measurements.', required=True)
    parser.add_argument('--geno', '-g', type=str, help='VCF file containing genotypes',
                        required=True)
    parser.add_argument('--out', '-o', type=str, help='Output directory for results', required=True)
    parser.add_argument('--maf', type=float,
                        help='Minimum minor allele frequency (between 0 and 1)')
    parser.add_argument('--mac', type=int, help='Minimum of minor allele count of a genotype')
    parser.add_argument('--version', action='version', version='gwas-tools ' + VERSION)

    args = parser.parse_args()
    if args.pheno == '-h' or args.geno == '-h' or args.out == '-h':
        parser.print_help()
    elif not os.path.isfile(args.geno) or not os.path.isfile(args.pheno):
        print("Invalid file path.")
        parser.print_help()
    elif not os.path.isdir(args.out):
        print("Invalid output directory.")
        parser.print_help()
    else:
        # read genotype file
        geno = read_geno(args.geno)
        # read phenotype file
        pheno = read_pheno(args.pheno)

        #Patch output directory
        if args.out[-1] != '/' and args.out[-1] != '\\':
            args.out += '/'
        # create output directory if it doesn't exist
        if not os.path.exists(args.out):
            os.makedirs(args.out)

        # run GWAS
        gwas_stats = run_gwas(pheno, geno, args.out, args.maf, args.mac)



if __name__ == '__main__':
    main()
