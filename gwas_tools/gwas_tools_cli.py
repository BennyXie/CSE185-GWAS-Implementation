from gwas_tools.gwas_tools import *
import argparse
import os
from gwas_tools.__version__ import __version__

def main():
    parser = argparse.ArgumentParser(prog='gwas_tools-cli',
                                     description='Perform GWAS analysis on phenotypes and genotypes.')
    parser.add_argument('--pheno', '-p', type=str, help='CSV file containing phenotypes. The first column must be \
        sample ID and the second column must be numeric phenotype measurements.', required=True)
    parser.add_argument('--geno', '-g', type=str, help='VCF file containing genotypes',
                        required=True)
    parser.add_argument('--out', '-o', type=str, help='Output directory for results', required=True)
    parser.add_argument('--maf', type=float,
                        help='Minimum minor allele frequency (between 0 and 1)')
    parser.add_argument('--mac', type=int, help='Minimum of minor allele count of a genotype')
    parser.add_argument('--version', action='version', version='gwas_tools ' + __version__)
    parser.add_argument('--debug', action='store_true', help='Print debug messages')
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
        #Patch output directory
        if args.out[-1] != '/' and args.out[-1] != '\\':
            args.out += '/'
        # create output directory if it doesn't exist
        if not os.path.exists(args.out):
            os.makedirs(args.out)

        # read genotype file
        try :
            geno = read_geno(args.geno)
        except InvalidFileFormatError:
            print("Invalid genotype file. It must be a VCF file.")
            parser.print_help()
            exit(-1)
        # read phenotype file
        try:
            pheno = read_pheno(args.pheno)
        except InvalidFileFormatError:
            print("Invalid phenotype file. It must be a CSV file with the first column being sample ID and the second column being numeric phenotype measurements.")
            parser.print_help()
            exit(-1)


        # run GWAS
        if args.debug:
            run_gwas(pheno, geno, args.out, args.maf, args.mac)
            exit(0)
        try:
            run_gwas(pheno, geno, args.out, args.maf, args.mac)
        except ValueError:
            print("Invalid input, likely to be malformed maf or mac value.") #TODO: This may ignore some errors and print incorrect message
            parser.print_help()
            exit(-1)
    exit(0)


if __name__ == '__main__':
    main()
