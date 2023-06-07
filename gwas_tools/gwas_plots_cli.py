from gwas_tools.gwas_tools import *
import argparse
import os
from gwas_tools.__version__ import __version__

def main():
    parser = argparse.ArgumentParser(prog='gwas-plots-cli',
                                     description='Plot GWAS results with calculated statistics.')
    parser.add_argument('--input', '-i', type=str, help='input file containing tab separated statistics', required=True)
    parser.add_argument('--version', '-v', action='version', version='gwas-plots-cli ' + __version__)
    parser.add_argument('--out', '-o', type=str, help='Output filename for results', required=True)
    parser.add_argument('--manhattan', action='store_true', help='Plot manhattan plot')
    parser.add_argument('--qq', action='store_true', help='Plot QQ plot')
    args = parser.parse_args()
    if args.input == '-h' or args.out == '-h':
        parser.print_help()
    else:
        print("Reading statistics file...")
        df = pd.read_csv(args.input, sep='\t')
        if args.manhattan:
            print("Plotting manhattan plot...")
            generate_manhattan_plot(df, args.out)
            exit(0)
        elif args.qq:
            print("Plotting QQ plot...")
            generate_qqplot(df, args.out)
            exit(0)
        else:
            print("Statistics Load dry run successful, please specify a plot type.")
            parser.print_help()
    exit(0)


if __name__ == '__main__':
    main()
