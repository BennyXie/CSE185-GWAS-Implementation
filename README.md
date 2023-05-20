# CSE185-GWAS-Implementation

An implementation of GWAS, a bioinformatics tool

## Installation


```bash
git clone --depth 1 https://github.com/BennyXie/CSE185-GWAS-Implementation
cd ./CSE185-GWAS-Implementation
pip install -r requirements.txt
pip install .
```

## Usage

```bash
Usage: gwas-tools --pheno <phenotypes.csv> --geno <genotypes.vcf> --out <output_directory> [--maf <min_maf>] [--count <sample_count>] [--help] [--version]

gwas-tools(version 0.0.1) by ___________ 

Perform GWAS analysis on phenotypes and genotypes.

Required parameters:
  --pheno, -p <phenotypes.csv>       CSV file containing phenotypes.
                                     The first column must be sample ID and the
                                     second column must be numeric phenotype
                                     measurements.
  --geno, -g <genotypes.vcf>         VCF file containing genotypes
  --out, -o <output_directory>       Output directory for results

Optional parameters:
  --maf <min_maf>                    Minimum minor allele frequency (between 0 and 1)
  --count <sample_count>             Number of samples to include

Other options:
  --help, -h                         Show this help message and exit
  --version                          Show version information

Note: Required parameters must be provided.
```