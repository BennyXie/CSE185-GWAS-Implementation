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
usage: gwas-tools-cli [-h] --pheno PHENO --geno GENO --out OUT [--maf MAF]
                      [--count COUNT] [--version]
Perform GWAS analysis on phenotypes and genotypes.
options:
  -h, --help            show this help message and exit
  --pheno PHENO, -p PHENO
                        CSV file containing phenotypes. The first column must
                        be sample ID and the second column must be numeric
                        phenotype measurements.
  --geno GENO, -g GENO  VCF file containing genotypes
  --out OUT, -o OUT     Output directory for results
  --maf MAF             Minimum minor allele frequency (between 0 and 1)
  --count COUNT         Number of samples to include
  --version             show program's version number and exit
Process finished with exit code 0

```