# CSE185-GWAS-Implementation

An implementation of GWAS, a bioinformatics tool

## Installation


```bash
git clone --depth 1 https://github.com/BennyXie/CSE185-GWAS-Implementation
cd ./CSE185-GWAS-Implementation
pip install -r requirements.txt
pip install . # For user installation, add --user
```



## Usage

```bash
usage: gwas_tools-cli [-h] --pheno PHENO --geno GENO --out OUT [--maf MAF] [--mac MAC] [--version] [--debug] [--threads THREADS] [--process]

Perform GWAS analysis on phenotypes and genotypes.

options:
  -h, --help            show this help message and exit
  --pheno PHENO, -p PHENO
                        CSV file containing phenotypes. The first column must be sample ID and the second column must be numeric phenotype measurements.
  --geno GENO, -g GENO  VCF file containing genotypes
  --out OUT, -o OUT     Output directory for results
  --maf MAF             Minimum minor allele frequency (between 0 and 1)
  --mac MAC             Minimum of minor allele count of a genotype
  --version             show program's version number and exit
  --debug               Print debug messages
  --threads THREADS, -t THREADS
                        Number of threads to use
  --process             Use process instead of thread, faster but more memory usage
```

## Report Error

If there is any error in using our tool, please help us fix the error by leaving a 
comment in the "Issues" section in our repo. We will get to you as soon as possible.
