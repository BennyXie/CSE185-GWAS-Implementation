import gwas_tools.gwas_tools as gwas_tools

def main():
    geno = gwas_tools.read_geno('./testfiles/vcf/gwas_test.vcf')
    pheno = gwas_tools.read_pheno('./testfiles/pheno/lab3_gwas.phen')
    gwas_tools.run_gwas(geno, pheno, './test_temp', 0.05, 10)
