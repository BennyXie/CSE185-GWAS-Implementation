import os.path
import unittest
import pandas as pd
from gwas_tools.gwas_tools import *
import numpy as np
import subprocess

class TestPlots(unittest.TestCase):
    def test_generate_qqplot_uniform(self):
        data = pd.DataFrame({
            'PVALUE': np.random.uniform(0, 1, 1000)
        })

        # Call the generate_plot function
        generate_qqplot(data)

        # Assert that no errors occurred during the plot generation
        self.assertTrue(True)

    def test_generate_qqplot_exponential(self):
        data = pd.DataFrame({
            'PVALUE': np.random.exponential(scale=1, size=1000)
        })

        # Call the generate_plot function
        generate_qqplot(data)

        # Assert that no errors occurred during the plot generation
        self.assertTrue(True)

    def test_generate_qqplot_polynomial(self):
        np.random.seed(1)

        # Generate x-values from 1 to 10
        x = np.linspace(1, 10, 100)
        coefficients = [2, 3, -5, 10]
        data = pd.DataFrame({
            'PVALUE': np.polyval(coefficients, x)
        })

        # Call the generate_plot function
        generate_qqplot(data, "./test_temp/qqplot.png")

        # Assert that no errors occurred during the plot generation
        self.assertTrue(os.path.isfile("./test_temp/qqplot.png"))
    def test_generate_manhattan_plot(self):
        # Create example DataFrame
        data = pd.DataFrame({
            'CHROM': np.random.randint(1, 10, size=1000), 'PVALUE': np.random.uniform(0,1,1000)
        })

        # Call the generate_plot function
        generate_manhattan_plot(data, "./test_temp/manhattan_plot.png")

        # Check if a plot is created
        self.assertTrue(os.path.isfile("./test_temp/manhattan_plot.png"))

    def test_generate_manhattanplot_morechromo(self):
        # Create example DataFrame
        data = pd.DataFrame({
            'PVALUE': np.random.uniform(0,1,1000),
            'CHROM': np.random.randint(1, 20, size=1000)
        })

        # Call the generate_plot function
        generate_manhattan_plot(data)

        # Check if a plot is created
        self.assertTrue(True)
    
    def test_generate_manhattanplot_bigger_size(self):
        # Create example DataFrame for p-values
        data = pd.DataFrame({
            'PVALUE': np.random.uniform(0,1,10000),
            'CHROM': np.random.randint(1, 8, size=10000)
        })

        # Call the generate_plot function
        generate_manhattan_plot(data)

        # Check if a plot is created
        self.assertTrue(True)



class TestReadGeno(unittest.TestCase):
    def test_read_vcf_single_id(self):
        out = read_geno("./testfiles/vcf/single_id.vcf")
        out_dict = out.to_dict()
        self.assertIsInstance(out['POS'][0], numpy.int32)
        self.assertIsInstance(out['A'][0], numpy.int8)
        self.assertEqual(out_dict,
                         {'CHROM': {0: 1},
                          'POS': {0: 2},
                          'ID': {0: 'rs1'},
                          'REF': {0: 'A'},
                          'ALT': {0: '.'},
                          'QUAL': {0: '.'},
                          'FILTER': {0: '.'},
                          'INFO': {0: '.'},
                          'FORMAT': {0: 'GT'},
                          'A': {0: 1}}
                         )

    def test_read_vcf_multiple_id(self):
        out = read_geno("./testfiles/vcf/multi_id.vcf")
        out_dict = out.to_dict()
        self.assertIsInstance(out['POS'][0], numpy.int32)
        self.assertIsInstance(out['A'][0], numpy.int8)
        self.assertIsInstance(out['B'][0], numpy.int8)
        self.assertIsInstance(out['C'][0], numpy.int8)
        self.assertEqual(out_dict,
                         {'CHROM': {0: 1},
                          'POS': {0: 2},
                          'ID': {0: 'rs1'},
                          'REF': {0: 'A'},
                          'ALT': {0: '.'},
                          'QUAL': {0: '.'},
                          'FILTER': {0: '.'},
                          'INFO': {0: '.'},
                          'FORMAT': {0: 'GT'},
                          'A': {0: 1},
                          'B': {0: 2},
                          'C': {0: 3}}
                         )
    def test_read_vcf_multi_line(self):
        out = read_geno("./testfiles/vcf/multi_line.vcf")
        out_dict = out.to_dict()

        self.assertIsInstance(out['POS'][0], numpy.int32)
        self.assertIsInstance(out['A'][0], numpy.int8)
        self.assertIsInstance(out['B'][0], numpy.int8)
        self.assertIsInstance(out['C'][0], numpy.int8)
        self.assertEqual(out_dict,
                         {'CHROM': {0: 1, 1: 2, 2: 3},
                          'POS': {0: 1, 1: 2, 2: 3},
                          'ID': {0: 'rs1', 1: 'rs2', 2: 'rs3'},
                          'REF': {0: 'C', 1: 'G', 2: 'G'},
                          'ALT': {0: '.', 1: '.', 2: '.'},
                          'QUAL': {0: '.', 1: '.', 2: '.'},
                          'FILTER': {0: '.', 1: '.', 2: '.'},
                          'INFO': {0: '.', 1: '.', 2: '.'},
                          'FORMAT': {0: 'GT', 1: 'GT', 2: 'GT'},
                          'A': {0: 1, 1: 3, 2: 2},
                          'B': {0: 1, 1: 3, 2: 2},
                          'C': {0: 1, 1: 3, 2: 2}}
                         )


class TestReadPheno(unittest.TestCase):
    def test_read_pheno(self):
        out = read_pheno("./testfiles/pheno/test.phen")
        out_dict = out.to_dict()

        self.assertIsInstance(out['a']['PHENO'], float)

        self.assertEqual(out_dict,
                         {'a':{'ID':'a', 'PHENO': 1.0}, 'b':{'ID':'b', 'PHENO': 2.0}, 'c':{'ID':'c', 'PHENO': 3.0}}
                         )

class TestFilterMaf(unittest.TestCase):
    def test_filter_maf(self):
        # Reading in the vcf file to out
        vcf_file = "./testfiles/vcf/multi_line_after_readin.vcf"
        with open(vcf_file, 'r') as file:
            column_names = file.readline().rstrip('\n').split("\t")
        out = pd.read_csv(vcf_file, \
            comment="#", sep="\t", names=column_names)
        geno = filter_maf(out, 0.3).to_dict()
        # The expected output
        expected = {0: 1, 3: 4}
        self.assertEqual(geno["POS"], expected)
        
class TestFilterCount(unittest.TestCase):
    def test_filter_count(self):
        # Reading in the vcf file to out
        vcf_file = "./testfiles/vcf/multi_line_after_readin.vcf"
        with open(vcf_file, 'r') as file:
            column_names = file.readline().rstrip('\n').split("\t")
        out = pd.read_csv(vcf_file, \
            comment="#", sep="\t", names=column_names)
        
        geno = filter_mac(out, 2).to_dict()
        # The expected output
        expected = {0: 1, 1 : 2, 3: 4, 4: 5}
        self.assertEqual(geno["POS"], expected)

class GWASTestCase(unittest.TestCase):
    def test_gwas_analysis(self):
        print(os.getcwd())
        command = "python -m gwas_tools.gwas_tools_cli --pheno ./testfiles/pheno/lab3_gwas.phen --geno ./testfiles/vcf/gwas_test.vcf --out ./test_temp --maf 0.05 --mac 10"
        result = subprocess.run(command, shell=True, capture_output=True, text=True)
        
        # Assert that the command ran successfully
        self.assertEqual(result.returncode, 0)

        # Define the expected output file names
        expected_output_files = ["manhattan.png","qqplot.png", "stats.csv"]
        files_in_dir = os.listdir("./test_temp")
        
        # Assert that each expected output file is in the expected directory
        for file_name in expected_output_files:
            self.assertIn(file_name, files_in_dir)


if __name__ == '__main__':
    unittest.main()