"""
This file contains the unit tests for the gwas_tools package.

run with: `python -m gwas_tools.test_gwas_tools` from project root
"""

import os.path
import unittest
import subprocess
import shutil
import random
#from gwas_tools import *
from gwas_tools.gwas_tools import *

if os.path.exists("./test_temp"):
    shutil.rmtree("./test_temp")
os.mkdir("./test_temp")

class TestPlots(unittest.TestCase):
    def test_generate_qqplot_uniform(self):
        data = pd.DataFrame({
            'PVALUE': np.random.uniform(0, 1, 1000)
        })
        file_name = "./test_temp/qq-plot-0.png"
        # Call the generate_plot function
        generate_qqplot(data, file_name)

        # Assert that no errors occurred during the plot generation
        self.assertTrue(os.path.isfile(file_name))

    def test_generate_qqplot_exponential(self):
        data = pd.DataFrame({
            'PVALUE': np.random.exponential(scale=1, size=1000)
        })

        file_name = "./test_temp/qq-plot-1.png"
        # Call the generate_plot function
        generate_qqplot(data, file_name)

        # Assert that no errors occurred during the plot generation
        self.assertTrue(os.path.isfile(file_name))

    def test_generate_qqplot_polynomial(self):
        np.random.seed(1)

        # Generate x-values from 1 to 10
        x = np.linspace(1, 10, 100)
        coefficients = [2, 3, -5, 10]
        data = pd.DataFrame({
            'PVALUE': np.polyval(coefficients, x)
        })

        file_name = "./test_temp/qq-plot-3.png"
        # Call the generate_plot function
        generate_qqplot(data, file_name)

        # Assert that no errors occurred during the plot generation
        self.assertTrue(os.path.isfile(file_name))

    def test_generate_manhattan_plot(self):
        # Create example DataFrame
        data = pd.DataFrame({
            'CHROM': np.random.randint(1, 10, size=1000), 'PVALUE': np.random.uniform(0, 1, 1000)
        })
        data['POS'] = data['CHROM'] * 10000 + np.random.randint(1, 10000, size=1000)

        file_name = "./test_temp/manhattan-plot-0.png"
        # Call the generate_plot function
        generate_manhattan_plot(data, file_name)

        # Assert that no errors occurred during the plot generation
        self.assertTrue(os.path.isfile(file_name))

    def test_generate_manhattanplot_morechromo(self):
        # Create example DataFrame
        data = pd.DataFrame({
            'PVALUE': np.random.uniform(0, 1, 1000),
            'CHROM': np.random.randint(1, 20, size=1000),
        })

        data['POS'] = data['CHROM'] * 10000 + np.random.randint(1, 10000, size=1000)

        file_name = "./test_temp/manhattan-plot-1.png"
        # Call the generate_plot function
        generate_manhattan_plot(data, file_name)

        # Assert that no errors occurred during the plot generation
        self.assertTrue(os.path.isfile(file_name))

    def test_generate_manhattanplot_bigger_size(self):
        # Create example DataFrame for p-values
        data = pd.DataFrame({
            'PVALUE': np.random.uniform(0, 1, 10000),
            'CHROM': np.random.randint(1, 8, size=10000)
        })
        data['POS'] = data['CHROM'] * 10000 + np.random.randint(1, 10000, size=10000)

        file_name = "./test_temp/manhattan-plot-2.png"
        # Call the generate_plot function
        generate_manhattan_plot(data, file_name)

        # Assert that no errors occurred during the plot generation
        self.assertTrue(os.path.isfile(file_name))


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
                         {'a': {'ID': 'a', 'PHENO': 1.0}, 'b': {'ID': 'b', 'PHENO': 2.0},
                          'c': {'ID': 'c', 'PHENO': 3.0}}
                         )


class TestFilterMaf(unittest.TestCase):
    def test_filter_maf(self):
        # Reading in the vcf file to out
        test_vcf = pd.DataFrame({'CHROM': [1, 1, 1, 1, 1],
                                 'POS': [1, 2, 3, 4, 5],
                                 'ID': ['rs1', 'rs2', 'rs3', 'rs4', 'rs5'],
                                 'REF': ['A', 'A', 'A', 'A', 'A'],
                                 'ALT': ['C', 'C', 'C', 'C', 'C'],
                                 'QUAL': ['.', '.', '.', '.', '.'],
                                 'FILTER': ['.', '.', '.', '.', '.'],
                                 'INFO': ['.', '.', '.', '.', '.'],
                                 'FORMAT': ['GT', 'GT', 'GT', 'GT', 'GT'],
                                 'A': [1, 1, 1, 1, 2],
                                 'B': [3, 1, 1, 1, 2],
                                 'C': [3, 2, 1, 1, 2],
                                 'D': [3, 2, 3, 1, 2],
                                 'E': [3, 2, 3, 2, 1],
                                 })
        geno = filter_maf_mac(test_vcf, maf=0.3).to_dict()
        # The expected output
        expected = ['rs2', 'rs3', 'rs5']
        self.assertEqual(expected, list(geno['ID'].values()))


class TestFilterCount(unittest.TestCase):
    def test_filter_count(self):
        test_vcf = pd.DataFrame({'CHROM': [1, 1, 1, 1, 1],
                                 'POS': [1, 2, 3, 4, 5],
                                 'ID': ['rs1', 'rs2', 'rs3', 'rs4', 'rs5'],
                                 'REF': ['A', 'A', 'A', 'A', 'A'],
                                 'ALT': ['C', 'C', 'C', 'C', 'C'],
                                 'QUAL': ['.', '.', '.', '.', '.'],
                                 'FILTER': ['.', '.', '.', '.', '.'],
                                 'INFO': ['.', '.', '.', '.', '.'],
                                 'FORMAT': ['GT', 'GT', 'GT', 'GT', 'GT'],
                                 'A': [1, 1, 1, 1, 2],
                                 'B': [3, 1, 1, 1, 2],
                                 'C': [3, 2, 1, 1, 2],
                                 'D': [3, 2, 3, 1, 2],
                                 'E': [3, 2, 3, 2, 1],
                                 })
        geno = filter_maf_mac(test_vcf, mac = 3).to_dict()
        # The expected output
        expected = ['rs2', 'rs3', 'rs5']
        self.assertEqual(expected, list(geno['ID'].values()))


class GWASTestCase(unittest.TestCase):
    def test_gwas_analysis(self):
        command = "python -m gwas_tools.gwas_tools_cli --pheno ./testfiles/pheno/lab3_gwas.phen --geno ./testfiles/vcf/gwas_test.vcf --out ./test_temp --maf 0.05 --mac 10"
        file_count_before = len(os.listdir("./test_temp"))
        result = subprocess.run(command, shell=True, capture_output=True, text=True)

        # Assert that the command ran successfully
        self.assertEqual(0,result.returncode)

        # Define the expected output file names
        expected_output_files = ["manhattan-plot.png", "qq-plot.png", "stats.csv"]
        files_in_dir = os.listdir("./test_temp")
        file_count_after = len(files_in_dir)
        self.assertTrue(file_count_after - file_count_before == len(expected_output_files))

    def test_gwas_analysis_no_args(self):
        command = "python -m gwas_tools.gwas-tools-cli"
        result = subprocess.run(command, shell=True, capture_output=True, text=True)

        # Assert that the command ran successfully
        self.assertEqual(1,result.returncode)
    def test_gwas_analysis_no_pheno(self):
        command = "python -m gwas_tools.gwas-tools-cli --geno testfiles/vcf/gwas_test.vcf --out ./test_temp"
        result = subprocess.run(command, shell=True, capture_output=True, text=True)

        # Assert that the command ran successfully
        self.assertEqual(result.returncode, 1)

class CalcStatsTestCase(unittest.TestCase):
    def test_null_case(self):
        genotype_data = {
            'CHROM': [1, 1, 1, 1, 1],
            'POS': [1, 2, 3, 4, 5],
            'ID': ['rs1', 'rs2', 'rs3', 'rs4', 'rs5'],
            'REF': ['A', 'A', 'A', 'A', 'A'],
            'ALT': ['C', 'C', 'C', 'C', 'C'],
            'QUAL': ['.', '.', '.', '.', '.'],
            'FILTER': ['.', '.', '.', '.', '.'],
            'INFO': ['.', '.', '.', '.', '.'],
            'FORMAT': ['GT', 'GT', 'GT', 'GT', 'GT'],
            'A': [random.randint(1, 3) for _ in range(5)],
            'B': [random.randint(1, 3) for _ in range(5)],
            'C': [random.randint(1, 3) for _ in range(5)],
            'D': [random.randint(1, 3) for _ in range(5)],
            'E': [random.randint(1, 3) for _ in range(5)],
        }
        genotypes = pd.DataFrame(genotype_data)

        # Create PHEN-formatted phenotype dataframe
        phenotypes_data = {
            'A': ['A', random.randint(1, 10)],
            'B': ['B', random.randint(1, 10)],
            'C': ['C', random.randint(1, 10)],
            'D': ['D', random.randint(1, 10)],
            'E': ['E', random.randint(1, 10)]
        }
        phenotypes = pd.DataFrame(phenotypes_data)
        result = calc_stats(genotypes, phenotypes)
    
        self.assertTrue(sum(result['SLOPE'])<1 and sum(result['SLOPE'])>-1)


    def test_slope_case(self):
        # Generate random numbers for genotypes and phenotypes with a correlation
        genotype_data = {
            'CHROM': [1, 1, 1, 1, 1],
            'POS': [1, 2, 3, 4, 5],
            'ID': ['rs1', 'rs2', 'rs3', 'rs4', 'rs5'],
            'REF': ['A', 'A', 'A', 'A', 'A'],
            'ALT': ['C', 'C', 'C', 'C', 'C'],
            'QUAL': ['.', '.', '.', '.', '.'],
            'FILTER': ['.', '.', '.', '.', '.'],
            'INFO': ['.', '.', '.', '.', '.'],
            'FORMAT': ['GT', 'GT', 'GT', 'GT', 'GT'],
            'A': [1, 1, 1, 1, 2],
            'B': [3, 1, 1, 1, 2],
            'C': [3, 2, 1, 1, 2],
            'D': [3, 2, 3, 1, 2],
            'E': [3, 2, 3, 2, 1],
        }
        genotypes = pd.DataFrame(genotype_data)

        # Create PHEN-formatted phenotype dataframe
        phenotypes_data = {
            'A':['A',2], 'B':['B',2], 'C':['C',2], 'D':['D',2], 'E': ['E',4]
        }
        phenotypes = pd.DataFrame(phenotypes_data)

        result = calc_stats(genotypes, phenotypes)

        # Check that the slope is non-zero
        self.assertTrue(any(result['SLOPE'] != 0.0))


if __name__ == '__main__':
    unittest.main()
