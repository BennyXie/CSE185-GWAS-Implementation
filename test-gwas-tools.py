import unittest
import pandas as pd
from gwas_tools import generate_qqplot
import numpy as np

class MyTestCase(unittest.TestCase):
    def test_something(self):
        self.assertEqual(False, False)  # add assertion here
    def test_generate_qqplot_uniform(self):
        data = pd.DataFrame({
            'pvalues': np.random.uniform(0,1,1000)
        })

        # Call the generate_plot function
        generate_qqplot(data)

        # Assert that no errors occurred during the plot generation
        self.assertTrue(True)
    def test_generate_qqplot_exponential(self):
        data = pd.DataFrame({
            'pvalues': np.random.exponential(scale=1,size=1000)
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
            'pvalues': np.polyval(coefficients,x)
        })

        # Call the generate_plot function
        generate_qqplot(data)

        # Assert that no errors occurred during the plot generation
        self.assertTrue(True)

if __name__ == '__main__':
    unittest.main()
