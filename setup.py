from setuptools import setup, find_packages
from gwas_tools.__version__ import __version__
setup(
    name='gwas_tools',
    version=__version__,
    author='Yifei Ding, Benny Xie, Marina Xu',
    author_email='yifeiding@pm.me, max003@ucsd.edu, ',
    description='A simplified GWAS tool for analyzing VCF files',
    packages=find_packages(),
    install_requires=[
        'pandas',
        'numpy',
        'matplotlib',
        'statsmodels',
    ],
    classifiers=[
        'Programming Language :: Python :: 3.10',
        'License :: OSI Approved :: MIT License',
        'Operating System :: Linux',
    ],
    entry_points={
        "console_scripts": [
            "gwas-tools-cli=gwas_tools.gwas_tools_cli:main",
            "gwas-plots-cli=gwas_tools.gwas_plots_cli:main"
        ],
    },
)
