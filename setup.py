from setuptools import setup, find_packages

setup(
    name='gwas_tools',
    version='0.0.3',
    author='1',
    author_email='1',
    description='A GWAS tool',
    packages=find_packages(),
    install_requires=[
        '',
    ],
    classifiers=[
        'Programming Language :: Python :: 3.10',
    ],
    entry_points={
        "console_scripts": [
            "gwas-tools-cli=gwas_tools.gwas_tools_cli:main"
        ],
    },
)
