from setuptools import setup, find_packages

setup(
    name='pdbParser',
    packages=find_packages(),
    version=1.0,
    author="Ozge Yoluk",
    author_email="ozgeyoluk@proteinart.net",
    description="Python package for preparing ensembles for eBDIMs projections",
    url="https://github.com/ozyo/pdbParser",
    install_requires=["numpy","biopython","requests","argparse"],
    classifiers=[
        "Programming Language :: Python :: 2.7",
        "LICENSE :: OSI APPROVED :: GNU GENERAL PUBLIC LICENSE V3 (GPLV3)",
        "Operating System :: OS Independent",
    ]
)
