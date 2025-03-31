from setuptools import setup, find_packages

setup(
    name="ASODesignPipeline",
    version="0.1",
    packages=find_packages(),
    install_requires=[
        "biopython",
        "pybloomfiltermmap3",  # Use the correct package name
        "numpy",
        "polars",
        "gget",
        
        # Add other dependencies your project needs
    ],
    author="Seyedarash Ayatollahi",
    author_email="ayatosey@b-tu.de",
    description="Pipeline for Antiesense Oligonucleotide Design",
)