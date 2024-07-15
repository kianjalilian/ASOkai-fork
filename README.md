# ASO-design-pipeline

This will be a software package that helps the ASO design process to find optimal target sites based on a gene, it's sequence properties, thermodynamic and kinetic parameters.

## Packages & tools

* pyensembl: python library to connect to ENSEMBL (https://pyensembl.readthedocs.io/en/latest/pyensembl.html)
* Bowtie2: Used for read alignment against genome
* RNAcofold (ViennaRNA): Calculating dG (affinity)


## Status Quo

`main.py` takes a gene id from ENsembl and a k-mer length as input. Genome and annotation release are configurable. Paths to bowtie and pyensemble are hard coded as of right now. Output is a csv file containing a list of distinct oligo candidates (reverse complement of target sites) targeting the input gene.

e.g.: KRAS 16 mers should result in 7402 candidate oligos.


# Binding Multiplicity

an oligo candidate has a reasonable (e.g. ddG=5) binding affinity to some other location on the target gene plus a large enough substring in the middle in common (RNASE H1 activity starts at 5 consecutive matches, best activity 8-10, Review Crooke 2021)

# Multiple Binding sites within one mRNA

Pedersen 2020

Test saturday
