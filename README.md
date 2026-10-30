# Table of Contents
- [Getting started](#getting-started)
- [Installation](#installation)
  - [Obtaining variantopia](#obtaining-variantopia)
  - [Prerequisites](#prerequisites)
  - [Installation suggestions](#installation-suggestions)
- [How to cite](#how-to-cite)

# Getting started
```
# Download this repository [making sure you have prerequisite Python packages]
git clone https://github.com/zkstewart/variantopia.git

# Assess an MSA to understand variants indicated by alignment
python /location/of/variantopia/variantopia.py msa plot \
    -i $VCF -o msa_plot.png \
    -s mac # Minor Allele count
python /location/of/variantopia/variantopia.py msa report \
    -i $VCF -o report.tsv \
    -f per_variant

# Plot VCF statistics to understand variants in genes and/or chromosomes
python /location/of/variantopia/variantopia.py plot \
    -i $VCF -o vcf_plot.png \
    -s mac -f chomosomes -w 10000 \
    --genome $GENOME

# Reformat a VCF for manual inspection or compatibility with other software
python /location/of/variantopia/variantopia.py reformat table \
    -i $VCF -o table.tsv

# Produce a comprehensive statistics report of a VCF
python /location/of/variantopia/variantopia.py stats \
    -i $VCF -o stats.tsv
```

# Installation
## Obtaining variantopia
Download variantopia by cloning the repository as below. It is available as a collection of Python scripts, so no further installation or compilation is necessary.

```
git clone --recurse-submodules https://github.com/zkstewart/variantopia.git
```

## Prerequisites
variantopia is written in Python and requires a modern version 3. It makes use of several other Python packages which include:
- biopython (https://biopython.org/)
- cyvcf2 (https://github.com/brentp/cyvcf2)
- matplotlib (https://matplotlib.org/)
- ncls (https://pypi.org/project/ncls)
- numpy (https://numpy.org/)
- pandas (https://pandas.pydata.org/)
- plotly (https://plotly.com/python/)
- pyCirclize (https://github.com/moshi4/pyCirclize)
- pyMSAviz (https://github.com/moshi4/pyMSAviz)
- SciPy (https://scipy.org/)

It also calls upon the external programs **tabix** and **bcftools** (https://www.htslib.org) for various operations involving VCF files.

The `plink` and `vcf to plink` modes call upon **plink2** (https://www.cog-genomics.org/plink/2.0/).

The `saige` mode calls upon **SAIGE** (https://saigegit.github.io/SAIGE-doc/), which comes embedded within this repository.

Plotting within `plink` and `saige` modes calls upon **R** (https://www.r-project.org/) with the packages optparse, data.table, ggplot2, and png.

## Installation suggestions
If you are interested in running variantopia, you should ideally set up an Anaconda or Miniconda environment containing a recent Python 3 version. All Python packages listed above are available through community channels including conda-forge.

The external programs can be installed through the bioconda and conda-forge channels, or you may opt to install them yourself by following any instructions detailed at the provided websites.

See the [Installing variantopia wiki page](https://github.com/zkstewart/variantopia/wiki/Installing-variantopia) for more information.

# How to use variantopia
On the command line, you can always ask variantopia to provide help information for its different functions by doing:

```
python /location/of/variantopia.py msa -h
python /location/of/variantopia.py plot -h
python /location/of/variantopia.py reformat table -h
... etc ...
```

Otherwise, refer to the [variantopia wiki](https://github.com/zkstewart/variantopia/wiki) for more detailed information on the program.

# How to cite
There are no plans to publish this repository as a scientific package. Specific aspects of it may be used in publications, and if you find yourself following along with their Methodology by using variantopia, please cite that paper.

If you find no publications involving this author that describe something in variantopia, and you do find variantopia useful in your own work, please feel free to link to this repository in your manuscript.
