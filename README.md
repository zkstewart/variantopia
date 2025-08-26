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
git clone https://github.com/zkstewart/variantopia.git
```

## Prerequisites
variantopia is written in Python and requires a modern version 3. It makes use of several other Python packages which include:
- biopython (https://biopython.org/)
- numpy (https://numpy.org/)
- pandas (https://pandas.pydata.org/)
- matplotlib (https://matplotlib.org/)
- cyvcf2 (https://github.com/brentp/cyvcf2)
- ncls (https://pypi.org/project/ncls)

It also calls upon the external programs **tabix** and **bcftools** (https://www.htslib.org)

variantopia has been developed on Linux and within Windows Subsystem for Linux (WSL); it is not currently compatible with plain Windows.

## Installation suggestions
If you are interested in running variantopia, you should ideally set up an Anaconda or Miniconda environment containing a recent Python 3 version. All Python packages listed above are available through community channels including conda-forge.

The external programs can be installed through the bioconda channel, or you may opt to install them yourself by following any instructions detailed at the provided websites.

See the [Installing variantopia wiki page](https://github.com/zkstewart/variantopia/wiki/Installing-variantopia) for more information.

# How to use variantopia
On the command line, you can always ask variantopia to provide help information for its different functions by doing:

```
python /location/of/variantopia.py msa -h
python /location/of/variantopia.py plot -h
python /location/of/variantopia_post.py reformat table -h
... etc ...
```

Otherwise, refer to the [variantopia wiki](https://github.com/zkstewart/variantopia/wiki) for more detailed information on the program.

# How to cite
There are no plans to publish any of the code in this repository. However, if you do find any of this useful in your own work, please feel free to link to this repository in your manuscript.
