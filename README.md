# HaploGone: pipeline for loss of heterozygosity (LOH) detection

This tool allows the user to detect regions with loss of heterozygosity in the genome. This pipeline is based on circular binary segmentation algorithm[1] and its oython realization[2]. 

## Installation

To install this pipeline run:
```bash
git clone https://github.com/wwoskie/HaploGone.git && \
cd HaploGone && \
git checkout wip-2-add-pycirularize-graphics && \
cd code &&\
git clone https://github.com/wwoskie/cbs.git
```
**NOTE:** please, pay attention to submodule initialization:

git submodule init && \
git submodule update

## Requirements

The user can install requirements via conda (1) or pip (2):
1. `conda env create -f haplogone.yml`
2. `pip install -r requirements.txt`

## Quickstart

Install this pipeline (see **Installation**), import it and run it with your data. You can see **examples.ipynb** as an example of quick start and for more detailed information of using it. 

## Guideline

This pipeline contains a class VCF where different functions for operating on the input VCF data are collected.
Attributes of the class:
```input_file``` - the path to the input VCF file.
```segment_size_threshold``` - a threshold for filtering segments by size, default 1e6.
```path_to_centromeres_bed``` a path to the centromeres data, default "centromeres.bed".
```segmentation_shuffles``` the number of random permutations, default 1000.
```segmentation_p``` - p-value for segmentation, default 0.01.
```validation_shuffles``` ... , default 1000.
```validation_p``` validation p-value, default 0.01.

Methods of the class:
```read``` - converts the input VCF file into pd.DataFrame. The input is a path to the VCF file.
```count_baf``` - calculates the B-allele frequency from the variant data of the input dataset.
```segment_baf``` - segments the input data and can filter the segment obtaint by size.
```plot_chromosomes``` - vizualizes the result of the data segmentation for each chromosome.
```read_bed``` - converts bed file into pd.DataFrame.
```create_bed``` - creates a bed pd.DataFrame from a VCF pd.DataFrame with bounds of filtered segments with a given threshold.

## Example of working

An example is in the folder **code** with the name **examples.ipynb**. 

## References

1. Olshen AB, Venkatraman ES, Lucito R, Wigler M. Circular binary segmentation for the analysis of array-based DNA copy number data. Biostatistics. 2004 Oct;5(4):557-72. doi: 10.1093/biostatistics/kxh008. PMID: 15475419.
2. https://github.com/jeremy9959/cbs
