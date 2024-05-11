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
**NOTE:** please, pay attention to submodule initialization

## Requirements

The user can install requirements via conda (1) or pip (2):
1. `conda env create -f haplogone.yml`
2. `pip install -r requirements.txt`

## Quickstart

Install this pipeline (see **Installation**), import it and run it with your data. You can see **examples.ipynb** as an example of quick start and for more detailed information of using it. 

## Guideline



## Example of working

An example is in the folder **code** with the name **examples.ipynb**. 

## References

1. Olshen AB, Venkatraman ES, Lucito R, Wigler M. Circular binary segmentation for the analysis of array-based DNA copy number data. Biostatistics. 2004 Oct;5(4):557-72. doi: 10.1093/biostatistics/kxh008. PMID: 15475419.
2. https://github.com/jeremy9959/cbs
