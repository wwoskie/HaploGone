# HaploGone: pipeline for loss of heterozygosity (LOH) detection

## Installation

To install this pipeline run:
```bash
git clone https://github.com/wwoskie/HaploGone.git && \
cd HaploGone && \
git checkout wip-1-optimize-vcf-parcer && \
cd code &&\
git clone https://github.com/wwoskie/cbs.git
```
**NOTE:** please, pay attention to submodule initialization

## Requirements

One can install requirements via conda (1) or pip (2):
1. `conda env create -f haplogone.yml`
2. `pip install -r requirements.txt`