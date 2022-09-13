# TEQUILA-seq: Data Analysis Software

This repository contains scripts for processing and analyzing TEQUILA-seq (**T**ranscript **E**nrichment and **Q**uantification **U**tilizing **I**sothermally **L**inear-**A**mplified probes in conjunction with long-read **seq**uencing) data.

## Table of Contents

* [Overview](#overview)
* [Dependencies](#dependencies)
* [Installation](#installation)
* [Usage](#usage)

## Overview

The scripts contained in this repository include those for processing TEQUILA-seq data as well as those for data visualization and analysis (see below). 

Our data processing scripts work with raw Oxford Nanopore (ONT) signal data (FAST5 format) as input, and these scripts encompass the following steps:
1. **Basecalling**: Raw ONT signal data (FAST5 format) is basecalled into nucleotide sequences (FASTQ format) using Guppy in fast mode.
2. **Alignment**: Basecalled reads are mapped to a user-supplied reference genome using [minimap2](https://github.com/lh3/minimap2) (Li H., *Bioinformatics* 2018) together with user-supplied reference transcript annotations.
3. **Transcript isoform discovery and quantification**: Full-length transcript isoforms are discovered and quantified from long-read RNA-seq alignment files using [ESPRESSO](https://github.com/Xinglab/espresso).

We also have scripts designed to visualize and further characterize transcript isoforms identified from TEQUILA-seq data. These scripts can perform the following tasks:
1. **Visualize discovered transcript isoforms**: Given a collection of samples subjected to TEQUILA-seq, we can visualize the structures and relative abundances of all transcript isoforms discovered for a given gene.
2. **Detect group-specific and sample-specific transcript isoforms**: For a collection of samples subjected to TEQUILA-seq, if the samples can be partitioned into different groups, we can identify transcript isoforms with group-specific expression and usage. Similarly, we can also identify transcript isoforms with sample-specific expression and usage.
3. **Characterize alternative splicing events underlying discovered transcript isoforms**: Local differences in transcript structure between a given isoform and the canonical isoform of the corresponding gene are classified into different alternative splicing patterns, including exon skipping, alternative 5' and 3' splice sites, mutually exclusive exons, retained introns, alternative first or last exons, and complex splicing. 
4. **Predict NMD-targeted transcript isoforms**: Transcript isoforms targeted by mRNA nonsense-mediated decay (NMD) are predicted from the set of isoforms discovered from TEQUILA-seq data based on the 50 nt rule.

## Dependencies

To run our scripts, the following dependencies will need to be installed and available on `$PATH`:

* [Snakemake](https://snakemake.readthedocs.io) (v5.31.1) 
* Guppy (must be downloaded manually from the [ONT software download page](https://community.nanoporetech.com/downloads) since a login is required
* [minimap2](https://github.com/lh3/minimap2) (v2.17)
* [SAMtools](http://samtools.sourceforge.net) (v1.9)
* [BLAST](https://www.ncbi.nlm.nih.gov/blast/) (v2.10.1)
* [HMMER](http://hmmer.org/) (v3.3.1)
* [UCSC KentUtils](http://hgdownload.soe.ucsc.edu/admin/exe/)
  + bedGraphToBigWig
  + faToTwoBit
  + twoBitInfo
* [Python](https://www.python.org/) 3.8
  + [NumPy](https://numpy.org/) (v1.20.1)
  + [BeautifulSoup4](https://pypi.org/project/beautifulsoup4/) (v4.8.2)
* [R](https://www.r-project.org/) (v4.0.5)
* [Perl](https://www.perl.org/) (v5.26.2)


## Installation

## Usage

* [Data process](#data-process)
  + [Install](#install)
  + [Usage](#usage)
  + [Configuration](#configuration)
* [Data analysis and visualization](#data-analysis-and-visualization)


## Data process

#### Install

[Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html) must be installed first. Then the detailed steps are as follows:

1. Enter `TEQUILA-seq` folder.
2. Modify `CONDA_ENV_PREFIX` paths in [set_env_vars.sh](set_env_vars.sh)
3. Install other dependencies using conda: `./install`.
4. Manually install following packages:
    + argparse: `conda install -c conda-forge configargparse`
5. Download reference genome sequence file and put it into the 'references' folder:         https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh37_mapping/gencode.v34lift37.annotation.gtf.gz
6. Modify some absolute file paths in [snakemake_config.yaml](snakemake_config.yaml).
    + In this project, we mainly have three sample groups: 'brain sample + SIRVset4', 'SH-SY5Y cell + SIRVset4' and 'BRCA cell lines'.
    + Configure files of all sample groups are pre-defined in the `Config_for_each_system`. When performing the analysis, the desired `snakemake_config.yaml` should be copied into main folder.
    + The `visualization_path` and `conda_wrapper` of `snakemake_config.yaml` should be modified accordingly.

#### Usage

Note: all commands should be run under the conda environment: 
1. Modify `snakemake_config.yaml` accordingly.
2. Submit the job with `sbatch --time=7-0:0:0 ./run` (Base-calling, alignment and isoform identification based on all samples)


#### Configuration

[snakemake_config.yaml](snakemake_config.yaml)
* `long_read_samples`
  + each sample has its own entry which defines the `sample_name`
  + if starting from fast5 files then under `sample_name` set: `fast5_dir` and `guppy_config`: e.g.
```
    samples:
      RNA_sample_1:
        - guppy_config: 'dna_r9.4.1_450bps_fast.cfg'
          fast5_dir: '~/fast5/'
```
  + if starting from a fastq file then under `sample_name` set: `fastq`: e.g.
```
    samples:
      RNA_sample_1:
        - fastq: '~/combined.fastq'
```
* `guppy_bin_path: '/path/to/guppy/bin/'` (only needed if starting from fast5 files)
* Specify the reference files to use:
  + Either manually put each file in `references/` or provide a url to download the (potentially gzipped) file:
  + `gtf_name: 'file_name.gtf'`
  + `gff3_name: 'file_name.gff3'`
  + `fasta_name: 'file_name.fasta'`
* `conda_wrapper:`: based on current folder
* `espresso_path:`: based on current folder
```
reference_files:
  file_name.gtf.gz:
    url: 'protocol://url/for/file_name.gtf.gz'
```


## Data analysis and visualization

### IMPACT genes panel on BRCA cell lines 
Note: all commands should be run under the conda environment: 
1. Manually install following packages:
  + R: `conda config --add channels conda-forge`
  + `conda config --set channel_priority strict`
  + `conda create -n r_plot r-essentials r-base=4.0.5`
  + `conda activate r_plot`
  + ggplot2: `conda install -c r r-ggplot2`
  + tidyverse: `conda install -c r r-tidyverse`
  + ggplotify: `conda install -c conda-forge r-ggplotify`
  + numpy for python: `conda install -c anaconda numpy`
  + scales: `conda install -c r r-scales`
  + forcats: `conda install -c conda-forge r-forcats`
  
2. Generate figures for selected examples.
  + Enter [Examples_visualization](Examples_visualization) folder.
  + Run `sh  Examples_visualization.sh`.
  + Check figures generated in [Example_res](./Examples_visualization/Example_res) folder.
  
