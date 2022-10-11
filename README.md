# TEQUILA-seq: Data Analysis Software

This repository contains scripts for processing and analyzing TEQUILA-seq (**T**ranscript **E**nrichment and **Q**uantification **U**tilizing **I**sothermally **L**inear-**A**mplified probes in conjunction with long-read **seq**uencing) data.

## Table of Contents

* [Overview](#overview)
* [Dependencies](#dependencies)
* [Usage](#usage)
  + [Data processing](#data-processing)
  + [Data visualization and analysis](#data-visualization-and-analysis)
    + [Transcript isoform visualization](#transcript-isoform-visualization)
    + [Detection of subtype-specific transcript isoforms](#detection-of-subtype-specific-transcript-isoforms)
    + [Detection of sample-specific transcript isoforms](#detection-of-sample-specific-transcript-isoforms)
    + [Characterization of alternative splicing events underlying discovered transcript isoforms](#characterization-of-alternative-splicing-events-underlying-discovered-transcript-isoforms)
    + [Prediction of NMD-targeted transcript isoforms](#prediction-of-nmd-targeted-transcript-isoforms)

## Overview

The scripts contained in this repository include those for processing TEQUILA-seq data as well as those for data visualization and analysis.

<img src="./files/TEQUILA-seq_Analysis_Workflow.png" width="800"/>

Our data processing scripts are designed to work with raw Oxford Nanopore (ONT) signal data (FAST5 format) as input. However, these scripts can work with data from any long-read sequencing platform (e.g., PacBio) as long as the data is provided in FASTQ/SAM/BAM format. These scripts encompass the following steps:
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
  + [pandas](https://pandas.pydata.org/) (v1.1.4)
  + [SciPy](https://scipy.org/) (v1.5.4)
  + [statsmodels](https://www.statsmodels.org/) (v0.12.2)
  + [NetworkX](https://networkx.org/) (v2.6.3)
  + [BeautifulSoup4](https://pypi.org/project/beautifulsoup4/) (v4.8.2)
  + [ConfigArgParse](https://pypi.org/project/ConfigArgParse/)
* [R](https://www.r-project.org/) (v4.0.5)
  + [ggplot2](https://ggplot2.tidyverse.org/)
  + [tidyverse](https://www.tidyverse.org/)
  + [ggplotify](https://cran.r-project.org/package=ggplotify)
  + [scales](https://scales.r-lib.org/)
  + [forcats](https://forcats.tidyverse.org/)
* [Perl](https://www.perl.org/) (v5.26.2) built with threading enabled
  + Check for thread support with `perl -e 'use threads; print("ok\n")'`

The source code for ESPRESSO (v1.2.2), which is another dependency, is available in the folder [ESPRESSO_alpha1.2.2](./ESPRESSO_alpha1.2.2/). For details on memory requirements (i.e., as a function of input data size and number of threads), please refer to the [ESPRESSO GitHub page](https://github.com/Xinglab/espresso).

## Usage

### Data processing

A Snakemake workflow is provided for running the data processing scripts. The workflow can start from raw ONT signal data (FAST5 format), sequencing reads data (FASTQ format), or read-to-genome alignment data (SAM/BAM format). **Note:** A different configuration will be required depending on what type of data you are starting with. The workflow can also be configured to run multiple samples and each sample can have multiple inputs. Set the configuration by editing the [snakemake_config.yaml](./snakemake_config.yaml) and [snakemake_profile/config.yaml](./snakemake_profile/config.yaml) files. 

Prior to running the data processing scripts, please make sure that [conda](https://conda.io/) is installed. Then, `./install` can be run to install all dependencies required for data processing (except Guppy, which requires a manual installation). Specifically, `./install` will create a conda environment with the required dependencies and set some absolute file paths in [snakemake_config.yaml](./snakemake_config.yaml)

Details on configuring the [snakemake_config.yaml](./snakemake_config.yaml) file are as follows:

* Set the amount of resources to allocate for each task as follows:
  + `{job_name}_threads: {num_threads}`
  + `{job_name}_mem_gb: {num_GBs}`
  + `{job_name}_time_hr: {num_hours}`
* Specify the path to the folder containing the ESPRESSO source code as follows:
  + `espresso_path: /path/to/ESPRESSO/src`
* If starting with raw ONT signal data, specify the path to the `bin` folder containing Guppy as follows:
  + `guppy_bin_path: /path/to/guppy/bin/`
* Specify the reference genome sequence (FASTA format) and reference transcript annotations (GTF format) to be used:
  + You can provide download URLs for these files as follows:
    + `gtf_url: 'protocol://url/for/some_file.gtf.gz'`
    + `gtf_name: 'some_file.gtf'`
    + `fasta_url: 'protocol://url/for/some_file.fasta.gz'`
    + `fasta_name: 'some_file.fasta'`
  + You can also place these files in the [references](./references) folder and just set the `gtf_name` and `fasta_name` fields. (Use '' for `gtf_url` and `fasta_url`)
* For each input sample, create a corresponding config entry as follows:
  + For samples with raw ONT signal data (FAST5 format), please provide the Guppy config file and the directory of FAST5 files as follows:
    + `guppy_config: 'the_guppy.cfg'`
    + `fast5_dir: '/path/to/fast5/dir'`
  + For samples with sequencing reads data (FASTQ format), please provide the appropriate file path as follows:
    + `fastq: '/path/to/the.fastq'`
  + Fpr samples with read-to-genome alignment data (SAM/BAM format), please provide **one** of the following fields as appropriate:
    + `sam: '/path/to/the.sam'`
    + `bam: '/path/to/the.bam'`
* Lastly, the following config values can be set to `true` or `false` as appropriate:
  + `use_annotated_junctions_with_minimap2`: Uses splice junctions recorded in the user-provided GTF as input to `minimap2`
  + `keep_espresso_c_temp`: Keep temporary files generated by the `C` step of `ESPRESSO`
  + `output_compatible_isoforms`: Generate a file (named `samples_N2_R0_compatible_isoform.tsv`) that maps long read IDs to their compatible transcript isoforms
  + `enable_visualization`: Generates files for visualizing transcript isoforms discovered by ESPRESSO. This requires setting other config values under "Visualization options" (**Note:** We recommend setting `enable_visualization` to false as we have separate scripts in this repository dedicated to generating visualizations for discovered transcript isoforms).

Edit the files in the folder [snakemake_profile](./snakemake_profile) to configure how jobs should be run in a cluster environment as follows:
* [config.yaml](./snakemake_profile/config.yaml): Sets various Snakemake parameters, such as whether jobs should be submitted to a cluster.
* [cluster_submit.py](./snakemake_profile/cluster_submit.py): Script to submit jobs.
* [cluster_status.py](./snakemake_profile/cluster_status.py): Script to check job status.
* [cluster_commands.py](./snakemake_profile/cluster_commands.py): Script to run commands specific to the cluster management system being used. The default implementation is for Slurm, but other cluster environments can be used by changing this file. For example, [cluster_commands_sge.py](./snakemake_profile/cluster_commands_sge.py) can be used if working with an SGE cluster.

For reference, examples of pre-configured `snakemake_config.yaml` files can be found in the folder [Config_for_each_system](./Config_for_each_system). Once the `snakemake_config.yaml` file has been appropriately configured, the Snakemake workflow can be run with `./run`

Upon completion, the following files will be generated in the folder `espresso_out/work_dir`:
* `samples_N2_R0_abundance.esp`: a tab-delimited file describing the expression levels of discovered transcript isoforms across input samples
  + Each discovered transcript isoform is reported on a separate line.
  + The first three columns are: `transcript_ID`, `transcript_name`, and `gene_ID`. Additional columns correspond to input samples and show the number of reads from a given sample that were counted towards each isoform.
  + Isoform read counts are assigned by expectation maximization, such that each read contributes at most 1 count, either to a single isoform or distributed as fractional counts to multiple isoforms.
* `samples_N2_R0_updated.gtf`: a transcript annotation file (GTF format) describing the coordinates of all discovered transcript isoforms
  + The `source` column indicates whether each transcript is a `novel_isoform` or an `annotated_isoform`
* `samples_N2_R0_compatible_isoform.tsv`: an optional tab-delimited file that describes compatible isoforms for each read in a given sample
  + The columns are `read_id`, `sample_name`, `read_classification`, and `compatible_isoforms`. Possible classifications for each read under the column `read_classification` include:
    + **Full Splice Match (FSM)**: indicates that the read carries a combination of splice junctions that is consistent with that of a known transcript
    + **Incomplete Splice Match (ISM)**: indicates that the read carries a combination of splice junctions that is part of a known transcript
    + **Novel In Catalog/Novel Not in Catalog (NIC/NNC)**: indicates that the read carries at least one splice junction involving a novel combination of known splice sites or novel splice sites, respectively.
    + **Not Completely Determined (NCD)**: indicates the read carries at least one splice junction that could not be corrected by ESPRESSO
    + **Single-exon**: indicates that the read does not carry any splice junctions

The Snakemake workflow will also produce log files that are named after the rules contained in the [Snakefile](./Snakefile). Specifically, there will be `{rule_name}_log.out` and `{rule_name}_log.err` files containing the stdout and stderr, respectively, of the command run for that rule. There will also be `.cluster.out`, `.cluster.err`, and `.cluster.usage` files if a rule was submitted to the cluster using [cluster_submit.py](./snakemake_profile/cluster_submit.py).

### Data visualization and analysis

After data processing, we can use the following scripts to further visualize and characterize the transcript isoforms discovered from our samples. 

#### Transcript isoform visualization

To visualize the isoform among given samples, the GENCODE annotation file is required, which should be downloaded and moved to [files](./scripts/Translation_scripts/files/) folder.

If GRCH37/hg19 is the expected genome version, we recommend this [GENCODE annotation file](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh37_mapping/gencode.v34lift37.annotation.gtf.gz). And if GRCH38/hg38 is the expected genome version, we recommend this [GENCODE annotation file](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.annotation.gtf.gz).

Our script can be run as follows (take DNMT3B and TP53 as examples): 
1. Enter [scripts](./scripts/) folder.
2. Run [Examples_visualization.sh](./scripts/Examples_visualization.sh).
3. Check the result figures generated in [Example_res](./scripts/Example_res/).


#### Detection of subtype-specific transcript isoforms

We have developed a script, [SubtypeSpecificIsoforms.py](./scripts/SubtypeSpecificIsoforms.py) that identify transcript isoforms with enriched usage in samples in given subtype relative to the samples from other subtypes within tumor. Briefly, our script uses the following approach to detect transcript isoforms with subtype-specific usage:
1. For each transcript identified by ESPRESSO, run a two-sided student-t test (FDR < user-specified threshold & Isoform proportion change (delta_proportion) >= user-specified threshold) between each subtype group and the corresponding background group (samples from the rest of subtypes) to identify isoform that has a significantly higher proportion in the given subtype group.

Our script can be run as follows:

```
python SubtypeSpecificIsoforms.py [-h] -i /path/to/isoform/proportion/matrix -a /path/to/sample-subtype_match_table -t <number of worker threads> \
    -c <FDR threshold> -d <delta_proportion threshold> -o /path/to/output/file

script arguments:
    -h, --help                                          show this message and exit

    -i /path/to/isoform_proportion_matrix               processed isoform proportion matrix

    -t <number of worker threads>                       number of worker threads (integer)

    -c <FDR threshold>                                  FDR threshold for calling transcript isoforms as subtype-specific 
                                                        (float between 0 and 1)

    -d <delta_proportion threshold>                     isoform proportion change threshold (float between 0 and 100(%))

    -a /path/to/sample-subtype_match_table              sample-subtype match table

    -o /path/to/output/file                             path to tsv file generated by script summarizing all transcript
                                                        isoforms prioritized as being subtype-specific
```

Our script will generate a tsv file summarizing all transcript isoforms prioritized as being subtype-specific at the specified FDR threshold and isoform proportion change threshold. Each line of the tsv file has the following output fields:

* **Field 1**: Gene symbol corresponding to gene of a subtype-specific isoform
* **Field 2**: Transcript ID of the subtype-specific isoform
* **Field 3**: Gene ID corresponding to gene of a subtype-specific isoform
* **Field 4**: Subtype pairs that have been tested
* **Field 5**: Average isoform proportion in correspoding subtype group
* **Field 6**: Raw p-value for subtype-specificity
* **Field 7**: FDR-adjusted p-value for subtype-specificity
* **Field 8**: Status of whether it is identified as significant


#### Detection of sample-specific transcript isoforms

We have developed a script, [SampleSpecificIsoforms.py](./scripts/SampleSpecificIsoforms.py) that identify transcript isoforms with enriched usage in a single sample relative to a collection of samples processed by ESPRESSO. Briefly, our script uses the following approach to detect transcript isoforms with sample-specific usage:
1. For each gene, generate an m-by-n contigency table comprised of read counts (rounded to the nearest integer) for m identified isoforms across n samples:
2. Run a chi-square test of homogeneity (FDR < user-specified threshold) on the contingency table to identify genes in which isoform proportions are not homogeneous across the considered samples.
3. For genes prioritized by the chi-square test with FDR < user-specified threshold, run a post-hoc one-tailed binomial test (FDR < user-specified threshold) to identify sample-isoform pairs in which the isoform proportion in the given sample is significantly higher than the global isoform proportion taken across all samples.

Our script can be run as follows:

```
python SampleSpecificIsoforms.py [-h] -i /path/to/ESPRESSO/isoform/abundance/matrix -t <number of worker threads> \
    -c <FDR threshold> -o /path/to/output/file

script arguments:
    -h, --help                                          show this message and exit

    -i /path/to/ESPRESSO/isoform/abundance/matrix       path to tsv file generated by ESPRESSO summarizing all detected transcript 
                                                        isoforms across samples and their assigned read counts

    -t <number of worker threads>                       number of worker threads (integer)

    -c <FDR threshold>                                  FDR threshold for calling transcript isoforms as sample-specific 
                                                        (float between 0 and 1)

    -o /path/to/output/file                             path to tsv file generated by script summarizing all transcript isoforms 
                                                        prioritized as being sample-specific
```

Our script will generate a tsv file summarizing all transcript isoforms prioritized as being sample-specific at the specified FDR threshold. Each line of the tsv file has the following output fields:

* **Field 1**: Gene ID corresponding to gene of a sample-specific isoform
* **Field 2**: Transcript ID of the sample-specific isoform
* **Field 3**: Name of the corresponding sample
* **Field 4**: Read counts assigned by ESPRESSO to sample-specific isoform in given sample
* **Field 5**: Global isoform proportion (sum of read counts of the isoform over all samples divided by the sum of read counts for all isoforms across all samples)
* **Field 6**: Isoform proportion in the given sample
* **Field 7**: Raw p-value for sample-specificity
* **Field 8**: FDR-adjusted p-value for sample-specificity

#### Characterization of alternative splicing events underlying discovered transcript isoforms

We wrote a script, [FindAltTSEvents.py](./scripts/FindAltTSEvents.py), that can enumerate all transcript structure differences between any given pair of transcript isoforms. These differences in transcript structure can then be classified into the following seven simple categories:

<img src="./files/AS_Patterns_Basic_2018_AJHG.jpg" width="800"/>

* Exon skipping (SE)
* Alternative 5'-splice site usage (A5SS)
* Alternative 3'-splice site usage (A3SS)
* Mutually exclusive exons (MXE)
* Intron retention (RI)
* Alternative first exon usage (AFE)
* Alternative last exon usage (ALE)

Any alternative transcript structure event that fails to fall into any of the seven categories listed above will, by default, be classified as complex (COMPLEX). **Note:** It is possible to have combinations of alternative transcript structure events for any given pair of transcript isoforms.

Our script can be run as follows:

```
python FindAltTSEvents.py [-h] -i /path/to/input/GTF -o /path/to/output/file

script arguments:
    -h, --help                                          show this message and exit

    -i /path/to/input/GTF                               path to GTF file describing structures of two transcript isoforms

    -o /path/to/output/file                             path to output file
```

Our script will subsequently generate a tab-delimited file consisting of four fields:
* **Field 1**: ID for transcript isoform 1
* **Field 2**: ID for transcript isoform 2
* **Field 3**: Identified alternative transcript structure event
* **Field 4**: Genomic coordinates for alternative transcript structure event.

**Note:** The designation for transcript isoforms 1 and 2 is completely arbitrary. Moreover, if the two transcript isoforms contained in the input GTF file exhibit a combination of multiple alternative transcript structure events, each event will be reported as its own line in the output file.

#### Prediction of NMD-targeted transcript isoforms
  
To determine NMD-targeted transcript, we firstly obtain open reading frame (ORF) for each identified transcript, either by adopting the annotated ORF for those annotated protein-coding transcripts with 'basic' tag (based on GENCODE annotation) or adopting its longest ORF as the predicted ORF for others.
For each transcript with predicted ORF, we determined it to contain a premature stop codon (PTC) and as the NMD-target transcript if: (i) it was longer than >= 200nts, (ii) it contained at least two exons and (iii) the last stop codon resided >=50nts upstream of the last exon-exon junction [Lindeboom, 2017, PMC5045715]. 

To run [Translation.py](./scripts/Translation_scripts/Translation.py), the genome sequence file and the corresponding GENCODE annotation file should be downloaded and moved to [files](./scripts/Translation_scripts/files/) folder.

If GRCH37/hg19 is the expected genome version, we recommend this [Genome fasta file](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh37_mapping/GRCh37.primary_assembly.genome.fa.gz) and this [GENCODE annotation file](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh37_mapping/gencode.v34lift37.annotation.gtf.gz).

If GRCH38/hg38 is the expected genome version, we recommend this [Genome fasta file](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/GRCh38.p13.genome.fa.gz) and this [GENCODE annotation file](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.annotation.gtf.gz).


Our script can be run as follows:

```
python Translation.py [-h] --mode long-read --trans_gtf /path/to/ESPRESSO/isoform/generated/gtf \
 --abundance_inf /path/to/ESPRESSO/isoform/abundance/matrix --genome_version GRCH37/GRCH38 \
 --ref_gtf /path/to/corresponding/GENCODE/annotation --out_dir /path/to/output/directory --out_file /prefix/name/of/output/file

script arguments:
    -h, --help                                          show this message and exit

    --mode                                              long-read/short-read data, please input 'long-read' here

    --abundance_inf                                     path to tsv file generated by ESPRESSO summarizing all detected transcript 
                                                        isoforms across samples and their assigned read counts

    --trans_gtf                                         path to gtf file generated by ESPRESSO summarizing all detected transcript isoforms                                                 

    --genome_version                                    genome sequence version, GRCH37/GRCH38

    --ref_gtf                                           corresponding GENCODE annotation file

    --out_dir                                           path to output directory

    --out_file                                          prefix of the output file  
```

Our script will generate protein fasta files for both protein-coding transcripts and NMD-targeted transcripts as the results.