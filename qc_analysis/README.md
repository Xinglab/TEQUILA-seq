# TEQUILA-seq Quality Control

## Dependencies

The dependencies can be installed to a conda environment by running:
```
conda create --prefix ./conda_env
conda activate ./conda_env
conda install -c conda-forge -c bioconda --file ./conda_requirements.txt
```

## Usage

Before running, make sure that all `.bam` files are sorted and indexed:
```
samtools sort -o sample_name_sorted.bam sample_name_unsorted.bam
samtools index sample_name_sorted.bam
```

Example input files are available at: https://xinglabtrackhub.research.chop.edu/tequila/qc_example/
* `SY5Y_TEQUILA.bam`
* `SY5Y_TEQUILA.bam.bai`
* `SY5Y_TEQUILA.fastq.gz`
* `SY5Y_WTS.bam`
* `SY5Y_WTS.bam.bai`
* `SY5Y_WTS.fastq.gz`
* `target.bed`

The files are based on data from:

DeBruyne N, Wang F, Xu Y, Lin L. Evaluating the potential and limitations of nanopore adaptive sampling for targeted transcriptome sequencing. Genome Biol. 2025 Oct 9;26(1):349. doi: 10.1186/s13059-025-03813-1. PMID: 41068925; PMCID: PMC12509409.

The `.bam` files were created using minimap2 v2.25. The minimap2 package includes paftools.js which can be used to create a `.bed` input annotation file. GENCODE v44 primary assembly files were used as input

Prepare the reference files:
```
curl -L 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.primary_assembly.annotation.gtf.gz' -O
curl -L 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz' -O
gunzip gencode.v44.primary_assembly.annotation.gtf.gz
gunzip GRCh38.primary_assembly.genome.fa.gz
paftools.js gff2bed gencode.v44.primary_assembly.annotation.gtf > gencode.v44.primary_assembly.annotation.bed
```

Align with minimap2 (about 6 minutes and 22GB of memory):
```
minimap2 -ax splice -ub --secondary=no -t 8 -k 14 -w 4 --junc-bed gencode.v44.primary_assembly.annotation.bed GRCh38.primary_assembly.genome.fa SY5Y_TEQUILA.fastq.gz > SY5Y_TEQUILA.sam
```

Convert the `.sam` to a sorted and indexed `.bam`:
```
samtools sort -o SY5Y_TEQUILA.bam SY5Y_TEQUILA.sam
samtools index SY5Y_TEQUILA.bam
```

### plot_read_attributes.py

In order to run `plot_read_attributes.py`, create a file `./attributes_in.tsv` with no header and one row for each sample where the columns are:
* sample name
* path to a fastq file
* path to a bam file
* a bed file defining the target regions

```
SY5Y_WTS	SY5Y_WTS.fastq.gz	SY5Y_WTS.bam	target.bed
SY5Y_TEQUILA	SY5Y_TEQUILA.fastq.gz	SY5Y_TEQUILA.bam	target.bed
```

Then run:
```
python plot_read_attributes.py --mapping-file ./attributes_in.tsv --outprefix ./attributes_out --title 'TEQUILA-seq: Test' --threads 2
```

The output files will be:
* [attributes_out_read_lengths_boxplot.pdf](./example_out/attributes_out_read_lengths_boxplot.pdf)
* [attributes_out_read_qualities_boxplot.pdf](./example_out/attributes_out_read_qualities_boxplot.pdf)
* `attributes_out_read_attributes.tsv`
* `attributes_out_summary.tsv`

A larger test dataset with two samples each with about 3M reads (3.5GB bam files) using 2 threads took 6 minutes and used 3GB of memory

```
python plot_read_attributes.py -h

usage: plot_read_attributes.py [-h] --mapping-file MAPPING_FILE
                               --outprefix OUTPREFIX [--title TITLE]
                               [--threads THREADS]

Make a plot of read lengths and qualities

options:
  -h, --help            show this help message and exit
  --mapping-file MAPPING_FILE
                        TSV file without a header and with the columns: name,
                        fastq, bam, bed. The bed column is optional and
                        defines the target regions for that row
  --outprefix OUTPREFIX
                        Prefix for output files _read_attributes.tsv,
                        _summary.tsv, _read_lengths_boxplot.pdf,
                        _read_qualities_boxplot.pdf
  --title TITLE
  --threads THREADS     Number of parallel threads
```

### plot_on_target_rates.py

In order to run `plot_on_target_rates.py`, create a file `./rates_in.tsv` with no header and one row for each sample where the columns are:
* sample name
* path to a bam file
* a bed file defining the target regions
* group name

```
SY5Y_WTS	SY5Y_WTS.bam	target.bed	SY5Y_WTS
SY5Y_TEQUILA	SY5Y_TEQUILA.bam	target.bed	SY5Y_TEQUILA
```

Then run:
```
python plot_on_target_rates.py --mapping-file ./rates_in.tsv --outprefix ./rates_out --title 'TEQUILA-seq: Test' --threads 2
```

The output files will be:
* [rates_out_mapping.pdf](./example_out/rates_out_mapping.pdf)
* [rates_out_ontarget.pdf](./example_out/rates_out_ontarget.pdf)
* `rates_out.tsv`

A larger test dataset with two samples each with about 3M reads (3.5GB bam files) using 2 threads took 2 minutes and used 0.6GB of memory

```
python plot_on_target_rates.py -h

usage: plot_on_target_rates.py [-h] --mapping-file MAPPING_FILE
                               --outprefix OUTPREFIX [--title TITLE]
                               [--threads THREADS]

Get mapping and on-target rates for each sample

options:
  -h, --help            show this help message and exit
  --mapping-file MAPPING_FILE
                        TSV file without a header and with the columns: name,
                        bam, bed, group. The group column is optional
  --outprefix OUTPREFIX
                        Prefix for output files .tsv and _mapping.pdf and
                        _ontarget.pdf
  --title TITLE
  --threads THREADS     Number of parallel threads
```
