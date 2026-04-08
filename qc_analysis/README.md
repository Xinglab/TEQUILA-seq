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

### plot_read_attributes.py

In order to run `plot_read_attributes.py`, create a file `./attributes_in.tsv` with no header and one row for each sample where the columns are:
* sample name
* path to a fastq file
* path to a bam file
* a bed file defining the target regions

```
sample_1	sample_1.fastq.gz	sample_1.bam	targets.bed
sample_2	sample_2.fastq.gz	sample_2.bam	targets.bed
```

Then run:
```
python plot_read_attributes.py --mapping-file ./attributes_in.tsv --outprefix ./attributes_out --title 'Test Run' --threads 2
```

The output files will be:
* `attributes_out_read_attributes.tsv`
* `attributes_out_read_lengths_boxplot.pdf`
* `attributes_out_read_qualities_boxplot.pdf`
* `attributes_out_summary.tsv`

Running with two samples each with about 3M reads (3.5GB bam files) using 2 threads took 6 minutes and used 3GB of memory

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

`./rates_in.tsv`
```
sample_1	sample_1.bam	targets.bed	group_1
sample_2	sample_2.bam	targets.bed	group_2
```

Then run:
```
python plot_on_target_rates.py --mapping-file ./rates_in.tsv --outprefix ./rates_out --title 'Test Run' --threads 2
```

The output files will be:
* `rates_out.tsv`
* `rates_out_mapping.pdf`
* `rates_out_ontarget.pdf`

Running with two samples each with about 3M reads (3.5GB bam files) using 2 threads took 2 minutes and used 0.6GB of memory

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
