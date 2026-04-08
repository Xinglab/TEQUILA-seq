#!/usr/bin/env python3

import argparse
import gzip
from concurrent.futures import ProcessPoolExecutor, as_completed
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import pandas as pd
import numpy as np
import pysam
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42

def parse_args():
    """ Parse command-line arguments for read length and quality plotting """

    parser = argparse.ArgumentParser(description="Make a plot of read lengths and qualities")
    parser.add_argument(
        "--mapping-file",
        required=True,
        help=("TSV file without a header and with the columns:"
              " name, fastq, bam, bed. The bed column is optional and"
              " defines the target regions for that row"))
    parser.add_argument(
        "--outprefix",
        required=True,
        help=("Prefix for output files _read_attributes.tsv,"
              " _summary.tsv, _read_lengths_boxplot.pdf, _read_qualities_boxplot.pdf"))
    parser.add_argument('--title')
    parser.add_argument("--threads", type=int, default=1, help="Number of parallel threads")
    return parser.parse_args()

def get_mapped_readIDs(bam):
    """Get mapped read IDs from a BAM file."""

    mapped_readIDs = set()
    unmapped_readIDs = set()

    with pysam.AlignmentFile(bam, "rb") as bamfile:
        for read in bamfile.fetch(until_eof=True):
            if read.is_unmapped:
                unmapped_readIDs.add(read.query_name)
            else:
                mapped_readIDs.add(read.query_name)

    return mapped_readIDs, unmapped_readIDs

def get_target_type(bam, bed):
    """Get on-target and off-target read IDs from a BAM file based on regions defined in a BED file."""

    mapped_readIDs, unmapped_readIDs = get_mapped_readIDs(bam)

    on_target_readIDs = set()
    off_target_readIDs = set()
    with pysam.AlignmentFile(bam, "rb") as bamfile:
        with open(bed) as bedfile:
            for line in bedfile:
                if line.strip():
                    chrom, start, end = line.split()[:3]
                    for read in bamfile.fetch(chrom, int(start), int(end)):
                        if not (read.is_unmapped or read.is_secondary):
                            on_target_readIDs.add(read.query_name)

    off_target_readIDs = mapped_readIDs - on_target_readIDs

    return on_target_readIDs, off_target_readIDs, unmapped_readIDs

def get_read_info(name, fastq_file, on_target_readIDs, off_target_readIDs, mapped_readIDs, unmapped_readIDs):
    """Extract read length, mean Phred quality, and target type from FASTQ."""

    records = []
    open_func = gzip.open if fastq_file.endswith(".gz") else open

    with open_func(fastq_file, "rt") as handle:
        for read_id, seq, qual in FastqGeneralIterator(handle):
            rid = read_id.split()[0]
            seq_len = len(seq)

            if seq_len > 0:
                phred_scores = np.frombuffer(qual.encode(), dtype=np.uint8) - 33
                mean_phred = phred_scores.mean()
            else:
                mean_phred = 0

            if rid in on_target_readIDs:
                target_type = "on_target"
            elif rid in off_target_readIDs:
                target_type = "off_target"
            elif rid in mapped_readIDs:
                target_type = "mapped"
            elif rid in unmapped_readIDs:
                target_type = "unmapped"
            else:
                target_type = "unknown"

            records.append({
                "Sample": name,
                "ReadLength": seq_len,
                "Phred": mean_phred,
                "TargetType": target_type
            })

    return records

def process_sample(name, fastq, bam, bed):
    """Process one sample; returns list of dicts"""

    if bed is None:
        mapped_readIDs, unmapped_readIDs = get_mapped_readIDs(bam)
        on_target_readIDs = set()
        off_target_readIDs = set()
    else:
        on_target_readIDs, off_target_readIDs, unmapped_readIDs = get_target_type(bam, bed)
        mapped_readIDs = set()

    records = get_read_info(name, fastq, on_target_readIDs, off_target_readIDs, mapped_readIDs, unmapped_readIDs)

    return records

def main():
    """ Main function """

    args = parse_args()

    df_mapping = pd.read_csv(args.mapping_file, sep="\t", header=None)
    if df_mapping.shape[1] == 3:
        df_mapping.columns = ["name", "fastq", "bam"]
        df_mapping["bed"] = None
    else:
        df_mapping.columns = ["name", "fastq", "bam", "bed"]
    sample_order = df_mapping["name"].tolist()
    samples = df_mapping.to_dict(orient="records")

    print(f"Processing {len(df_mapping)} samples using {args.threads} threads...")

    dfs = []

    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        futures = {executor.submit(process_sample, s["name"], s["fastq"], s["bam"], s["bed"]): s["name"] for s in samples}
        for future in as_completed(futures):
            name = futures[future]
            try:
                result = future.result()
                dfs.append(pd.DataFrame(result))
                print(f"Finished {name}")
            except Exception as e:
                print(f"❌ Error in {name}: {e}")

    print("Merging attribute info...")
    df_all = pd.concat(dfs, ignore_index=True)

    df_all = df_all.sort_values('Sample', key=lambda x: x.map({s: i for i, s in enumerate(sample_order)}))
    df_all.to_csv(f"{args.outprefix}_read_attributes.tsv", sep="\t", index=False)
    print(f"\n✅ Attributes saved to {args.outprefix}_read_attributes.tsv")

    print("\nGenerating summary statistics...")

    summary_overall = (
        df_all
        .groupby("Sample")
        .agg(
            MeanReadLength=("ReadLength", "mean"),
            MedianReadLength=("ReadLength", "median"),
            MeanPhred=("Phred", "mean"),
            MedianPhred=("Phred", "median"),
            TotalReads=("ReadLength", "count")
        )
        .reset_index()
    )
    summary_overall["TargetType"] = "ALL"

    summary_by_type = (
        df_all
        .groupby(["Sample", "TargetType"])
        .agg(
            MeanReadLength=("ReadLength", "mean"),
            MedianReadLength=("ReadLength", "median"),
            MeanPhred=("Phred", "mean"),
            MedianPhred=("Phred", "median"),
            TotalReads=("ReadLength", "count")
        )
        .reset_index()
    )

    summary = pd.concat([summary_overall, summary_by_type], ignore_index=True)

    default_target_order = ["on_target", "off_target", "mapped", "unmapped", "unknown"]
    target_order = default_target_order + [
        t for t in summary['TargetType'].drop_duplicates()
        if t not in default_target_order
    ]
    sample_order = [s for s in sample_order if s in summary['Sample'].unique()]
    target_order = [t for t in target_order if t != "ALL" and t in summary['TargetType'].unique()]
    summary['Sample_order'] = summary['Sample'].map({s: i for i, s in enumerate(sample_order)})
    summary['Target_order'] = summary['TargetType'].map({t: i for i, t in enumerate(target_order)})
    summary_sorted = summary.sort_values(['Sample_order', 'Target_order']).drop(columns=['Sample_order', 'Target_order'])

    summary_file = f"{args.outprefix}_summary.tsv"
    summary_sorted.to_csv(summary_file, sep="\t", index=False)

    print(f"✅ Summary statistics saved to {summary_file}")

    df_ds = (
        df_all
        .groupby(["Sample", "TargetType"], group_keys=True)
        .apply(lambda d: d.sample(min(len(d), 100_000), random_state=0))
    )

    palette_dict = dict(zip(target_order, sns.color_palette("pastel")))

    print("\nPlotting read length boxplot...")
    fig, ax = plt.subplots(figsize=(max(6, df_all['Sample'].nunique() * 0.5), 6))
    sns.boxplot(data=df_ds, x="Sample", y="ReadLength", hue="TargetType",
                order=sample_order, hue_order=target_order,
                palette=palette_dict, showfliers=False, width=0.7, ax=ax)
    ymax = ax.get_ylim()[1]
    plt.xticks(rotation=45, ha="right")
    plt.ylim(0, ymax * 1.05)
    ax.margins(y=0.05)
    if args.title: fig.suptitle(args.title)
    fig.tight_layout()
    fig.savefig(f"{args.outprefix}_read_lengths_boxplot.pdf")
    print(f"✅ Read length box plot saved to {args.outprefix}_read_lengths_boxplot.pdf")

    print("\nPlotting read quality boxplot...")
    fig, ax = plt.subplots(figsize=(max(6, df_all['Sample'].nunique() * 0.5), 6))
    sns.boxplot(data=df_ds, x="Sample", y="Phred", hue="TargetType",
                order=sample_order, hue_order=target_order,
                palette=palette_dict, showfliers=False, width=0.7, ax=ax)
    ymax = ax.get_ylim()[1]
    plt.xticks(rotation=45, ha="right")
    plt.ylim(0, ymax * 1.05)
    ax.margins(y=0.05)
    if args.title: fig.suptitle(args.title)
    fig.tight_layout()
    fig.savefig(f"{args.outprefix}_read_qualities_boxplot.pdf")
    print(f"✅ Read quality box plot saved to {args.outprefix}_read_qualities_boxplot.pdf")

if __name__ == "__main__":
    main()
