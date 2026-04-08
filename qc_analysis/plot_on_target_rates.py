#!/usr/bin/env python3

import argparse
import os
import concurrent.futures
import pandas as pd
import pysam
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42

def parse_args():
    parser = argparse.ArgumentParser(description='Get mapping and on-target rates for each sample')
    parser.add_argument(
        "--mapping-file",
        required=True,
        help=("TSV file without a header and with the columns:"
              " name, bam, bed, group. The group column is optional"))
    parser.add_argument(
        "--outprefix",
        required=True,
        help="Prefix for output files .tsv and _mapping.pdf and _ontarget.pdf")
    parser.add_argument('--title')
    parser.add_argument("--threads", type=int, default=1, help="Number of parallel threads")
    return parser.parse_args()

def count_reads(bam, bed):
    ids=set()
    with pysam.AlignmentFile(bam,"rb") as f:
        with open(bed) as b:
            for line in b:
                if line.strip() and not line.startswith('#'):
                    c,s,e=line.split()[:3]
                    for r in f.fetch(c,int(s),int(e)):
                        if not r.is_unmapped and not r.is_secondary:
                            ids.add(r.query_name)
    return len(ids)

def count_mapped(bam):
    m,u=set(),set()
    with pysam.AlignmentFile(bam,"rb") as f:
        for r in f.fetch(until_eof=True):
            (u if r.is_unmapped else m).add(r.query_name)
    return len(m),len(u)

def get_counts(sample,group,bam,bed):
    m,u=count_mapped(bam)
    t=count_reads(bam,bed)
    return sample,group,m+u,m,u,t

def fmt(n):
    return f'{n/1e6:.1f}M' if n>=1e6 else f'{n//1e3}k' if n>=1e3 else str(n)

def rate_label(r):
    return '100%' if r==100 else '>99%' if r>=99 else f"{r:.0f}%" if r>=10 else f"{r:.1f}%" if r>=0.1 else '<0.1%'

def main():
    args=parse_args()

    df=pd.read_csv(args.mapping_file,sep="\t",header=None)
    while df.shape[1]<4: df[df.shape[1]]=None
    df.columns=["sample","bam","bed","group"]
    df["_order"]=range(len(df))

    missing=[f for f in df['bam'] if not os.path.exists(f)]
    if missing:
        raise FileNotFoundError("Missing BAM files:\n"+"\n".join(missing))

    df=df.assign(_s=df['bam'].map(os.path.getsize)).sort_values('_s',ascending=False).drop(columns='_s')

    print(f"Processing {len(df)} samples using {args.threads} threads...")

    results = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.threads) as ex:
        futures = {ex.submit(get_counts, r['sample'], r['group'], r['bam'], r['bed']): idx
                for idx, r in df.iterrows()}

        for f in concurrent.futures.as_completed(futures):
            res = f.result()
            idx = futures[f]
            results.append((idx, *res))  # prepend index
            sample, group, total, mapped, unmapped, target = res
            print(f"Finished: {sample} | total={total}, mapped={mapped}, target={target}")

    results.sort(key=lambda x: x[0])
    df = pd.DataFrame([r[1:] for r in results], columns=['sample','group','total','mapped','unmapped','target'])

    df['mapping_rate']=df['mapped']/df['total']*100
    df['on_target_rate']=df['target']/df['mapped']*100

    out_tsv=f"{args.outprefix}.tsv"
    df.to_csv(out_tsv,sep="\t",index=False)
    print(f"Saved TSV: {out_tsv}")

    groups=list(df['group'].dropna().unique())
    colors = [
        '#6997B9',  # blue
        '#BB6A68',  # red
        '#70A677',  # green
        '#D48653',  # orange
        '#A783A3',  # purple
    ]
    color_dict={g:colors[i%len(colors)] for i,g in enumerate(groups)} if groups else {}
    bar_colors=df['group'].map(color_dict).fillna('#6997B9')  # blue

    width=min(18,max(len(df)/2,6))

    # -------- FIGURE 1 --------
    fig,ax=plt.subplots(2,figsize=(width,10))

    ax[0].bar(range(len(df)),df['mapping_rate'],color="#FFD676")  # yellow
    ax[0].bar(range(len(df)),100-df['mapping_rate'],bottom=df['mapping_rate'],color='#C4C4C4')  # grey
    for i,r in df.iterrows():
        ax[0].text(i,r['mapping_rate']+1,rate_label(r['mapping_rate']),ha='center',fontsize=8)
    ax[0].set_ylim(0,105)
    ax[0].set_yticks([0, 20, 40, 60, 80, 100])
    ax[0].set_ylabel("Mapping Rate (%)")
    ax[0].set_xticks(range(len(df)))
    ax[0].set_xticklabels([])

    ax[1].bar(range(len(df)), df['mapped']/1e6, color="#FFD676")  # yellow
    ax[1].bar(range(len(df)), df['unmapped']/1e6, bottom=df['mapped']/1e6, color='#C4C4C4')  # grey
    offset=df['total'].max()*0.02
    for i,r in df.iterrows():
        ax[1].text(i,(r['total']+offset)/1e6,f"{fmt(r['total'])}",ha='center',fontsize=8)
    ax[1].set_ylim(0,df['total'].max()*1.05/1e6)
    ax[1].set_ylabel("Read Count (M)")
    ax[1].set_xticks(range(len(df)))
    ax[1].set_xticklabels(df['sample'],rotation=45,ha='right')

    if args.title: fig.suptitle(args.title)
    plt.tight_layout()
    out1=f"{args.outprefix}_mapping.pdf"
    plt.savefig(out1,bbox_inches='tight')
    print(f"Saved figure: {out1}")

    # -------- FIGURE 2 --------
    fig,ax=plt.subplots(2,figsize=(width,10))

    ax[0].bar(range(len(df)),df['on_target_rate'],color=bar_colors)
    ax[0].bar(range(len(df)),100-df['on_target_rate'],bottom=df['on_target_rate'],color='#C4C4C4')  # grey
    for i,r in df.iterrows():
        ax[0].text(i,r['on_target_rate']+1,rate_label(r['on_target_rate']),ha='center',fontsize=8)
    ax[0].set_ylim(0,105)
    ax[0].set_yticks([0, 20, 40, 60, 80, 100])
    ax[0].set_ylabel("On-target Rate (%)")
    ax[0].set_xticks(range(len(df)))
    ax[0].set_xticklabels([])

    off=df['mapped']-df['target']
    ax[1].bar(range(len(df)), df['target']/1e6, color=bar_colors)
    ax[1].bar(range(len(df)), off/1e6, bottom=df['target']/1e6, color='#C4C4C4')  # grey
    offset=df['mapped'].max()*0.02
    for i,r in df.iterrows():
        ax[1].text(i,(r['mapped']+offset)/1e6,f"{fmt(r['mapped'])}",ha='center',fontsize=8)
    ax[1].set_ylim(0,df['mapped'].max()*1.05/1e6)
    ax[1].set_ylabel("Read Count (M)")
    ax[1].set_xticks(range(len(df)))
    ax[1].set_xticklabels(df['sample'],rotation=45,ha='right')

    if args.title: fig.suptitle(args.title)
    plt.tight_layout()
    out2=f"{args.outprefix}_ontarget.pdf"
    plt.savefig(out2,bbox_inches='tight')
    print(f"Saved figure: {out2}")

if __name__=="__main__":
    main()
