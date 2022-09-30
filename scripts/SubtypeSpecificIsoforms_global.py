#!/usr/bin/env python3

'''
This is a script to identify subtype-specific isoforms based on
an isoform read count matrix generated by ESPRESSO (rows are 
detected isoforms and columns are samples)
'''

# Load required libraries
import argparse
import numpy as np
import pandas as pd
import concurrent.futures as cf
import scipy
import re,sys,os
from scipy.stats import binom, chi2, chi2_contingency
from statsmodels.stats.multitest import multipletests
from collections import defaultdict


ID2name_dict = defaultdict()
with open('../files/Target_IMPACT_gene_list_onco_suppre.txt','r') as gene_type_inf:
    for index,line in enumerate(gene_type_inf):
        if index == 0: continue
        arr = line.strip().split('\t')
        ID2name_dict[arr[1]] = arr[0]

def t_test(sublist, cutoff):
    [sub_1, sub_2] = sublist
    res2p_dict = defaultdict()
    for each_trans in isoform_proportion_dict:
        each_gene = trans2gene_dict[each_trans]
        sub_1_list = np.array(isoform_proportion_dict[each_trans][sub_1])
        if sub_2 in isoform_proportion_dict[each_trans]:
            sub_2_list = np.array(isoform_proportion_dict[each_trans][sub_2])
        elif sub_2 == 'global':
            sub_2_list = []
            for each_sub in isoform_proportion_dict[each_trans]:
                if each_sub != sub_1:
                    sub_2_list += isoform_proportion_dict[each_trans][each_sub]
            sub_2_list = np.array(sub_2_list)
        mean_sub_1 = np.nanmean(sub_1_list)
        mean_sub_2 = np.nanmean(sub_2_list)
        if list(map(str,set(sub_1_list)))==['nan'] or list(map(str,set(sub_2_list)))==['nan']: 
            continue
        #if float(mean_sub_1) > float(mean_sub_2):
        p_value = scipy.stats.ttest_ind(sub_1_list, sub_2_list, equal_var=False, nan_policy='omit', alternative='two-sided')[1]
        if str(p_value) == "nan": continue
        if str(p_value) == "--": continue
        res = ';'.join([each_trans, each_gene, sub_1+','+sub_2, str(mean_sub_1)+','+str(mean_sub_2), str(p_value)])
        res2p_dict[res] = p_value

    res_list = []
    sorted_res2p = sorted(res2p_dict.items(), key=lambda x:float(x[1]))
    for i in range(0, len(sorted_res2p)):
        rank = float(i + 1)
        p_adjust = float(sorted_res2p[i][1]) * len(sorted_res2p) / rank
        print (sorted_res2p[i][0], p_adjust, len(sorted_res2p), rank)
        #print (rank, sorted_res2p[i][0].split(';')[0], float(sorted_res2p[i][1]), p_adjust)
        #delta_mean = float(sorted_res2p[i][0].split(';')[3].split(',')[0]) - float(sorted_res2p[i][0].split(';')[3].split(',')[1])
        significance = 'Not_significant'
        if float(p_adjust) <= cutoff:
            significance = 'Significant'        
        res_list.append(sorted_res2p[i][0]+';'+str(p_adjust)+';'+significance)
    return res_list


def SubtypeSpecificIsoforms(infile, anno_table, threads, cutoff, outfile):
    # Load sample-subtype pair
    print('Load sample-subtype pair...', flush=True)
    code2subtype_dict = defaultdict()
    cell2subtype_dict = defaultdict()
    code2cell_dict = defaultdict()
    with open(anno_table, 'r') as anno_table_inf:
        for index, line in enumerate(anno_table_inf):
            if index == 0: continue
            arr = line.strip().split('\t')
            code2cell_dict[arr[0]] = arr[1]
            code2subtype_dict[arr[0]] = arr[2]
            cell2subtype_dict[arr[1]] = arr[2]

    # Parse isoform read count matrix and extract list of genes
    print('Parsing isoform read count matrix...', flush=True)
    value_col = 2
    subtype_list = []
    sample_list = []
    global isoform_proportion_dict
    isoform_proportion_dict = defaultdict()
    global trans2gene_dict
    trans2gene_dict = defaultdict()
    with open(infile, 'r') as inf:
        for index, line in enumerate(inf):
            arr = line.strip().split('\t')
            if index == 0:
                sample_list = arr[value_col:len(arr)]
                for i in range(value_col, len(arr)):
                    this_sample = '_'.join(arr[i].split('_')[0:-1])
                    this_subtype = cell2subtype_dict[this_sample]
                    if this_subtype not in subtype_list:
                        subtype_list.append(this_subtype)
            else:
                trans_ID = arr[0].split('.')[0]
                gene_ID = arr[1].split('.')[0]
                if trans_ID not in isoform_proportion_dict:
                    isoform_proportion_dict[trans_ID] = defaultdict()
                    trans2gene_dict[trans_ID] = gene_ID
                for i in range(value_col, len(arr)):
                    this_sample = '_'.join(sample_list[i-value_col].split('_')[0:-1])
                    this_subtype = cell2subtype_dict[this_sample]
                    if this_subtype not in isoform_proportion_dict[trans_ID]:
                        isoform_proportion_dict[trans_ID][this_subtype] = []
                    isoform_proportion_dict[trans_ID][this_subtype].append(float(arr[i]))
    outf= open(outfile,'w')
    outf.write('Gene_name\tTranscript\tGene\tSubtype_pair\tProportion\tP_value\tAdjust_P_value\tSignificance\n')
    outfile_significant = re.sub('.txt','_significant.txt', outfile)
    outf_significant = open(outfile_significant, 'w')
    outf_significant.write('Gene_name\tTranscript\tGene\tSubtype_pair\tProportion\tP_value\tAdjust_P_value\tSignificance\n')
    subtype_pair_list = []
    subtype_list_2 = subtype_list
    for i_1 in range(0, len(subtype_list)):
        subtype_1 = subtype_list[i_1]
        subtype_2 = 'global'
        subtype_pair_list.append([subtype_1, subtype_2])
        outf_list = t_test([subtype_1, subtype_2], cutoff)
        for each_outf_record in outf_list:
            Gene_name = ID2name_dict[each_outf_record.split(';')[1]]
            outf.write(Gene_name+'\t'+'\t'.join(each_outf_record.split(';'))+'\n')
            delta = float(each_outf_record.split(';')[3].split(',')[0]) - float(each_outf_record.split(';')[3].split(',')[1])
            if each_outf_record.split(';')[-1] == 'Significant' and delta >= 10:
                outf_significant.write(Gene_name+'\t'+'\t'.join(each_outf_record.split(';'))+'\n')
    print('Subtype pair:', subtype_pair_list)
    outf.close()
    outf_significant.close()

def main():
    moduleSummary = 'This is a script to call subtype-specific isoforms from an isoform read count matrix generated by ESPRESSO'
    parser = argparse.ArgumentParser(description=moduleSummary)

    # Add arguments
    parser.add_argument('-i', metavar='/path/to/read/count/matrix', required=True,
        help='path to isoform read count matrix generated by ESPRESSO')
    parser.add_argument('-a', metavar='###', required=True,
        help='path to sample-subtype match table')
    parser.add_argument('-t', metavar='###', required=True,
        help='number of worker threads')
    parser.add_argument('-c', metavar='###', required=True,
        help='FDR threshold (between 0 and 1)') 
    parser.add_argument('-o', metavar='/path/to/output/file', required=True,
        help='path to output file')
    
    # Parse command-line arguments
    args = parser.parse_args()
    infile, anno_table, threads, cutoff, outfile = args.i, args.a, int(args.t), float(args.c), args.o

    print('Isoform read count matrix: ' + infile, flush=True)
    print('Sample-subtype matched table: ' + anno_table, flush=True)
    print('Number of threads: ' + str(threads), flush=True)
    print('FDR cutoff: ' + str(cutoff), flush=True)
    print('Output file: ' + outfile, flush=True)

    # Run SubtypeSpecificIsoforms
    SubtypeSpecificIsoforms(infile, anno_table, threads, cutoff, outfile)

if __name__ == '__main__':
    main()
