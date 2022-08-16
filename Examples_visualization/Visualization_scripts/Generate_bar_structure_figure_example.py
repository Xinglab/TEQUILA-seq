####################################################
#Goal:      Generate bar and structure figures for examples
#Author:    Yang Xu
#E-mail:    yangax@pennmedicine.upenn.edu
####################################################

import os,sys,re
from collections import defaultdict
import configargparse

def parse_args():
	parser = configargparse.ArgParser(description='Generate bar and structure figures for examples')
	parser.add('-g', '--gene', dest='gene', type=str, help='Target gene ID', required=True)
	parser.add('-gn', '--gene_name', dest='gene_name', type=str, help='Target gene symbol', required=True)
	parser.add('-t', '--transcript', dest='transcript', type=str, help='Target transcript', required=True)
	parser.add('-a', '--abundance_CPM_original', dest='abundance_CPM_original', type=str, help='e.g. ../../samples_N2_R0_abundance_combined_with_tissue_CPM_more_than_1.esp', required=True)
	parser.add('-c', '--name_col', dest='name_col', type=int, help='the index of transcript ID column in the filename, 0-based', required=True)
	parser.add('-ac', '--abundance_CPM', dest='abundance_CPM', type=str, help='e.g. samples_N2_R0_abundance_combined_with_tissue_CPM_more_than_1_change_name_reshaped_merge_others.txt', required=True)
	parser.add('-ap', '--abundance_proportion', dest='abundance_proportion', type=str, help='e.g. samples_N2_R0_abundance_combined_with_tissue_CPM_proportion_more_than_1_change_name_reshaped_merge_others.txt', required=True)
	parser.add('-b', '--bedgraph', dest='bedgraph', type=str, help='sample BedGraph file', required=True)
	parser.add('-of', '--order', dest='order', type=str, default='no', help='order samples based on proportion?')
	parser.add('-od', '--out_dir', dest='out_dir', type=str, help='output directory', required=True)
	args = parser.parse_args()
	return parser, args

############ main script ##############
parser, args = parse_args()

gene = args.gene
gene_name = args.gene_name
transcript = args.transcript
name_col = args.name_col
abundance_CPM_original = args.abundance_CPM_original
abundance_CPM = args.abundance_CPM
abundance_proportion = args.abundance_proportion
bedgraph = args.bedgraph
out_dir = args.out_dir
order = args.order

#### 1. Generate proportion bar plot ####
command_1 = "Rscript ./Visualization_scripts/Bar_sample_isoform_single_color.R %s %s %s %s %s %s %s" % (gene, gene_name, transcript, out_dir, abundance_proportion, abundance_CPM, order)
print (command_1)
os.system(command_1)
print ('Step 1 completed: Generate proportion bar plot\n')


#### 1_5. Generate proportion bar plot ####
command_1_5 = "Rscript ./Visualization_scripts/Bar_sample_isoform_single_color_one_isoform.R %s %s %s %s %s %s %s" % (gene, gene_name, transcript, out_dir, abundance_proportion, abundance_CPM, order)
print (command_1_5)
os.system(command_1_5)
print ('Step 1.5 completed: Generate proportion bar plot (single transcript)\n')


#### 2. Generate transcript structure plot ####
command_2 = "python ./Visualization_scripts/Isoform_structure_generate_single_color.py %s %s %s %s %s %s" % (gene_name, transcript, out_dir, name_col, abundance_CPM_original, bedgraph)
print (command_2)
os.system(command_2)
print ('Step 2 completed: Generate transcript structure plot.\n')



