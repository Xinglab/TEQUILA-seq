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
	parser.add('-ac', '--abundance_CPM', dest='abundance_CPM', type=str, help='e.g. samples_N2_R0_abundance_combined_with_tissue_CPM_more_than_1_change_name_reshaped_merge_others.txt', required=True)
	parser.add('-ap', '--abundance_proportion', dest='abundance_proportion', type=str, help='e.g. samples_N2_R0_abundance_combined_with_tissue_CPM_proportion_more_than_1_change_name_reshaped_merge_others.txt', required=True)
	parser.add('-b', '--bedgraph', dest='bedgraph', type=str, help='sample BedGraph file', required=True)
	parser.add('-of', '--order', dest='order', type=str, default='no', help='order samples based on proportion?')
	parser.add('-od', '--out_dir', dest='out_dir', type=str, help='output directory', required=True)
	parser.add('-sg', '--sorted_group', dest='sorted_group', type=str, help='sorted group string, separated by \',\'', required=True)
	parser.add('-v', '--genome_version', dest='genome_version', type=str, default='hg38', help='genome version: hg19 or hg38')
	parser.add('-ct', '--canonical_transcript', dest='canonical_transcript', type=str, help='canonical transcript')
	parser.add('-bt', '--basic_transcript', dest='basic_transcript', type=str, help='basic transcript')
	parser.add('-ag', '--anno_gtf', dest='anno_gtf', type=str, help='annotated gencode gtf')
	parser.add('-fi', '--figures', dest='figures', nargs='+', help='Choose from [Isoform, Single_isoform, Structure]', default = 'Isoform Single_isoform Structure')
	args = parser.parse_args()
	return parser, args

############ main script ##############
parser, args = parse_args()

gene = args.gene
gene_name = args.gene_name
transcript = args.transcript
abundance_CPM_original = args.abundance_CPM_original
abundance_CPM = args.abundance_CPM
abundance_proportion = args.abundance_proportion
bedgraph = args.bedgraph
out_dir = args.out_dir
order = args.order
sorted_group = args.sorted_group
genome_version = args.genome_version
canonical_transcript = args.canonical_transcript
basic_transcript = args.basic_transcript
anno_gtf = args.anno_gtf
figures = args.figures
print (figures)

with open(abundance_CPM_original,'r') as inf:
	for index, line in enumerate(inf):
		arr = line.strip().split("\t")
		if index == 0:
			name_col = len(arr[3].split('_'))+1

if sorted_group == 'Melanoma':
	sorted_group = "\"Melanocytic,Transitory,Neural crest like,Undifferentiated\""
	

file_dir = os.path.dirname(os.path.realpath(__file__))

#### 1. Generate proportion bar plot ####
if 'Isoform' in figures:
	command_1 = "Rscript %s/2_2_1_Bar_sample_isoform_single_color.R %s %s %s %s %s %s %s %s %s %s" % (file_dir, gene, gene_name, transcript, out_dir, abundance_proportion, abundance_CPM, order, sorted_group, canonical_transcript, basic_transcript)
	print (command_1)
	os.system(command_1)
	print ('Step 1 completed: Generate proportion bar plot\n')


#### 1_5. Generate proportion bar plot ####
if 'Single_isoform' in figures:
	command_1_5 = "Rscript %s/2_2_2_Bar_sample_isoform_single_color_one_isoform.R %s %s %s %s %s %s %s %s %s %s" % (file_dir, gene, gene_name, transcript, out_dir, abundance_proportion, abundance_CPM, order, sorted_group, canonical_transcript, basic_transcript)
	print (command_1_5)
	os.system(command_1_5)
	print ('Step 1.5 completed: Generate proportion bar plot (single transcript)\n')


#### 2. Generate transcript structure plot ####
if 'Structure' in figures:
	command_2 = "python %s/2_2_3_Isoform_structure_generate_single_color.py %s %s %s %s %s %s %s %s" % (file_dir, gene_name, transcript, out_dir, name_col, abundance_CPM_original, bedgraph, genome_version, anno_gtf)
	print (command_2)
	os.system(command_2)
	print ('Step 2 completed: Generate transcript structure plot.\n')


