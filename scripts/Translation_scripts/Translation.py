####################################################
#Goal:      Translate transcripts into proteins
#Author:    Yang Xu
#E-mail:    yangax@pennmedicine.upenn.edu
####################################################

import os,sys,re
from collections import defaultdict
import configargparse

def parse_args():
	parser = configargparse.ArgParser(description='Translate transcripts into proteins')
	parser.add('-m', '--mode', dest='mode', type=str, help='\"long-read\" or \"short-read\"', required=True)
	parser.add('-i', '--trans_gtf', dest='trans_gtf', type=str, help='input gtf file of identified transcripts', required=True)
	parser.add('-a', '--abundance_inf', dest='abundance_inf', type=str, help='input file of transcript abundance, e.g. samples_N2_R0_abundance.esp', required=True)
	parser.add('-od', '--out_dir', dest='out_dir', type=str, help='output directory', required=True)
	parser.add('-of', '--out_file', dest='out_file', type=str, help='output file name prefix', required=True)
	parser.add('-r', '--ref_gtf', dest='ref_gtf', type=str, help='reference annotated gtf file', required=True)
	#parser.add('-g', '--genome', dest='genome', type=str, help='human reference genome fasta', required=True)
	parser.add('-v', '--genome_version', dest='genome_version', type=str, help='genome_version, such as \"hg19\" or \"hg38\"', required=True)
	args = parser.parse_args()
	return parser, args

############ main script ##############
parser, args = parse_args()

dir_path = os.path.dirname(os.path.realpath(__file__))

#genome = args.genome
if args.genome_version in ['hg19','GRCH37','grch37']:
	genome = '%s/files/hg19.fa' % dir_path
elif args.genome_version in ['hg38','GRCH38','grch38']:
	genome = '%s/files/GRCh38.primary_assembly.genome.fa' % dir_path

input_trans_gtf = args.trans_gtf
ref_gtf = args.ref_gtf
out_dir = args.out_dir
out_file = args.out_file
abun_inf = args.abundance_inf
mode = args.mode

#### 0. process reference gtf file to generate CDS annotation ####
command_0 = "python %s/0_process_gtf.py %s %s" % (dir_path, ref_gtf, out_dir)
print (command_0)
os.system(command_0)
print ('Step 0 completed: reference gtf has been converted to CDS annotation file.\n')

#### 1. convert input gtf file into bed file ####
command_1 = "python %s/1_convert_gtf2transcript_1_to_1_based.py %s %s/1_%s_gtf_processed.txt" % (dir_path, input_trans_gtf, out_dir, out_file)
print (command_1)
os.system(command_1)
print ('Step 1 completed: input gtf has been converted to bed file.\n')

#### 2. translate transcripts into proteins ####
command_2 = "python %s/2_seq_translate_coordinate.py --mode %s --trans_inf %s/1_%s_gtf_processed.txt --abundance_inf %s --gtf_inf %s --outf_name %s/2_%s --genome %s --cds %s/Gencode_converted_CDS_annotation.txt" % (dir_path, mode, out_dir, out_file, abun_inf, input_trans_gtf, out_dir, out_file, genome, out_dir)
print (command_2)
os.system(command_2)
print ('Step 2 completed: input transcripts have been translated into proteins.\n')

#### 3. classify transcripts and rename ####
command_3 = "python %s/3_classify_transcript.py %s/2_%s_protein.txt %s %s %s" % (dir_path, out_dir, out_file, out_dir, ref_gtf, input_trans_gtf)
print (command_3)
os.system(command_3)
print ('Step 3 completed: translated proteins have been classified into different categories.\n')


