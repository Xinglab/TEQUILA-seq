####################################################
#Goal:      Translation
#Author:    Yang Xu
#E-mail:    yangax@pennmedicine.upenn.edu
####################################################

import os,sys,re
from collections import defaultdict
import configargparse

def parse_args():
	parser = configargparse.ArgParser(description='Translation')
	parser.add_argument('-mo', '--mode', dest='mode', type=str, default = 'long-read', help='Long-read RNA-seq data mode', required=True)
	parser.add_argument('-tg', '--trans_gtf', dest='trans_gtf', type=str, help='generated gtf file', required=True)
	parser.add_argument('-ic', '--isoform_cpm_inf', dest='isoform_cpm_inf', type=str, help='Isoform CPM infile', required=True)
	parser.add_argument('-gv', '--genome_version', dest='genome_version', type=str, choices=['GRCH38','GRCH37','hg38','hg19'], help="choose from ['GRCH38','GRCH37','hg38','hg19']", required=True)
	parser.add_argument('-rg', '--ref_gtf', dest='ref_gtf', type=str, help='reference gencode annotation', required=True)
	parser.add_argument('-of', '--out_file', dest='out_file', type=str, help='prefix of the name of output file', required=True)
	parser.add_argument('-od', '--outf_dir', dest='outf_dir', type=str, help='output directory', required=True)
	args = parser.parse_args()
	return parser, args

############ main script ##############
parser, args = parse_args()

dir_path = os.path.dirname(os.path.realpath(__file__))

#genome = args.genome
if args.genome_version in ['hg19','GRCH37','grch37']:
	genome = '%s/files/GRCh37.primary_assembly.genome.fa' % dir_path
elif args.genome_version in ['hg38','GRCH38','grch38']:
	genome = '%s/files/GRCh38.primary_assembly.genome.fa' % dir_path
input_trans_gtf = args.trans_gtf
ref_gtf = args.ref_gtf
outf_dir = args.outf_dir
out_file = args.out_file
isoform_cpm_inf = args.isoform_cpm_inf
mode = args.mode

#### 0. process reference gtf file to generate CDS annotation ####
command_0 = "python %s/4_1_process_gtf.py %s %s" % (dir_path, ref_gtf, outf_dir)
print (command_0)
os.system(command_0)
print ('Step 0 completed: reference gtf has been converted to CDS annotation file.\n')

#### 1. convert input gtf file into bed file ####
command_1 = "python %s/4_2_convert_gtf2transcript_1_to_1_based.py %s %s/4_2_%s_gtf_processed.txt" % (dir_path, input_trans_gtf, outf_dir, out_file)
print (command_1)
os.system(command_1)
print ('Step 1 completed: input gtf has been converted to bed file.\n')

#### 2. translate transcripts into proteins ####
command_2 = "python %s/4_3_seq_translate_coordinate.py --mode %s --trans_inf %s/4_2_%s_gtf_processed.txt --abundance_inf %s --gtf_inf %s --outf_name %s/4_3_%s --genome %s --cds %s/4_1_Gencode_converted_CDS_annotation.txt" % (dir_path, mode, outf_dir, out_file, isoform_cpm_inf, input_trans_gtf, outf_dir, out_file, genome, outf_dir)
print (command_2)
os.system(command_2)
print ('Step 2 completed: input transcripts have been translated into proteins.\n')

#### 3. classify transcripts and rename ####
command_3 = "python %s/4_4_classify_transcript.py %s/4_3_%s_protein.txt %s %s %s %s" % (dir_path, outf_dir, out_file, outf_dir, ref_gtf, input_trans_gtf, out_file)
print (command_3)
os.system(command_3)
print ('Step 3 completed: translated proteins have been classified into different categories.\n')

