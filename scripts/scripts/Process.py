####################################################
#Goal:      Processing of ESPRESSO outputs
#Author:    Yang Xu
#E-mail:    yangax@pennmedicine.upenn.edu
####################################################

import os,sys,re
from collections import defaultdict
import configargparse

def parse_args():
	parser = configargparse.ArgParser(description='Processing of ESPRESSO outputs')
	parser.add_argument('-eg', '--espresso_gtf', dest='espresso_gtf', type=str, help='ESPRESSO gtf file', required=True)
	parser.add_argument('-ea', '--espresso_abundance', dest='espresso_abundance', type=str, help='ESPRESSO abundance file', required=True)
	parser.add_argument('-nm', '--normalized_mode', dest='normalized_mode', type=str, choices=['SAM','ESPRESSO'], help="Choose normalization mode from ['SAM','ESPRESSO']", required=True)
	parser.add_argument('-fs', '--folder_sam', dest='folder_sam', type=str, help='directory of corresponding sam files, sample names need to match abundance matrix', default = './')
	parser.add_argument('-od', '--outf_dir', dest='outf_dir', type=str, help='output directory', required=True)
	args = parser.parse_args()
	return parser, args

############ main script ##############
parser, args = parse_args()

dir_path = os.path.dirname(os.path.realpath(__file__))

espresso_gtf = args.espresso_gtf
espresso_abundance = args.espresso_abundance
folder_sam = args.folder_sam
normalized_mode = args.normalized_mode
outf_dir = args.outf_dir

# 1.1 convert gtf to bed file
cmd_process_1 = "python %s/1_1_gtf2bed_1_to_1_based.py %s %s %s" % (dir_path, espresso_gtf, espresso_abundance, outf_dir)
print (cmd_process_1)
os.system(cmd_process_1)
# 1.2 calculate CPM based on sam/ban files
cmd_process_2 = "python %s/1_2_calculate_CPM_based_on_bam.py %s %s %s" % (dir_path, espresso_abundance, folder_sam, normalized_mode)
print (cmd_process_2)
os.system(cmd_process_2)
# 1.3 merge to gene level, and generate files to filter out bad-quality genes
out_file_name = re.sub(".txt|.esp", "_CPM_%s.txt" % normalized_mode, espresso_abundance.split('/')[-1])
cmd_process_3 = "python %s/1_3_merge_to_gene_CPM.py %s/%s" % (dir_path, outf_dir, out_file_name)
print (cmd_process_3)
os.system(cmd_process_3)
# 1.4 calculate isoform proportion
cmd_process_4 = "python %s/1_4_calculate_isoform_proportion_vector.py %s/%s" % (dir_path, outf_dir, out_file_name)
print (cmd_process_4)
os.system(cmd_process_4)
