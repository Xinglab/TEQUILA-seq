####################################################
#Goal:      Generation of figures
#Author:    Yang Xu
#E-mail:    yangax@pennmedicine.upenn.edu
####################################################

import os,sys,re
from collections import defaultdict
import configargparse

def parse_args():
	parser = configargparse.ArgParser(description='Generation of figures')
	parser.add_argument('-ip', '--isoform_proportion_inf', dest='isoform_proportion_inf', type=str, help='Isoform proportion infile', required=True)
	parser.add_argument('-ic', '--isoform_cpm_inf', dest='isoform_cpm_inf', type=str, help='Isoform CPM infile', required=True)
	parser.add_argument('-gi', '--group_info_inf', dest='group_info_inf', type=str, help='Group information file', required=True)
	parser.add_argument('-rt', '--required_trans_inf', dest='required_trans_inf', type=str, help='Required transcripts information file', required=True)
	parser.add_argument('-be', '--bedgraph', dest='bedgraph', type=str, help='Generated bedgraph file for sample', required=True)
	parser.add_argument('-gv', '--genome_version', dest='genome_version', type=str, choices=['GRCH38','GRCH37','hg38','hg19'], help="choose from ['GRCH38','GRCH37','hg38','hg19']", default='hg38')
	parser.add_argument('-fi', '--figures', dest='figures', nargs='+', type=str, choices=['Isoform','Single_isoform','Structure'], default=['Isoform','Single_isoform','Structure'],  help='Choose from [Isoform, Single_isoform, Structure]', required=True)
	parser.add_argument('-of', '--order', dest='order', type=str, help='order samples by isoform proportion', default='no')
	parser.add_argument('-od', '--outf_dir', dest='outf_dir', type=str, help='output directory', required=True)

	args = parser.parse_args()
	return parser, args

############ main script ##############
parser, args = parse_args()

dir_path = os.path.dirname(os.path.realpath(__file__))

isoform_proportion_inf = args.isoform_proportion_inf
isoform_cpm_inf = args.isoform_cpm_inf
group_info_inf = args.group_info_inf
required_trans_inf = args.required_trans_inf
bedgraph = args.bedgraph
outf_dir = args.outf_dir
figures = args.figures
order = args.order
genome_version = args.genome_version

# 2.1 only show top 5 isoform's proportion and reshape the file format
cmd_figure_1 = "python %s/2_1_merge_isoforms_and_reshape_format.py %s %s %s %s %s" % (dir_path, isoform_proportion_inf, isoform_cpm_inf, group_info_inf, required_trans_inf, outf_dir)
print(cmd_figure_1)
os.system(cmd_figure_1)

prop_reshaped_inf_name = outf_dir+'/'+re.sub(".txt", "_reshaped_merge_others.txt", isoform_proportion_inf.split('/')[-1])
exp_reshaped_inf_name = outf_dir+'/'+re.sub(".esp|.txt", "_reshaped_merge_others.txt", isoform_cpm_inf.split('/')[-1])
sorted_group_list = []
with open(group_info_inf,'r') as inf:
	for index, line in enumerate(inf):
		if index == 0: continue
		arr = line.strip().split("\t")
		sorted_group_list.append(arr[0])
sorted_group_info = ','.join(sorted_group_list)
# 2.2 Generate bargraph and isoform structure figures [one example]
if os.path.exists("%s/Example_res" % outf_dir):
	pass
else:
	os.system("mkdir %s/Example_res" % outf_dir)
if genome_version in ['hg19','GRCH37']:
	cmd_figure_2 = "python %s/2_2_Generate_bar_structure_figure_example.py --gene [Ensembl_Gene_ID] --gene_name [Gene_Symbol] --transcript [Interested_Transcript_ID] --abundance_CPM_original %s --abundance_proportion %s --abundance_CPM %s --bedgraph %s --sorted_group %s --out_dir %s/Example_res --figures %s --canonical_transcript %s/files/Gencode_v39_canonical_isoform.txt --basic_transcript %s/files/gencode.v34lift37.annotation_basic_trans.txt --anno_gtf %s/files/gencode.v34lift37.annotation.gtf --genome_version %s --order %s" % (dir_path, isoform_cpm_inf, prop_reshaped_inf_name, exp_reshaped_inf_name, bedgraph, sorted_group_info, outf_dir, ' '.join(figures), dir_path, dir_path, dir_path, genome_version, order)
elif genome_version in ['hg38', 'GRCH38']:
	cmd_figure_2 = "python %s/2_2_Generate_bar_structure_figure_example.py --gene [Ensembl_Gene_ID] --gene_name [Gene_Symbol] --transcript [Interested_Transcript_ID] --abundance_CPM_original %s --abundance_proportion %s --abundance_CPM %s --bedgraph %s --sorted_group %s --out_dir %s/Example_res --figures %s --canonical_transcript %s/files/Gencode_v39_canonical_isoform.txt --basic_transcript %s/files/gencode.v39.annotation_basic_trans.txt --anno_gtf %s/files/gencode.v39.annotation.gtf --genome_version %s --order %s" % (dir_path, isoform_cpm_inf, prop_reshaped_inf_name, exp_reshaped_inf_name, bedgraph, sorted_group_info, outf_dir, ' '.join(figures), dir_path, dir_path, dir_path, genome_version, order)
outf_figure_2 = open("%s/Template_to_generate_figures.sh" % outf_dir,"w")
outf_figure_2.write(cmd_figure_2+'\n')
outf_figure_2.close()
