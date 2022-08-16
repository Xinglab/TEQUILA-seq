import os,re,sys
from collections import defaultdict
import numpy as np

bed_dict = defaultdict()
trans_list = []
bed_inf_name = sys.argv[1]
gene_name = sys.argv[2]
name_col = sys.argv[3] # 2
with open(bed_inf_name, 'r') as list_inf:
	for list_line in list_inf:
		file_name = list_line.strip().split('/')[-1].split('.')[0]
		trans_ID = '_'.join(file_name.split('_')[int(name_col):len(file_name.split('_'))])
		bed_dict[trans_ID] = list_line.strip()
		if trans_ID not in trans_list:
			trans_list.append(trans_ID)

header_list = []
sample_list = []
CPM_group_dict = defaultdict()
CPM_total_dict = defaultdict()
abundance_inf_name = sys.argv[4]
gtf_inf_name = '/'.join(abundance_inf_name.split('/')[0:-1])+'/samples_N2_R0_updated_combined_with_tissue.gtf'
if os.path.exists(gtf_inf_name):
	pass
else:
	gtf_inf_name = '/'.join(abundance_inf_name.split('/')[0:-1])+'/samples_N2_R0_updated.gtf'

with open(abundance_inf_name, 'r') as inf:
	for index,line in enumerate(inf):
		arr = line.strip().split('\t')
		if index == 0: 
			header_list = arr
			sample_list = arr[3:len(arr)]
			continue
		trans_ID = arr[0].split('.')[0]
		if trans_ID not in trans_list: continue
		for i in range(3,len(arr)):
			sample = header_list[i]
			if trans_ID not in CPM_group_dict:
				CPM_group_dict[trans_ID] = defaultdict()
			CPM_group_dict[trans_ID][sample] = float(arr[i])
		CPM_total_dict[trans_ID] = np.nanmean(list(map(float,arr[3:len(arr)])))


#### detect isoforms (anno+novel) ####
common_strand = ''
strand_dict = defaultdict()
with open(gtf_inf_name, 'r') as gtf_inf:
	for line in gtf_inf:
		if line.startswith('#'): continue
		arr = line.strip().split('\t')
		if arr[2] == 'transcript':
			trans_ID = re.findall('transcript_id \"(.+?)\";',arr[8])[0].split('.')[0]
			if trans_ID in trans_list:
				strand_dict[trans_ID] = arr[6]
				common_strand = arr[6]

for each_trans in trans_list:
	if each_trans not in CPM_total_dict:
		CPM_total_dict[each_trans] = 0
		CPM_group_dict[each_trans] = defaultdict()
		for each_sample in sample_list:
			CPM_group_dict[each_trans][each_sample] = 0
#	if each_trans not in strand_dict:
#		strand_dict[each_trans] = common_strand

##### undetected anno isoforms #####
tag_dict = defaultdict()
#with open('/mnt/isilon/xing_lab/aspera/xuy/snakemake_ESPRESSO_reference/Translation_scripts/files/gencode.v39.annotation.gtf','r') as gencode_gtf:
with open('/mnt/isilon/xing_lab/aspera/xuy/snakemake_ESPRESSO_reference/references/gencode.v34lift37.annotation.gtf','r') as gencode_gtf:
	for line in gencode_gtf:
		if line.startswith('#'): continue
		arr = line.strip().split('\t')
		if arr[2] == 'transcript':
			trans_ID = re.findall('transcript_id \"(.+?)\";',arr[8])[0].split('.')[0]
			if trans_ID in trans_list:
				if trans_ID not in strand_dict:
					strand_dict[trans_ID] = arr[6]
				trans_type = re.findall('transcript_type \"(.+?)\";',arr[8])[0]
				trans_tag = ';'.join(re.findall('tag \"(.+?)\";',arr[8]))
				if trans_type == 'nonsense_mediated_decay':
					trans_type_final = 'NMD'
				else:
					trans_type_final = trans_type
				if re.findall('basic',trans_tag):
					trans_tag_final = 'basic'
				elif re.findall('mRNA_start_NF', trans_tag):
					trans_tag_final = 'mRNA_start_NF'
				elif re.findall('mRNA_end_NF', trans_tag):
					trans_tag_final = 'mRNA_end_NF'
				else:
					trans_tag_final = trans_tag
				if trans_tag_final != '':
					tag_dict[trans_ID] = trans_type_final+'; '+trans_tag_final
				else:
					tag_dict[trans_ID] = trans_type_final



sorted_list = sorted(CPM_total_dict.items(), key=lambda x:x[1], reverse=True)
sorted_list_ascend = sorted(CPM_total_dict.items(), key=lambda x:x[1], reverse=False)

count_novel = 0
selected_novel_list = []
for each_item in sorted_list:
	each_trans_ID = each_item[0]
	if each_trans_ID.startswith('ESPRESSO'):
		flag = 'Novel'
		count_novel += 1
		#if count_novel > 5:
		#	continue
		selected_novel_list.append(each_trans_ID)


exp_outf_name = re.sub('bed_list.txt','exp_data.txt',bed_inf_name)
outf = open(exp_outf_name, 'w')
outf.write('Transcript\tGroup\tCPM\tType\n')
sorted_bed_name = re.sub('.txt','_sorted.txt', bed_inf_name)
outf_bed = open(sorted_bed_name, 'w')

for each_item in sorted_list_ascend:
	each_trans_ID = each_item[0]
	if each_trans_ID.startswith('ESPRESSO'):
		flag = 'Novel'
		if each_trans_ID not in selected_novel_list:
			continue
		else:
			tag_dict[each_trans_ID] = ''
	else:
		flag = 'Annotated'
	for each_sample in sample_list:
		outf.write(each_trans_ID.split('#')[0]+'\t'+each_sample+'\t'+str(CPM_group_dict[each_trans_ID][each_sample])+'\t'+flag+'\n')
	outf_bed.write(bed_dict[each_trans_ID]+'\t'+strand_dict[each_trans_ID]+'\t'+tag_dict[each_trans_ID]+'\n')
outf.close()
outf_bed.close()


'''
##### isoform structure plot #####
group_col = int(name_col) - 1
command = 'python /mnt/isilon/xing_lab/aspera/xuy/snakemake_ESPRESSO_reference/pipeline_test/scripts_figure/Isoform_structure_plot.py %s %s %s ESPRESSO' % (sorted_bed_name, gene_name, group_col)
#print (command)
#os.system(command)

##### exp bar plot ######
exp_plot_inf_name = exp_outf_name
exp_plot_outf_name = re.sub('.txt','.png',exp_plot_inf_name)
command_2 = 'Rscript --vanilla /mnt/isilon/xing_lab/aspera/xuy/snakemake_ESPRESSO_reference/pipeline_test/scripts_figure/Fig1E_exp_data.R %s %s' % (exp_plot_inf_name, exp_plot_outf_name)
#print (command_2)
#os.system(command_2)
'''
