import os,re,sys
from collections import defaultdict
import numpy as np
import itertools

query_sample = sys.argv[1]
query_gene = sys.argv[2]
query_sample_bedgraph = sys.argv[6]
query_transcript = []

if (re.findall('ENSG',query_gene)) or (re.findall('SIRV',query_gene)):
	query_gene_name = 'NA'
	query_gene_ID = query_gene
else:
	query_gene_name = query_gene
	query_gene_ID = 'NA'

exp_dict = defaultdict()
ave_exp_dict = defaultdict()
sample_list = []
header_list = []
value_column = 3
abundance_inf_name = sys.argv[3] #CD44_sample_list_N2_R0_abundance_CPM.esp
key_word = abundance_inf_name.split('/')[-1].split('_')[0]
object_path = '/'.join(os.path.abspath(abundance_inf_name).split('/')[0:-1])+'/statistics'
CPM_cut_off = float(sys.argv[4])  #default is 1.0

if len(sys.argv) > 5:
	object_path = sys.argv[5]


exp_inf_0 = open(abundance_inf_name,'r')
for line in exp_inf_0:
	arr = line.strip().split('\t')
	if re.findall('transcript',arr[0]): continue
	else:
		ENST_ID = arr[0].split('.')[0]
		ENSG_ID = arr[2].split('.')[0]
		if ENSG_ID == 'NA':
			continue
		if (query_gene_name != 'NA') and (query_gene_ID == 'NA') :
			if re.findall(query_gene_name+'-\d+',arr[1]):
				query_gene_ID = arr[2].split('.')[0]
				break
		elif (query_gene_ID != 'NA') and (query_gene_name == 'NA'):
			if ENSG_ID == query_gene_ID:
				query_gene_name = arr[1].split('-')[0]
				if re.findall('SIRV',query_gene_ID):
					query_gene_name = query_gene_ID
				break
exp_inf_0.close()
print (query_gene_name,query_gene_ID,'Cutoff:',CPM_cut_off)

exp_inf = open(abundance_inf_name,'r')
for line in exp_inf:
	arr = line.strip().split('\t')
	if re.findall('transcript',arr[0]):
		header_list = arr
		for i in range(value_column,len(arr)):
			sample = arr[i]
			exp_dict[sample] = defaultdict()
			sample_list.append(sample)
	else:
		ENST_ID = arr[0].split('.')[0]
		ENSG_ID = arr[2].split('.')[0]
		if ENSG_ID == 'NA':
			continue
		trans_name = arr[1]
		if ENSG_ID == query_gene_ID:
			query_transcript.append(ENST_ID)
			if (np.max(list(map(float,arr[value_column:len(arr)]))) >= CPM_cut_off) or (re.findall('ENST', ENST_ID)):
				for i in range(value_column,len(arr)):
					sample = header_list[i]
					exp_dict[sample][ENST_ID] = float(arr[i])
				ave_exp_dict[ENST_ID] = np.mean(list(map(float,arr[value_column:len(arr)])))
exp_inf.close()

#print exp_dict[sample_list[0]]['ENST00000655278']

def get_region(chromStart,block_size,block_start):
	block_start_list = block_start.split(',')
	block_size_list = block_size.split(',')
	res_list = []
	for i in range(len(block_start_list)):
		block = str(int(chromStart)+int(block_start_list[i]))+'-'+str(int(chromStart)+int(block_start_list[i])+int(block_size_list[i]))
		res_list.append(block)
	return res_list	

#bed_inf = open(object_path+'/'+key_word+'_BedGraph.bed','r')
bed_inf = open(query_sample_bedgraph, 'r')
outf_name = query_sample+'_'+query_gene_name+'.bed'

if os.path.exists(object_path+'/target_genes'):
	pass
else:
	os.system('mkdir '+object_path+'/target_genes')

#bed_graph_outf = open(object_path+'/target_genes/'+outf_name,'w')
out_str_dict = defaultdict()

for line in bed_inf:
	arr = line.strip().split('\t')
	if re.findall('track name',line):
		continue
	if re.findall('chrUn',arr[0]):
		continue
	block_region = get_region(arr[1], arr[-2], arr[-1])	
	sample = query_sample
	ENST_ID = arr[3].split('.')[0]
	if ENST_ID not in query_transcript:
		continue
	if ENST_ID in exp_dict[sample]:
		out_str = ''
		outf_separate_name = query_sample+'_'+query_gene_name+'_'+ENST_ID+'.bed'
		bed_separate_outf = open(object_path+'/target_genes/'+outf_separate_name,'w')

		if (re.findall('ENST',ENST_ID)) or (re.findall('SIRV',query_gene)):
			out_str = 'track type=bedGraph name="%s" visibility=dense color=26,37,187\n' % (sample+'_'+ENST_ID)
		else:
			out_str = 'track type=bedGraph name="%s" visibility=dense color=211,37,17\n' % (sample+'_'+ENST_ID)
		for each_block in block_region:
			b_start = int(each_block.split('-')[0])
			b_end = int(each_block.split('-')[1])
			b_exp = round(float(exp_dict[sample][ENST_ID]),2)
			out_str = out_str + arr[0]+'\t'+str(b_start)+'\t'+str(b_end)+'\t'+str(b_exp)+'\n'
		#out_str_dict[sample+'_'+ENST_ID] = out_str
		bed_separate_outf.write(out_str)
		bed_separate_outf.close()
	
bed_inf.close()

##### get bed from annotated gtf #####
anno_gtf_inf_name = '../../Translation_scripts/files/Gencode_converted_all_exon_annotation_v39.txt'
with open(anno_gtf_inf_name, 'r') as anno_gtf:
	for line in anno_gtf:
		arr = line.strip().split('\t')
		if arr[3].split('.')[0] == query_gene_ID:
			ENST_ID = arr[0].split('.')[0]
			if ENST_ID not in exp_dict[query_sample]:
				out_str = ''
				outf_separate_name = query_sample+'_'+query_gene_name+'_'+ENST_ID+'.bed'
				bed_separate_outf_2 = open(object_path+'/target_genes/'+outf_separate_name,'w')

				if (re.findall('ENST',ENST_ID)) or (re.findall('SIRV',query_gene)):
					out_str = 'track type=bedGraph name="%s" visibility=dense color=26,37,187\n' % (sample+'_'+ENST_ID)
				else:
					out_str = 'track type=bedGraph name="%s" visibility=dense color=211,37,17\n' % (sample+'_'+ENST_ID)
				for each_block in arr[6].split(';'):
					b_start = int(each_block.split('#')[0])-1
					b_end = int(each_block.split('#')[1])
					b_exp = round(0)
					out_str = out_str + arr[0]+'\t'+str(b_start)+'\t'+str(b_end)+'\t'+str(b_exp)+'\n'
				#out_str_dict[sample+'_'+ENST_ID] = out_str
				bed_separate_outf_2.write(out_str)
				bed_separate_outf_2.close()

