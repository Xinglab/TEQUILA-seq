import os,re,sys
from collections import defaultdict
import numpy as np
import itertools

query_sample = sys.argv[1]
query_gene = sys.argv[2]
query_sample_bedgraph = sys.argv[6]
anno_gtf_inf_name = sys.argv[7]
query_transcript = []

if (re.findall('ENS',query_gene)) or (re.findall('SIRV',query_gene)):
	query_gene_name = 'NA'
	query_gene_ID = query_gene
else:
	query_gene_name = query_gene
	query_gene_ID = 'NA'

def get_right_ID(raw_ID):
	if re.findall('_PAR_Y', raw_ID):
		new_ID = raw_ID.split('.')[0] + '-PAR-Y'
	else:
		new_ID = raw_ID.split('.')[0]
	return new_ID

exp_dict = defaultdict()
ave_exp_dict = defaultdict()
sample_list = []
header_list = []
value_column = 3
abundance_inf_name = sys.argv[3] #CD44_sample_list_N2_R0_abundance_CPM.esp
key_word = abundance_inf_name.split('/')[-1].split('_')[0]
CPM_cut_off = float(sys.argv[4])  #default is 1.0
object_path = sys.argv[5]


exp_inf_0 = open(abundance_inf_name,'r')
for index,line in enumerate(exp_inf_0):
	arr = line.strip().split('\t')
	if index==0: continue
	else:
		ENST_ID = get_right_ID(arr[0])
		ENSG_ID = get_right_ID(arr[2])
		if ENSG_ID == 'NA':
			continue
		if (query_gene_name != 'NA') and (query_gene_ID == 'NA') :
			if re.findall(query_gene_name+'-\d+',arr[1]):
				query_gene_ID = ENSG_ID
				break
		elif (query_gene_ID != 'NA') and (query_gene_name == 'NA'):
			if ENSG_ID == query_gene_ID:
				if re.findall('SIRV',query_gene_ID):
					query_gene_name = query_gene_ID
				elif arr[1].startswith('ENS'):
					query_gene_name = query_gene_ID
				else:
					query_gene_name = arr[1].split('-')[0]
				break
exp_inf_0.close()
print (query_gene_name,query_gene_ID,'Cutoff:',CPM_cut_off)

exp_inf = open(abundance_inf_name,'r')
for index,line in enumerate(exp_inf):
	arr = line.strip().split('\t')
	if re.findall('transcript',arr[0]):
		header_list = arr
		for i in range(value_column,len(arr)):
			sample = arr[i]
			exp_dict[sample] = defaultdict()
			sample_list.append(sample)
	else:
		ENST_ID = get_right_ID(arr[0])
		ENSG_ID = get_right_ID(arr[2])
		if ENSG_ID == 'NA': continue
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
	ENST_ID = get_right_ID(arr[3].split('|')[0])
	if ENST_ID not in query_transcript:
		continue
	if ENST_ID in exp_dict[sample]:
		out_str = ''
		outf_separate_name = query_sample+'_'+query_gene_name+'_'+ENST_ID+'.bed'
		bed_separate_outf = open(object_path+'/target_genes/'+outf_separate_name,'w')

		if (re.findall('ENS',ENST_ID)) or (re.findall('SIRV',query_gene)):
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

#### get bed from annotated gtf: processing of annotated gtf #####
anno_gtf_outf_name = re.sub(".gtf", "_processed.gtf", anno_gtf_inf_name)
if os.path.exists(anno_gtf_outf_name):
	pass
else:
	anno_gtf_outf = open(anno_gtf_outf_name, 'w')
	coding_flag = 0
	CDS_list = []
	last_ENST_ID = ''
	gene_range_dict = defaultdict()
	trans2info_dict = defaultdict()
	with open(anno_gtf_inf_name, 'r') as gtf_inf:
		for line in gtf_inf:
			arr = line.strip().split('\t')
			if line.startswith('#'): continue
			if arr[2] == 'gene':
				ENSG_ID = get_right_ID(re.findall('gene_id \"(.+?)\";',arr[8])[0])
				gene_range = arr[0]+'#'+arr[6]+'#'+str(arr[3])+'#'+str(arr[4])
				gene_range_dict[ENSG_ID] = gene_range
				continue
			elif arr[2]=='transcript':
				if coding_flag == 1: #last transcript has not been done 
					info_arr = trans2info_dict[last_ENST_ID]
					if (info_arr[4] == '-') and (len(CDS_list)>=2):
						first_start = int(CDS_list[0].split('#')[0])
						second_start = int(CDS_list[1].split('#')[0])
						if first_start > second_start:
							CDS_list = CDS_list[::-1]
					CDS_region = ';'.join(CDS_list)
					CDS_list = []
					coding_flag = 0
					info_arr = trans2info_dict[last_ENST_ID]
					anno_gtf_outf.write(last_ENST_ID+'\t'+info_arr[2]+'\t'+info_arr[3]+'\t'+info_arr[0]+'\t'+info_arr[1]+'\t'+gene_range_dict[info_arr[0]]+'\t'+CDS_region+'\n')
				# new transcript
				ENST_ID = get_right_ID(re.findall('transcript_id \"(.+?)\";',arr[8])[0])
				ENSG_ID = get_right_ID(re.findall('gene_id \"(.+?)\";',arr[8])[0])
				gene_name = '-'
				if re.findall('gene_name \"(.+?)\";',arr[8]):
					gene_name = re.findall('gene_name \"(.+?)\";',arr[8])[0]
				trans_type = '-'
				if re.findall('transcript_type \"(.+?)\";',arr[8]):
					trans_type = re.findall('transcript_type \"(.+?)\";',arr[8])[0]
				trans_tag = ';'.join(re.findall('tag \"(.+?)\";',arr[8]))
				strand = arr[6]
				coding_flag = 1
				trans2info_dict[ENST_ID] = [ENSG_ID, gene_name, trans_type, trans_tag, strand]
				last_ENST_ID = ENST_ID
				continue
			#elif (coding_flag==1) and (arr[2]=='CDS'):  #only show CDS resgion	
			elif (coding_flag==1) and (arr[2]=='exon'): #show whole transcript
				CDS_exon = str(arr[3])+'#'+str(arr[4])
				CDS_list.append(CDS_exon)
	anno_gtf_outf.close()


##### get bed from annotated gtf #####
with open(anno_gtf_outf_name, 'r') as anno_gtf:
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

