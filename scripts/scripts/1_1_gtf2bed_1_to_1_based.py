import os,re,sys
from collections import defaultdict

dir_path = os.path.dirname(os.path.realpath(__file__))
inf_name = sys.argv[1] #'CD44_sample_list_N2_R0_updated.gtf'
inf_abundance_name = sys.argv[2] # 'abundance.esp'
if re.findall("encode", inf_name.split('/')[-1]):
	key_word = 'gencode'
else:
	key_word = inf_name.split('/')[-1].split('_')[0]

if len(sys.argv)>2:
	dire_path = sys.argv[3]
else:
	dire_path = '/'.join(os.path.abspath(inf_name).split('/')[0:-1])

outf_name = dire_path+'/'+key_word+'_BedGraph.bed'
outf = open(outf_name, 'w')
outf.write('track name="%s" visibility=dense\n' % key_word)

def get_right_ID(raw_ID):
	if re.findall('_PAR_Y', raw_ID):
		return_ID = raw_ID.split('.')[0]+'-PAR-Y'
	else:
		return_ID = raw_ID.split('.')[0]
	return return_ID

def convert_list_to_str(trans_range,exon_list,trans_anno):
	chrom = trans_range.split('$')[0]
	strand = trans_range.split('$')[1]
	chromStart = int(trans_range.split('$')[2]) - 1 # 0-based
	chromEnd = int(trans_range.split('$')[3])
	thickStart = chromStart
	thickEnd = chromStart
	itemRgb = '0'
	name = exon_list[0].split('$')[0] #ENST_ID
	score = 0
	blockCount = len(exon_list)
	blockSizes = []
	blockStarts = []
	if len(exon_list) > 1:
		first_exon_start = int(exon_list[0].split('$')[3]) - 1
		second_exon_start = int(exon_list[1].split('$')[3]) - 1
		if first_exon_start > second_exon_start:
			exon_list = exon_list[::-1]
	for each_exon in exon_list:
		exon_start = int(each_exon.split('$')[3]) - 1
		exon_end = int(each_exon.split('$')[4])
		blockSizes.append(str(exon_end - exon_start))
		blockStarts.append(str(exon_start - chromStart))
	res = chrom+'\t'+str(chromStart)+'\t'+str(chromEnd)+'\t'+name+'\t'+str(score)+'\t'+strand+'\t'+str(thickStart)+'\t'+str(thickEnd)+'\t'+itemRgb+'\t'+str(blockCount)+'\t'+','.join(blockSizes)+'\t'+','.join(blockStarts)
	return res

#### load some annotation ####
trans2gene_dict = defaultdict()
'''
with open("%s/files/gencode.v39.annotation.gtf" % dir_path,"r") as gtf_hg38:
	for line in gtf_hg38:
		if line.startswith('#'): continue
		arr = line.strip().split('\t')
		if arr[2] != 'transcript': continue
		ENST_ID = get_right_ID(re.findall('transcript_id \"(.+?)\"',arr[8])[0])
		ENSG_ID = get_right_ID(re.findall('gene_id \"(.+?)\"',arr[8])[0])
		if ENST_ID not in trans2gene_dict:
			trans2gene_dict[ENST_ID] = ENSG_ID
'''
with open("%s/files/gencode.v34lift37.annotation.gtf" % dir_path,"r") as gtf_hg19:
	for line in gtf_hg19:
		if line.startswith('#'): continue
		arr = line.strip().split('\t')
		if arr[2] != 'transcript': continue
		ENST_ID = get_right_ID(re.findall('transcript_id \"(.+?)\"',arr[8])[0])
		ENSG_ID = get_right_ID(re.findall('gene_id \"(.+?)\"',arr[8])[0])
		if ENST_ID not in trans2gene_dict:
			trans2gene_dict[ENST_ID] = ENSG_ID

with open(inf_abundance_name, 'r') as abun_inf:
	for index, line in enumerate(abun_inf):
		if index == 0: continue
		arr = line.strip().split('\t')
		ENST_ID = get_right_ID(arr[0])
		ENSG_ID = get_right_ID(arr[2])
		if ENST_ID not in trans2gene_dict:
			trans2gene_dict[ENST_ID] = ENSG_ID

######################
key_list = []
prev_transcript_annotate = ''
prev_transcript_range = ''

first_exon_flag = 0
file_len_str = os.popen('wc -l %s' % inf_name).readlines()[0]
file_len = int(file_len_str.split(' ')[0])
inf = open(inf_name,'r')
for i,line in enumerate(inf):
	if line.startswith('#'): continue
	arr = line.strip().split('\t')
	if arr[2] == 'transcript':
		if len(key_list)>0:
			out_str = convert_list_to_str(prev_transcript_range,key_list,prev_transcript_annotate)
			outf.write(out_str+'\n')
		key_list = []
		prev_transcript_annotate = arr[1]
		prev_transcript_range = arr[0]+'$'+arr[6]+'$'+str(int(arr[3]))+'$'+str(arr[4])
		#first_exon_flag = 1
		continue
	elif arr[2] == 'exon':
		ENST_ID = get_right_ID(re.findall('transcript_id \"(.+?)\"',arr[8])[0])
		if re.findall("gene_id", arr[8]):
			ENSG_ID = get_right_ID(re.findall('gene_id \"(.+?)\"',arr[8])[0])
		elif ENST_ID in trans2gene_dict:
			ENSG_ID = trans2gene_dict[ENST_ID]
		else:
			continue
		key = ENST_ID+'|'+ENSG_ID+'$'+arr[0]+'$'+arr[6]+'$'+str(int(arr[3]))+'$'+str(arr[4])
		key_list.append(key)
		if i == file_len-1: # the end of entire file
			out_str = convert_list_to_str(prev_transcript_range,key_list,prev_transcript_annotate)
			outf.write(out_str+'\n')
inf.close()
outf.close()	
