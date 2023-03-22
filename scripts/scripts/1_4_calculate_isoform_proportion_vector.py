import re,sys,os
from collections import defaultdict
import numpy as np

def get_right_ID(raw_ID):
	if re.findall('_PAR_Y', raw_ID):
		return raw_ID.split('.')[0]+'-PAR-Y'
	else:
		return raw_ID.split('.')[0]

gene_CPM_dict = defaultdict(lambda: defaultdict(lambda: 0))
sample_list = []
inf_name = sys.argv[1]
with open(inf_name,'r') as inf:
	for index,line in enumerate(inf):
		arr = line.strip().split('\t')
		if index == 0:
			sample_list = arr[3:len(arr)] 
			gene_CPM_dict = defaultdict(lambda: np.array([0.00001]*(len(arr)-3)))
			continue
		if arr[2] == 'NA': continue
		if re.findall(',', arr[2]): continue
		trans_ID = get_right_ID(arr[0])
		gene_ID = get_right_ID(arr[2])
		max_CPM = np.max(list(map(float,arr[3:len(arr)])))
		if max_CPM <= 0: continue
		CPM_list = np.array(list(map(float, arr[3:len(arr)])))
		gene_CPM_dict[gene_ID] += CPM_list

outf_name = re.sub(".esp$|.txt$",'_proportion.txt',inf_name)
outf = open(outf_name, 'w')
with open(inf_name,'r') as inf2:
	for index,line in enumerate(inf2):
		arr = line.strip().split('\t')
		if index == 0:
			outf.write('Transcript_ID\tGene_ID\t'+'\t'.join(arr[3:len(arr)])+'\n')
			continue
		if arr[2] == 'NA': continue
		if re.findall(',', arr[2]): continue
		trans_ID = get_right_ID(arr[0])
		gene_ID = get_right_ID(arr[2])
		CPM_list = np.array(list(map(float, arr[3:len(arr)])))
		prop_list = CPM_list*100 / gene_CPM_dict[gene_ID]
		outf.write(trans_ID+'\t'+gene_ID)
		for i in range(len(sample_list)):
			outf.write('\t'+str(round(prop_list[i],2)))
		outf.write('\n')
outf.close()


