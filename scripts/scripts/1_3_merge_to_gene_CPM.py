import os,re,sys
from collections import defaultdict
import numpy as np

gene_CPM = defaultdict(lambda: defaultdict(lambda: 0))
anno_CPM = defaultdict(lambda: defaultdict(lambda: 0))

value_column = 3
inf_name = sys.argv[1]
sample_list = []
inf_2 = open(inf_name)
for index,line in enumerate(inf_2):
	arr = line.strip().split('\t')
	if index == 0: 
		sample_list = arr[value_column:len(arr)]
		continue
	else:
		if arr[2] == 'NA': continue
		if re.findall('_PAR_Y', arr[2]):
			gene_ID = arr[2].split('.')[0]+'-PAR-Y'
		else:
			gene_ID = arr[2].split('.')[0]
		for i in range(value_column, len(arr)):
			each_sample = sample_list[i-value_column]
			gene_CPM[gene_ID][each_sample] += float(arr[i])
			if arr[0].startswith("ENS"):
				anno_CPM[gene_ID][each_sample] += float(arr[i])
inf_2.close()

dire_path = '/'.join(os.path.abspath(inf_name).split('/')[0:-1])
outf_name = dire_path+'/'+re.sub(".esp$|.txt$","_gene.txt",inf_name.split('/')[-1])
outf_name_anno_percentage = dire_path+'/'+re.sub(".esp$|.txt$","_gene_annotated_isoform_contribution.txt",inf_name.split('/')[-1])
outf_3 = open(outf_name,'w')
outf_3.write('Gene_ID\t'+'\t'.join(sample_list)+'\n')
outf_4 = open(outf_name_anno_percentage,'w')
outf_4.write('Gene_ID\t'+'\t'.join(sample_list)+'\n')
for each_gene in gene_CPM:
	outf_3.write(each_gene)
	outf_4.write(each_gene)
	for each_sample in sample_list:
		outf_3.write('\t'+str(gene_CPM[each_gene][each_sample]))
		if float(gene_CPM[each_gene][each_sample]) == 0:
			anno_percent = 0
		else:
			anno_percent = round(float(anno_CPM[each_gene][each_sample])*100/float(gene_CPM[each_gene][each_sample]), 2)
		outf_4.write('\t'+str(anno_percent))
	outf_3.write('\n')
	outf_4.write('\n')
outf_3.close()
outf_4.close()

