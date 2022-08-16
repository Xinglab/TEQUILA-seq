import os,re,sys
from collections import defaultdict

abundance_index_dict = defaultdict()
bar_plot_isoform_list = []
isoform2index = defaultdict()
top_isoform = sys.argv[2].rstrip('.txt').split('/')[-1].split('_')[4]
with open(sys.argv[2],'r') as bar_plot_inf:
	for index,line in enumerate(bar_plot_inf):
		arr = line.strip().split('\t')
		if arr[0] != 'Others':
			bar_plot_isoform_list.append(arr[0])
			isoform2index[arr[0]] = index
		abundance_index_dict[arr[0]] = int(arr[1])


rank_isoform_list = ['']*len(bar_plot_isoform_list)
with open(sys.argv[1],'r') as inf:
	for line in inf:
		arr = line.strip().split('\t')
		trans_ID = arr[0].split('.bed')[0].split('_')[-1]
		if trans_ID in isoform2index:
			line_str = '\t'.join(arr)
			rank_isoform_list[isoform2index[trans_ID]] = str(abundance_index_dict[trans_ID]) +'\t'+ line_str

outf_name = re.sub('.txt','_'+top_isoform+'_rank.txt',sys.argv[1])
outf = open(outf_name,'w')
for each_line_str in rank_isoform_list[::-1]:
	outf.write(each_line_str+'\n')
outf.close()
