import os,re,sys
from collections import defaultdict
import numpy as np

Total_e_reads_list = []
value_column = 3
sample_list = []
inf_name = sys.argv[1]
sam_folder = sys.argv[2]

mode = "ESPRESSO"
if len(sys.argv) > 3: 
	mode = sys.argv[3]

with open(inf_name) as inf_1:
	for index, line in enumerate(inf_1):
		arr = line.strip().split('\t')
		if index == 0:
			sample_list = arr[value_column:len(arr)]
			Total_e_reads_list = np.array([0.0]*(len(arr)-value_column))
		else:
			Total_e_reads_list += np.array(list(map(float,arr[value_column:len(arr)])))

if mode == "ESPRESSO":
	#print ('Total_e_reads_list (ESPRESSO):', Total_e_reads_list)
	Total_reads_list = np.array(list(map(float, Total_e_reads_list)))
elif mode in ['Sam','sam','SAM']:
	mode = 'sam'
	Total_s_reads_list = []
	for each_sample in sample_list:
		sample_sam_file = sam_folder+'/'+each_sample+'/0/aligned.sorted.sam'
		command_2 = "awk -F'\t' '{if($5>=1) print $1}' %s |sort -u| wc -l" % sample_sam_file
		reads_genome = os.popen(command_2).readlines()[0].strip().split(' ')[0]
		Total_s_reads_list.append(reads_genome)
	Total_reads_list = np.array(list(map(float, Total_s_reads_list)))
	#print ('Total_s_reads_list (sam):', Total_s_reads_list)



dire_path = '/'.join(os.path.abspath(inf_name).split('/')[0:-1])
outf_name = dire_path+'/'+re.sub(".esp$|.txt$", '_CPM_%s.txt'%mode, inf_name.split('/')[-1])
outf_3 = open(outf_name,'w')
inf_3 = open(inf_name,'r')
with open(inf_name) as inf_3:
	for index,line in enumerate(inf_3):
		arr = line.strip().split('\t')
		if index == 0:
			outf_3.write(line)
			continue
		else:
			value_list = np.array(list(map(float,arr[value_column:len(arr)])))*1000000/Total_reads_list
			outf_3.write('\t'.join(arr[0:value_column])+'\t'+'\t'.join(list(map(str,value_list)))+'\n')	
outf_3.close()

