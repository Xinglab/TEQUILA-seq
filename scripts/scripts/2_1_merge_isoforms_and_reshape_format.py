import os,re,sys
import numpy as np
from collections import defaultdict

proportion_inf_name = sys.argv[1]
CPM_inf_name = sys.argv[2]
sample_type_inf_name = sys.argv[3]
required_trans_inf_name = sys.argv[4]
outf_dir = sys.argv[5]

file_dir = os.path.dirname(os.path.realpath(__file__))
if len(sys.argv) > 6:
	canonical_trans_inf_name = sys.argv[6]
else:
	canonical_trans_inf_name = file_dir+"/files/Gencode_v39_canonical_isoform.txt"

sample_type_list = []
subtype_list = []
########## load number of samples for each subtype #############
with open(sample_type_inf_name,'r') as sample_type_inf:
	for index,line in enumerate(sample_type_inf):
		arr = line.strip().split('\t')
		if index == 0: continue
		sample_type_list += [arr[0]]*int(arr[1])
		subtype_list.append(arr[0])

## add some required lowly expressed isform ###
#required_dict["ENSG00000127990"]=["ESPRESSO:chr7:10157:1061@Melanoma","ENST00000648936"] #SGCE
required_dict = defaultdict(lambda: [])
with open(required_trans_inf_name, 'r') as required_inf:
	for index,line in enumerate(required_inf):
		arr = line.strip().split('\t')
		if index == 0: continue
		for each_trans in arr[1].split(';'):
			required_dict[arr[0]].append(each_trans)

with open(canonical_trans_inf_name, "r") as canonical_trans_inf:
	for index,line in enumerate(canonical_trans_inf):
		arr = line.strip().split('\t')
		if index == 0: continue
		required_dict[arr[0]].append(arr[2])

#### change format ########
column = 2
cell2subtype_dict = defaultdict()
sample_list_2 = []
ave_pro_dict = defaultdict()
outf_2_name = re.sub(".txt", "_reshaped.txt", proportion_inf_name.split('/')[-1])
outf_2 = open(outf_dir+'/'+outf_2_name, "w")
outf_2.write("Gene_ID\tTranscript_ID\tSample\tProportion\tSubtype\n")
with open(proportion_inf_name, "r") as inf_2:
	for index,line in enumerate(inf_2):
		arr = line.strip().split("\t")
		if index == 0:
			sample_list_2 = arr[column:len(arr)]
			for i in range(len(sample_list_2)):
				cell2subtype_dict[sample_list_2[i]] = sample_type_list[i]
			continue
		if arr[1] == "NA": continue
		if re.findall("_PAR_Y", arr[2]):
			gene_ID = arr[1].split(".")[0]+"-PAR-Y"
		else:
			gene_ID = arr[1].split(".")[0]
		if re.findall("_PAR_Y", arr[0]):
			trans_ID = arr[0].split(".")[0]+"-PAR-Y"
		else:
			trans_ID = arr[0].split(".")[0]
		if gene_ID not in ave_pro_dict:
			ave_pro_dict[gene_ID] = defaultdict()
		ave_pro_dict[gene_ID][trans_ID] = np.nanmean(np.array(list(map(float,arr[column:len(arr)]))))
		for i in range(column,len(arr)):
			this_sample = sample_list_2[i-column]
			this_subtype = cell2subtype_dict[this_sample]
			outf_2.write(gene_ID+"\t"+trans_ID+"\t"+this_sample+"\t"+str(arr[i])+"\t"+this_subtype+"\n")
outf_2.close()

sample_list_defined = sample_list_2 # could reorder samples

top_exp_trans_dict = defaultdict(lambda: [])
for each_gene in ave_pro_dict:
	if each_gene in required_dict:
		for each_value in required_dict[each_gene]:
			if each_value not in top_exp_trans_dict[each_gene]:
				top_exp_trans_dict[each_gene].append(each_value)
	sorted_key = sorted(ave_pro_dict[each_gene].items(), key=lambda x:float(x[1]), reverse = True)
	for index,each_item in enumerate(sorted_key):
		if len(top_exp_trans_dict[each_gene]) >= 5:
			break
		else:
			if each_item[0] not in top_exp_trans_dict[each_gene]:
				top_exp_trans_dict[each_gene].append(each_item[0])

############### merged low expressed isoforms ###########
subtype2color_dict = defaultdict()
color_list = ["\"#FF0018\"", "\"#e68e19\"", "\"#02a620\"", "\"#0000F9\"","\"#86007D\""] * 10
#subtype2color_dict['Tissue'] = "\"#86007D\""
for i in range(len(subtype_list)):
	#if subtype_list[i] == "Tissue": continue
	subtype2color_dict[subtype_list[i]] = color_list[i]
print (subtype_list, subtype2color_dict)

Others_subtype_prop_dict = defaultdict(lambda: defaultdict(lambda: []))
Top_pro_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: 0)))
Others_pro_dict = defaultdict(lambda: defaultdict(lambda: 0))
outf_3_name = re.sub(".txt", "_reshaped_merge_others.txt", proportion_inf_name.split('/')[-1])
outf_3 = open(outf_dir+'/'+outf_3_name, "w")
outf_3.write("Gene_ID\tTranscript_ID\tSample\tProportion\tSubtype\tSubtype_color\n")
with open(proportion_inf_name, "r") as inf_3:
	for index,line in enumerate(inf_3):
		arr = line.strip().split("\t")
		if index == 0: continue
		gene_ID = arr[1]
		trans_ID = arr[0]
		if trans_ID in top_exp_trans_dict[gene_ID]:
			for i in range(column,len(arr)):
				this_sample = sample_list_2[i-column]
				Top_pro_dict[gene_ID][this_sample][trans_ID] = float(arr[i])
		else:
			for i in range(column,len(arr)):
				this_sample = sample_list_2[i-column]
				Others_pro_dict[gene_ID][this_sample] += float(arr[i])

for each_gene in Top_pro_dict:
	for each_sample in sample_list_defined:
		this_subtype = cell2subtype_dict[each_sample]
		for each_trans in top_exp_trans_dict[each_gene]:
			this_pro = 0
			if each_trans in Top_pro_dict[each_gene][each_sample]:
				this_pro = Top_pro_dict[each_gene][each_sample][each_trans]
			outf_3.write(each_gene+"\t"+each_trans+"\t"+each_sample+"\t"+str(this_pro)+"\t"+this_subtype+"\t"+subtype2color_dict[this_subtype]+"\n")
for each_gene in Others_pro_dict:
	for each_sample in sample_list_defined:
		this_subtype = cell2subtype_dict[each_sample]
		this_pro = float(100-np.sum(list(map(float, Top_pro_dict[each_gene][each_sample].values()))))
		if this_pro < 0: this_pro = 0.0
		Others_subtype_prop_dict[each_gene][this_subtype].append(this_pro)
		#this_pro = Others_pro_dict[each_gene][each_sample]
		outf_3.write(each_gene+"\tOthers\t"+each_sample+"\t"+str(this_pro)+"\t"+this_subtype+"\t"+subtype2color_dict[this_subtype]+"\n")
outf_3.close()

### generate black gene list based on contribution from others isoforms ###
outf_3_other_name = re.sub(".txt", "_only_focus_others.txt", proportion_inf_name.split('/')[-1])
outf_3_other = open(outf_dir+'/'+outf_3_other_name,"w")
outf_3_other.write("Gene_ID\tSubtype\tMean_others_proportion\tMedian_others_proportion\n")
for each_gene in Others_subtype_prop_dict:
	for each_subtype in subtype_list:
		other_mean = np.nanmean(Others_subtype_prop_dict[each_gene][each_subtype])
		other_median = np.nanmedian(Others_subtype_prop_dict[each_gene][each_subtype])
		outf_3_other.write(each_gene+'\t'+each_subtype+'\t'+str(other_mean)+'\t'+str(other_median)+'\n')
outf_3_other.close()

### exp ###
exp_column = 3
sample_list_exp = []
Others_exp_dict = defaultdict(lambda: defaultdict(lambda: 0))
Top_exp_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: 0)))
outf_4_name = re.sub(".esp|.txt", "_reshaped_merge_others.txt", CPM_inf_name.split('/')[-1])
outf_4 = open(outf_dir+'/'+outf_4_name, "w")
outf_4.write("Gene_ID\tTranscript_ID\tSample\tCPM\tSubtype\tSubtype_color\n")
with open(CPM_inf_name, "r") as inf_4:
	for index,line in enumerate(inf_4):
		arr = line.strip().split("\t")
		if index == 0: 
			sample_list_exp = arr[exp_column:len(arr)]
			continue
		if len(arr[2].split(",")) > 1: continue # mapped to multiple gene
		if arr[2] == "NA": continue
		if re.findall("_PAR_Y", arr[2]):
			gene_ID = arr[2].split(".")[0]+"-PAR-Y"
		else:
			gene_ID = arr[2].split(".")[0]
		if re.findall("_PAR_Y", arr[0]):
			trans_ID = arr[0].split(".")[0]+"-PAR-Y"
		else:
			trans_ID = arr[0].split(".")[0]
		if gene_ID not in top_exp_trans_dict: continue
		if trans_ID in top_exp_trans_dict[gene_ID]:
			for i in range(exp_column,len(arr)):
				this_sample = sample_list_exp[i-exp_column]
				Top_exp_dict[gene_ID][this_sample][trans_ID] = float(arr[i])
		else:
			for i in range(exp_column,len(arr)):
				this_sample = sample_list_exp[i-exp_column]
				Others_exp_dict[gene_ID][this_sample] += float(arr[i])
	
for each_gene in Top_exp_dict:
	for each_sample in sample_list_defined:
		this_subtype = cell2subtype_dict[each_sample]
		for each_trans in top_exp_trans_dict[each_gene]:
			this_CPM = 0
			if each_trans in Top_exp_dict[each_gene][each_sample]:
				this_CPM = Top_exp_dict[each_gene][each_sample][each_trans]
			outf_4.write(each_gene+"\t"+each_trans+"\t"+each_sample+"\t"+str(this_CPM)+"\t"+this_subtype+"\t"+subtype2color_dict[this_subtype]+"\n")
for each_gene in Others_exp_dict:
	for each_sample in sample_list_defined:
		this_subtype = cell2subtype_dict[each_sample]
		this_CPM = Others_exp_dict[each_gene][each_sample]
		outf_4.write(each_gene+"\tOthers\t"+each_sample+"\t"+str(this_CPM)+"\t"+this_subtype+"\t"+subtype2color_dict[this_subtype]+"\n")
outf_4.close()
