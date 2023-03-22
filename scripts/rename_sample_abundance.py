import re,sys,os
from collections import defaultdict
import numpy as np
import pandas as pd

ID2sample_dict = defaultdict()
sample2subtype_dict = defaultdict()
with open("/home/xuy2/scratch/snakemake_1.2.2_Target_1115_BRCA_IMPACT_mix_test_batch3_rescue_FC1_8_final/BRCA_cell_lines.txt","r") as ref:
	for index, line in enumerate(ref):
		if index == 0: continue
		arr = line.strip().split('\t')
		ID2sample_dict[arr[0]] = arr[1]
		sample2subtype_dict[arr[1]] = arr[2]

matrix = []
with open("./samples_N2_R0_abundance_original.esp", 'r') as inf:
	for index, line in enumerate(inf):
		arr = line.strip().split("\t")
		if index == 0:
			this_list = arr[0:3]
			for i in range(3,len(arr)):
				code_ID = arr[i].split('_')[1]
				sample = ID2sample_dict[code_ID]+'_'+arr[i].split('_')[2]
				this_list.append(sample)
			matrix.append(this_list)
		else:
			matrix.append(arr)

pd_matrix = pd.DataFrame(matrix)
pd_matrix = pd_matrix.rename(columns=pd_matrix.iloc[0]).drop(pd_matrix.index[0]).reset_index(drop=True)

sorted_sample_list = ["transcript_ID","transcript_name","gene_ID","BT-549_1","BT-549_2","HCC1395_1","HCC1395_2","Hs-578T_1","Hs-578T_2","MDA-MB-157_1","MDA-MB-157_2","MDA-MB-231_1","MDA-MB-231_2","MDA-MB-436_1","MDA-MB-436_2","DU4475_1","DU4475_2","BT-20_1","BT-20_2","HCC38_1","HCC38_2","HCC70_1","HCC70_2","HCC1187_1","HCC1187_2","HCC1569_1","HCC1569_2","HCC1599_1","HCC1599_2","HCC1806_1","HCC1806_2","HCC1937_1","HCC1937_2","HCC1954_1","HCC1954_2","HCC2157_1","HCC2157_2","MDA-MB-468_1","MDA-MB-468_2","AU-565_1","AU-565_2","BT-474_1","BT-474_2","HCC202_1","HCC202_2","HCC1419_1","HCC1419_2","MDA-MB-175-VII_1","MDA-MB-175-VII_2","MDA-MB-361_1","MDA-MB-361_2","MDA-MB-415_1","MDA-MB-415_2","MDA-MB-453_1","MDA-MB-453_2","SK-BR-3_1","SK-BR-3_2","UACC-893_1","UACC-893_2","ZR-75-30_1","ZR-75-30_2","HCC2218_1","HCC2218_2","MDA-kb2_1","MDA-kb2_2","BT-483_1","BT-483_2","CAMA-1_1","CAMA-1_2","HCC1428_1","HCC1428_2","HCC1500_1","HCC1500_2","MCF7_1","MCF7_2","MDA-MB-134-VI_1","MDA-MB-134-VI_2","T47D_1","T47D_2","UACC-812_1","UACC-812_2","ZR-75-1_1","ZR-75-1_2"]
pd_matrix = pd_matrix.reindex(columns = sorted_sample_list)
np.savetxt("./samples_N2_R0_abundance.esp", pd_matrix, fmt='%s', delimiter='\t', newline='\n', header='\t'.join(sorted_sample_list))

