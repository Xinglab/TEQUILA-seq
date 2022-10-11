import os,re,sys
from collections import defaultdict

CPM_inf_name = sys.argv[1]
target_gene_ID = sys.argv[2]
CPM_cut_off = 1.0
sample_bedgraph = sys.argv[5]

if len(sys.argv) > 3:
	CPM_cut_off = sys.argv[3]

key_word = os.path.abspath(CPM_inf_name).split('/')[-1].split('_')[0]
sample_list = []

CPM_inf = open(CPM_inf_name,'r')
for line_index,line in enumerate(CPM_inf):
	arr = line.strip().split('\t')
	if line_index == 0:
		#sample_list = arr[3:len(arr)]
		sample_list = arr[3:4]
		break
CPM_inf.close()

###### generate 1.* target genes' bedGraph sh ####
object_path = '/'.join(os.path.abspath(CPM_inf_name).split('/')[0:-1])+'/statistics'
if len(sys.argv) > 4:
	object_path = sys.argv[4]

outf_name = object_path+'/BedGraph_TargetGene_'+target_gene_ID+'.sh'
outf = open(outf_name,'w')
dire_path = '/'.join(os.path.abspath(__file__).split('/')[0:-1])

for each_sample in sample_list:
	cmd = 'python '+dire_path+'/1.5_generate_bedgraph_for_each_gene.py '+each_sample+' '+target_gene_ID+' '+CPM_inf_name+' '+CPM_cut_off+' '+object_path+' '+sample_bedgraph+'\n'
	outf.write(cmd)
	cmd_2 = 'ls --color=never '+object_path+'/target_genes/*'+target_gene_ID+'*.bed > '+object_path+'/Figure_'+target_gene_ID+'_bed_list.txt\n'
	outf.write(cmd_2)
outf.close()


###### run those bash scripts #######
command_2 = "sh %s/BedGraph_TargetGene_%s.sh" % (object_path, target_gene_ID)
os.system(command_2)
