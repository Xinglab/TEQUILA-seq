import os,re,sys

gene_name = sys.argv[1]
target_trans_ID = sys.argv[2]
object_dir = sys.argv[3]
name_col = int(sys.argv[4])
abun_CPM_original = sys.argv[5]
sample_bedgraph = sys.argv[6]
file_dir = os.path.dirname(os.path.realpath(__file__))

## attention: need to change the column index
## normal format is Group_GeneSymbol_TransID.bed
## Therefore, the parameter should be 2 in command_3, and 1 in command_5

#### 1. generate bedGraph files for each transcript #######
command_1 = "python %s/Structure_trans_script/1_create_target_gene_bedgraph_sh.py %s %s 0 %s %s" % (file_dir, abun_CPM_original, gene_name, object_dir, sample_bedgraph)
#print (command_1)
os.system(command_1)

command_2 = "python %s/Structure_trans_script/Figure_exp_data.py %s/Figure_%s_bed_list.txt %s %s %s %s" % (file_dir, object_dir, gene_name, gene_name, name_col, abun_CPM_original, target_trans_ID)
#print (command_2)
os.system(command_2)

command_3 = "python %s/Isoform_structure_plot_single_color.py %s/Figure_%s_bed_list_sorted_%s_rank.txt %s %s ESPRESSO" % (file_dir, object_dir, gene_name, target_trans_ID, gene_name, name_col-1)
#print (command_3)
os.system(command_3)



