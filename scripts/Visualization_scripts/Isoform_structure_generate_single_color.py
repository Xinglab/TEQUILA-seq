import os,re,sys

gene_name = sys.argv[1]
target_trans_ID = sys.argv[2]
object_dir = sys.argv[3]
name_col = int(sys.argv[4])
abun_CPM_original = sys.argv[5]
sample_bedgraph = sys.argv[6]

## attention: need to change the column index
## normal format is Group_GeneSymbol_TransID.bed
## Therefore, the parameter should be 2 in command_3, and 1 in command_5

command_1 = "python ./Visualization_scripts/Create_target_gene_bedgraph_sh.py %s %s 0 %s %s" % (abun_CPM_original, gene_name, object_dir, sample_bedgraph)
#print (command_1)
os.system(command_1)

command_2 = "sh %s/BedGraph_TargetGene_%s.sh" % (object_dir, gene_name)
#print (command_2)
os.system(command_2)

command_3 = "python ./Visualization_scripts/Display_exp_data.py %s/Display_%s_bed_list.txt %s %s %s" % (object_dir, gene_name, gene_name, name_col, abun_CPM_original)
#print (command_3)
os.system(command_3)

command_4 = "python ./Visualization_scripts/rank_isoform_by_bar_plot_proportion.py %s/Display_%s_bed_list_sorted.txt %s/Proportion_sample_isoform_%s_%s.txt" % (object_dir, gene_name, object_dir, gene_name, target_trans_ID)
#print (command_4)
os.system(command_4)

command_5 = "python ./Visualization_scripts/Isoform_structure_plot_single_color.py %s/Display_%s_bed_list_sorted_%s_rank.txt %s %s ESPRESSO" % (object_dir, gene_name, target_trans_ID, gene_name, name_col-1)
#print (command_5)
os.system(command_5)



