#!/usr/bin/env Rscript
####################################################
#Goal:		CPM and proportion bar plot for sample-specific isoforms
#Author:	Yang Xu
#E-mail:	yangax@pennmedicine.upenn.edu
####################################################

args = commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(dplyr)

#inf_name = './samples_N2_R0_abundance_combined_with_tissue_CPM_proportion_more_than_1_change_name_reshaped_merge_others.txt'

target_gene = args[1]
target_gene_name = args[2]
target_trans = args[3]
object_path = args[4]
proportion_inf_name = args[5]
CPM_inf_name = args[6]
order_flag = args[7]
subtype_order <- strsplit(args[8],',')[[1]]
uniprot_inf_name = args[9] # '/mnt/isilon/xing_lab/aspera/xuy/snakemake_ESPRESSO_reference/pipeline_test/scripts_example_figure/Gencode_v39_canonical_isoform.txt'
basic_inf_name = args[10] #'/mnt/isilon/xing_lab/aspera/xuy/snakemake_ESPRESSO_reference/pipeline_test/scripts_example_figure/gencode.v34lift37.annotation_basic_trans.txt'

subtype_order
fig_data <- read.table(proportion_inf_name, header = T,sep='\t')
outf_name = paste(object_path,'/Bar_sample_isoform_',target_gene_name,'_',target_trans,'.png',sep='')
outf_pdf_name = paste(object_path,'/Bar_sample_isoform_',target_gene_name,'_',target_trans,'.pdf',sep='')

gene_data <- fig_data[fig_data$Gene_ID == target_gene,]
## sort transcript based on proportion ##
mean_pro_data <- gene_data %>% group_by(Transcript_ID) %>% summarize(mean_Proportion = mean(Proportion, na.rm=TRUE))
mean_pro_data <- mean_pro_data[order(mean_pro_data$mean_Proportion, decreasing=TRUE),]
isoform_order_list <- as.vector(mean_pro_data$Transcript_ID)
isoform_order_list
abundance_rank <- c(1:length(isoform_order_list))
isoform_order_table <- as.data.frame(cbind(isoform_order_list, abundance_rank))


## sort samples based on proportion of given isoform ##
gene_data$Subtype <- factor(gene_data$Subtype, levels = subtype_order)
#gene_data$Subtype <- factor(gene_data$Subtype, levels=unique(as.vector(gene_data$Subtype)))
gene_data <- gene_data[order(gene_data$Subtype),]

trans_data <- gene_data[gene_data$Transcript_ID == target_trans,]
if (order_flag == 'yes'){
	trans_data <- trans_data[order(trans_data$Proportion, decreasing=TRUE),]
}
gene_data$Sample <- factor(gene_data$Sample, levels=as.vector(trans_data$Sample))
sub_color <- as.vector(trans_data$Subtype_color)


### sort isoforms ####
#uniprot_inf_name = '/mnt/isilon/xing_lab/aspera/xuy/snakemake_1.2.2_Target_1229_BRCA_IMPACT_mix_all_45/Target_IMPACT_gene_list.Ensembl_Release_100.Canonical_UniProt.txt'
uniprot_inf <- read.table(uniprot_inf_name, header = T,sep='\t')
canonical_trans <- as.character(uniprot_inf[uniprot_inf$GeneSymbol==target_gene_name, 3])

basic_inf <- read.table(basic_inf_name, header = T,sep='\t')
if (grepl('lift37', basic_inf_name)){
	basic_trans <- as.character(basic_inf[basic_inf$Gene_name==target_gene_name, 2])
}else{
	basic_trans <- as.character(basic_inf[basic_inf$Gene_ID==target_gene, 3])
}


canonical_trans
basic_trans
#isoform_order_list

if(target_trans == canonical_trans){
	print ('Aberrant transcript is canonical transcript')
	new_cano <- ''
	for(each_iso in isoform_order_list){
		if (each_iso != canonical_trans){
			if (grepl(each_iso, basic_trans)){
				new_cano <- as.character(each_iso)
				break
			}
		}
	}
	if (new_cano == ''){
		print ('The new basic reference isoform cannot be found')
		if ('Others' %in% isoform_order_list){
			isoform_order_list_rest <- isoform_order_list[!isoform_order_list %in% c(target_trans,'Others')]
			isoform_order_list <- c(target_trans,isoform_order_list_rest,'Others')
		} else{
			isoform_order_list_rest <- isoform_order_list[!isoform_order_list %in% c(target_trans)]
			isoform_order_list <- c(target_trans,isoform_order_list_rest)
		}	
	} else{
		if ('Others' %in% isoform_order_list){
			isoform_order_list_rest <- isoform_order_list[!isoform_order_list %in% c(target_trans,new_cano,'Others')]
			isoform_order_list <- c(target_trans,new_cano,isoform_order_list_rest,'Others')
		} else{
			isoform_order_list_rest <- isoform_order_list[!isoform_order_list %in% c(target_trans,new_cano)]
			isoform_order_list <- c(target_trans,new_cano,isoform_order_list_rest)
		}
	}
} else{
	print ('Aberrant transcript is not canonical transcript')
	if ('Others' %in% isoform_order_list){
		isoform_order_list_rest <- isoform_order_list[!isoform_order_list %in% c(canonical_trans,target_trans,'Others')]
		isoform_order_list <- c(target_trans,canonical_trans,isoform_order_list_rest,'Others')
	} else{
		isoform_order_list_rest <- isoform_order_list[!isoform_order_list %in% c(canonical_trans,target_trans)]
		isoform_order_list <- c(target_trans,canonical_trans,isoform_order_list_rest)
	}
}

isoform_order_list
isoform_reorder_table <- isoform_order_table[match(isoform_order_list, isoform_order_table$isoform_order_list),]
isoform_reorder_table


###### modify the isoform ID ####
display_target_trans = unlist(strsplit(target_trans,'@'))[1]
outf_Trans_list_name = paste(object_path,'/Proportion_sample_isoform_',target_gene_name,'_',target_trans,'.txt',sep='')
write.table(isoform_reorder_table, file = outf_Trans_list_name, sep = "\t", quote = FALSE, row.names = F, col.names=F)
#######


gene_data$Transcript_ID <- factor(gene_data$Transcript_ID, levels=isoform_order_list)



#trans_color_original_list <- c('#F31E08','#2171b5','#4292C6','#6baed6','#9ecae1','#97a1a6')
trans_color_original_list <- c('#FF5F42','#003E7F','#0068AF','#5495E1','#A1C2E8','#D9D9D9')
trans_color_table <- as.data.frame(cbind(trans_color_original_list,c(1:length(trans_color_original_list))))
trans_reorder_list <- as.integer(as.vector(isoform_reorder_table[,2]))
trans_reorder_color_table <- trans_color_table[match(trans_reorder_list, trans_color_table[,2]),]
#trans_color_list <- as.vector(trans_reorder_color_table[,1])
trans_color_list <- trans_color_original_list 

group_color_original_list <- c("#FF0018","#e68e19","#02a620","#0000F9","#86007D")
group_color_list <- group_color_original_list[0:length(unique(gene_data$Subtype))]
#group_color_list <- c(group_color_list[-length(group_color_list)], "#86007D")


font_size = 10
f_fig <- ggplot(gene_data, aes(x=Sample, y=Proportion, fill=Transcript_ID, group=Transcript_ID))
f_fig <- f_fig + geom_point(aes(y=NA_integer_, color = Subtype),size=6,shape=15)
f_fig <- f_fig + scale_colour_manual(name="Group",values = group_color_list)
f_fig <- f_fig + geom_bar(stat="identity", position = position_stack(vjust = 1, reverse = FALSE), width=0.9, color="black", size=0.5)
f_fig <- f_fig + scale_fill_manual(values = trans_color_list)
f_fig <- f_fig + theme(panel.grid.major = element_line(color = "white",size=0.01),panel.grid.minor = element_line(color = "white",size=0.01), axis.title.y=element_text(size = font_size), axis.text.y=element_text(size = font_size+2), axis.title.x=element_text(size = font_size), axis.text.x=element_text(size = font_size+1, colour = sub_color, angle=90, vjust=0.5, hjust=1), plot.title = element_blank(), legend.position = 'none', legend.title = element_text(size = font_size-1, hjust=0.5), legend.text = element_text(size = font_size-1.5), panel.border = element_rect(fill=NA, size=0.3, colour = "black"), panel.background = element_blank(), legend.box="horizontal", legend.direction="vertical", plot.margin = unit(c(0.5,5,0.5,5), "pt"))
f_fig <- f_fig + labs(x="", y="Isoform proportion (%)", title="")
#f_fig <- f_fig + scale_y_continuous(label = unit_format(unit = "K",sep="",accuracy=0.1))
#f_fig <- f_fig + guides(fill = guide_legend(reverse = FALSE, ncol=1, label.theme = element_text(size = font_size-1),keywidth=1, keyheight=1, order=1)) + guides(color = guide_legend(reverse = FALSE, ncol=1, label.theme = element_text(size = font_size-1), keywidth=1, keyheight=1, order=2))

####### Exp bar plot  ##########
#inf_name_exp = './samples_N2_R0_abundance_combined_with_tissue_CPM_more_than_1_change_name_reshaped_merge_others.txt'

fig_data_exp <- read.table(CPM_inf_name, header = T,sep='\t')
outf_exp_name = paste(object_path,'/Bar_sample_isoform_',target_gene_name,'_',target_trans,'_exp.png',sep='')
outf_pdf_exp_name = paste(object_path,'/Bar_sample_isoform_',target_gene_name,'_',target_trans,'_exp.pdf',sep='')

gene_data_exp <- fig_data_exp[fig_data_exp$Gene_ID == target_gene,]

## sort samples based on proportion of given isoform (based on above file) ##
gene_data_exp$Sample <- factor(gene_data_exp$Sample, levels=as.vector(trans_data$Sample))
sub_color <- as.vector(trans_data$Subtype_color)

## sort transcript based on proportion (based on above file) ##
gene_data_exp$Transcript_ID <- factor(gene_data_exp$Transcript_ID, levels=isoform_order_list)

## plot ##
f_fig_exp <- ggplot(gene_data_exp, aes(x=Sample, y=CPM, fill=Transcript_ID, group=Transcript_ID))
f_fig_exp <- f_fig_exp + geom_point(aes(y=NA_integer_, color = Subtype),size=6,shape=15)
f_fig_exp <- f_fig_exp + scale_colour_manual(name="Group",values = group_color_list)
f_fig_exp <- f_fig_exp + geom_bar(stat="identity", position = position_stack(vjust = 1, reverse = TRUE), width=0.9, color="black", size=0.5)
f_fig_exp <- f_fig_exp + scale_fill_manual(values = trans_color_list)
f_fig_exp <- f_fig_exp + theme(panel.grid.major = element_line(color = "white",size=0.01),panel.grid.minor = element_line(color = "white",size=0.01), axis.title.y=element_text(size = font_size), axis.text.y=element_text(size = font_size+2), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.text.x=element_blank(), plot.title = element_blank(), legend.position ='none', legend.title = element_text(size = font_size-1), legend.text = element_text(size = font_size-2), panel.border = element_rect(fill=NA, size=0.3, colour ="black"), panel.background = element_blank(), legend.box="horizontal", legend.direction="horizontal", plot.margin = unit(c(5,5,0.5,5), "pt"))
f_fig_exp <- f_fig_exp + labs(x="", y="Abundance (CPM)", title="")
#f_fig_exp <- f_fig_exp + guides(fill = guide_legend(reverse = FALSE, nrow=2, label.theme = element_text(size = 5),keywidth=1, keyheight=1, order=1)) + guides(color = guide_legend(reverse = FALSE, nrow=2, label.theme = element_text(size = 5), keywidth=1, keyheight=1, order=2))


#### combined plots ####
legend_a <- get_legend(f_fig + guides(fill = guide_legend(reverse = FALSE, nrow=2, title.position='top', label.theme = element_text(size = font_size-1.5), keywidth=1, keyheight=1, order=1)) + guides(color = guide_legend(reverse = FALSE, nrow=2, title.position='top', label.theme = element_text(size = font_size-1.5), keywidth=1, keyheight=1, order=2)) + theme(legend.position='bottom'))
box_plot <- plot_grid(f_fig_exp, f_fig, align = "v", ncol=1, rel_heights = c(1,1.8))
combined_plot <- plot_grid(box_plot, legend_a, ncol=1, align='v', axis='t', rel_heights=c(1,0.2))


def_width = max(0.25*length(unique(gene_data$Sample)) + 1, 8)
png(file=outf_exp_name, width=def_width, height=5.8, units='in',res=1500)
combined_plot
dev.off()
ggsave(outf_pdf_exp_name, plot=combined_plot, width=def_width, height=5.8, dpi=300)

