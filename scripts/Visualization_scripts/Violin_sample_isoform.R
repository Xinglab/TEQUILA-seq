#!/usr/bin/env Rscript
####################################################
#Goal:      Violin plot of isoform proportion for each subtype of cancer
#Author:    Yang Xu
#E-mail:    yangax@pennmedicine.upenn.edu
####################################################
args = commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(RColorBrewer)

inf_name = '/mnt/isilon/xing_lab/aspera/xuy/snakemake_1.2.2_Target_1229_BRCA_IMPACT_mix_all_45/samples_N2_R0_abundance_CPM_proportion_#_target_filtered_change_name_reshaped.txt'

fig_data <- read.table(inf_name, header = T,sep='\t')
target_gene = args[1]
target_gene_name = args[2]
object_path = args[3]
outf_name = paste(object_path,'/Violin_sample_isoform_',target_gene_name,'_',target_gene,'.png',sep='')
outf_pdf_name = paste(object_path,'/Violin_sample_isoform_',target_gene_name,'_',target_gene,'.pdf',sep='')
#target_gene = 'ENSG00000088305'

gene_data <- fig_data[fig_data$Transcript_ID == target_gene,]
gene_data$Subtype <- factor(gene_data$Subtype, levels=c('Luminal','HER2_Amp','Basal_A','Basal_B'))

font_size = 16
f_fig <- ggplot(gene_data, aes(x=Subtype, y=Proportion, color=Subtype, group=Subtype)) +scale_color_manual(name="Subtype",values=c("#0000F9","#02a620","#e68e19","#FF0018"),labels=c("Luminal","HER2 enriched","Basal A","Basal B"),breaks=c("Luminal","HER2_Amp","Basal_A","Basal_B")) 
#f_fig <- f_fig + scale_x_discrete(labels=c("Basal A","Basal B","HER2 Amp","Luminal","Non-tumorigenic"),breaks=c('Basal_A','Basal_B','HER2_Amp','Luminal','Non-tumorigenic'))
f_fig <- f_fig + scale_x_discrete(labels=c("Luminal","HER2\nenriched","Basal A","Basal B"),breaks=c('Luminal','HER2_Amp','Basal_A','Basal_B'))
f_fig <- f_fig + geom_violin() + geom_jitter(shape=16, position=position_jitter(0.2)) + geom_boxplot(width=0.1)
#f_fig <- f_fig + geom_text(aes(label=log), size=3, position = position_dodge(width = 0.5), hjust=-0.5)
f_fig <- f_fig + theme(panel.grid.major = element_line(color = "white",size=0.01),panel.grid.minor = element_line(color = "white",size=0.01), axis.title.y=element_text(size = font_size+2), axis.text.y=element_text(size = font_size+2), axis.title.x=element_blank(), axis.text.x=element_text(size = font_size-0.5, colour = 'black'), plot.title = element_text(hjust = 0.5, size=font_size+1), legend.position = 'none', legend.title = element_text(size = font_size), legend.text = element_text(size = font_size-0.5), panel.border = element_rect(fill=NA, size=0.3, colour = "black"), panel.background = element_blank(), plot.margin = unit(c(5,5,5,5), "pt"))
f_fig <- f_fig + labs(x="Subtypes", y="Isoform proportion (%)", title=paste(target_gene_name,': ',target_gene,sep=''))
f_fig <- f_fig + scale_y_continuous(limits=c(0,100))
#f_fig <- f_fig + scale_y_continuous(labels = scales::number_format(accuracy = 0.1), limits=c(0,15))
f_fig <- f_fig + guides(fill = guide_legend(reverse = FALSE))
f_fig


#png(file=outf_name, width=5, height=4,units='in',res=300)
png(file=outf_name, width=4.2, height=4,units='in',res=300)
f_fig
dev.off()
#ggsave(outf_pdf_name, plot=f_fig, width=5, height=4, dpi=300)
ggsave(outf_pdf_name, plot=f_fig, width=4.2, height=4, dpi=300)
