####################################################
#Goal:      Translate transcripts into proteins
#Author:    Yang Xu
#E-mail:    yangax@pennmedicine.upenn.edu
####################################################

python /mnt/isilon/xing_lab/aspera/xuy/snakemake_ESPRESSO_reference/Translation_scripts/Translation.py \
--mode long-read \
--trans_gtf /home/xuy2/scratch/snakemake_1.2.2_Target_1115_BRCA_IMPACT_mix_test_batch3_rescue_FC1_8_final/espresso_out_combined/work_dir/work_dir_40_cell_lines/samples_N2_R0_updated.gtf \
--abundance_inf /home/xuy2/scratch/snakemake_1.2.2_Target_1115_BRCA_IMPACT_mix_test_batch3_rescue_FC1_8_final/espresso_out_combined/work_dir/work_dir_40_cell_lines/samples_N2_R0_abundance.esp \
--genome_version GRCH37 \
--ref_gtf /home/xuy2/scratch/snakemake_1.2.2_Target_1115_BRCA_IMPACT_mix_test_batch3_rescue_FC1_8_final/references/gencode.v34lift37.annotation.gtf \
--out_dir ./ \
--out_file BRCA
