####################################################
#Goal:      Translate transcripts into proteins
#Author:    Yang Xu
#E-mail:    yangax@pennmedicine.upenn.edu
####################################################

python scripts/Translation.py \
--mode long-read \
--trans_gtf BRCA_Example/samples_N2_R0_updated.gtf \
--isoform_cpm_inf BRCA_Example/samples_N2_R0_abundance.esp \
--genome_version GRCH37 \
--ref_gtf scripts/files/gencode.v34lift37.annotation.gtf \
--outf_dir BRCA_Example \
--out_file BRCA
