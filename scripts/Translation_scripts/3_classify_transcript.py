import os,re,sys
from collections import defaultdict

inf_name = sys.argv[1] # '2_CD44_hg19_protein.txt'
outf_dir = sys.argv[2]
ref_gtf_name = sys.argv[3]
gtf_name = sys.argv[4]
key_word = inf_name.split('/')[-1].split('_')[1]

dir_path = os.path.dirname(os.path.realpath(__file__))

#out_pre = inf_name.split('_hg19')[1]
numbering_dict = defaultdict()
outf_ID = open(outf_dir+'/3_'+key_word+'_detailed_match_ID.txt','w')
outf_fasta = open(outf_dir+'/3_'+key_word+'_PC.fasta','w')
outf_NMD_fasta = open(outf_dir+'/3_'+key_word+'_NMD.fasta','w')
outf_PTC_fasta = open(outf_dir+'/3_'+key_word+'_PC_PTC.fasta','w')
outf_ID.write('Transcript_ID\tGene_ID\tGene_name\tChr\tStrand\tCDS\tUniProt_ID\tStatus\n')

ENST2uni_dict = defaultdict()
ENSG2name_dict = defaultdict()
uni2ENSG_dict = defaultdict()
## load uniprot information ##
with open('%s/files/UP000005640_9606.idmapping_coverted_2.txt' % dir_path,'r') as pro_inf:
	for line in pro_inf:
		arr = line.strip().split('\t')
		Uni_ID = arr[0]
		Gene_ID = arr[2]
		Gene_name = arr[3]
		if len(arr) < 5: continue
		Trans_ID_list = arr[4].split(';')
		uni2ENSG_dict[Uni_ID] = Gene_ID
		if (Gene_ID == '-') or (Gene_ID == ''): continue
		for each_trans in Trans_ID_list:
			ENST2uni_dict[each_trans] = Uni_ID

with open(ref_gtf_name, 'r') as ref_gtf_inf:
	for line in ref_gtf_inf:
		if line.startswith('#'): continue
		arr = line.strip().split('\t')
		if arr[2] != 'gene': continue
		raw_gene_ID = re.findall('gene_id \"(.+?)\";',arr[8])[0]
		if re.findall('_PAR_Y', raw_gene_ID):
			gene_ID = raw_gene_ID.split('.')[0]+'-PAR-Y'
		else:
			gene_ID = raw_gene_ID.split('.')[0]
		gene_name = re.findall('gene_name \"(.+?)\";',arr[8])[0]
		ENSG2name_dict[gene_ID] = gene_name


with open(gtf_name, 'r') as gtf_inf:
	for line in gtf_inf:
		if line.startswith('#'): continue
		arr = line.strip().split('\t')
		if arr[2] != 'gene': continue
		if re.findall('ref_gene_id', arr[8]):
			raw_gene_ID = re.findall('ref_gene_id \"(.+?)\";',arr[8])[0]
		else:
			raw_gene_ID = re.findall('gene_id \"(.+?)\";',arr[8])[0]
		if re.findall('STRG', raw_gene_ID):
			gene_ID = raw_gene_ID
		elif re.findall('_PAR_Y', raw_gene_ID):
			gene_ID = raw_gene_ID.split('.')[0]+'-PAR-Y'
		else:
			gene_ID = raw_gene_ID.split('.')[0]
		if gene_ID not in ENSG2name_dict:
			if re.findall('ref_gene_name', arr[8]):
				gene_name = re.findall('ref_gene_name \"(.+?)\"',arr[8])[0]
			else:
				gene_name = re.findall('gene_name \"(.+?)\"',arr[8])[0]
			ENSG2name_dict[gene_ID] = gene_name

flag = ''
rename_ID = ''
ID_remove_redundancy_dict = defaultdict() 

inf = open(inf_name,'r')
for line in inf:
	line = line.strip()
	if line.startswith('>'):
		# remove redundancy
		if line.lstrip('>') in ID_remove_redundancy_dict:
			continue
		else:
			ID_remove_redundancy_dict[line.lstrip('>')] = 'yes'
		arr = line.strip().lstrip('>').split('_')
		status = arr[-1]
		chrom = arr[0]
		strand = arr[1]
		ENST_ID = arr[5]
		ENSG_ID = arr[6]
		CDS_region = arr[7]
		if ENST_ID in ENST2uni_dict:
			Uni_ID =  ENST2uni_dict[ENST_ID]
			protein_description = ENSG_ID
		else:
			Uni_ID = ENST_ID
			protein_description = ENSG_ID
		if ENSG_ID in ENSG2name_dict:
			gene_name = ENSG2name_dict[ENSG_ID]
		else:
			gene_name = '-'
		protein_name = ENST_ID+'_HUMAN'
		outf_ID.write(ENST_ID+'\t'+ENSG_ID+'\t'+gene_name+'\t'+chrom+'\t'+strand+'\t'+CDS_region+'\t'+Uni_ID+'\t'+status+'\n')
		#outf_ID.write(Re_name+'_'+str(numbering_dict[Re_name])+'\t'+line.lstrip('>')+'\n')
		rename_ID = '>db|'+Uni_ID+'|'+protein_name+' '+protein_description+' OS=Homo sapiens OX=9606 GN='+gene_name
		if re.findall('PC',status):
			flag = 'PC'
		elif re.findall('NMD',status):
			flag = 'NMD'
	else:
		if flag == 'PC':
			if (re.findall('\*',line.upper())) or (re.findall('XXXXX',line.upper())):
				outf_PTC_fasta.write(rename_ID+'\n'+line.upper()+'\n')
			else:
				outf_fasta.write(rename_ID+'\n'+line.upper()+'\n')
		elif flag == 'NMD':
			outf_NMD_fasta.write(rename_ID+'\n'+line.upper()+'\n')
		flag = ''
		rename_ID = ''
inf.close()
outf_ID.close()
outf_fasta.close()
outf_PTC_fasta.close()
outf_NMD_fasta.close()
