import os,re,sys
from collections import defaultdict


stop_flag_list = ['five_prime_utr','stop_codon','transcript']
#gtf_inf = open('/home/xuy2/xuy2/Protein_Annotation_tool/Genome_Fasta_GTF/gencode.v34lift37.annotation.gtf')
#outf = open('1.5_transcript_annotation.txt','w')
gtf_inf = open(sys.argv[1],'r')
#gtf_inf_dire = '/'.join(os.path.abspath(gtf_inf).split('/')[0:-1])
outf_dir = sys.argv[2]
outf = open(outf_dir+'/Gencode_converted_CDS_annotation.txt','w')

coding_flag = 0
CDS_list = []
last_ENST_ID = ''
last_ENST_type = ''
last_ENSG_ID = ''
last_gene_name = ''
last_tag = ''
last_strand = ''
gene_range_dict = defaultdict()

#protein_coding_group = ['protein_coding','nonsense_mediated_decay','non_stop_decay']
protein_coding_group = ['protein_coding','nonsense_mediated_decay','non_stop_decay','retained_intron','processed_transcript']
for line in gtf_inf:
	arr = line.strip().split('\t')
	if line.startswith('#'): continue
	if arr[2] == 'gene':
		raw_ENSG_ID = re.findall('gene_id \"(.+?)\";',arr[8])[0]
		if re.findall('_PAR_Y', raw_ENSG_ID):
			ENSG_ID = raw_ENSG_ID.split('.')[0]+'-PAR-Y'
		else:
			ENSG_ID = raw_ENSG_ID.split('.')[0]
		gene_range = arr[0]+'#'+arr[6]+'#'+str(arr[3])+'#'+str(arr[4])
		gene_range_dict[ENSG_ID] = gene_range
		continue
	elif arr[2]=='transcript':
		if coding_flag == 1: #last transcript has not been done
			if (last_strand == '-') and (len(CDS_list)>=2):
				first_start = int(CDS_list[0].split('#')[0])
				second_start = int(CDS_list[1].split('#')[0])
				if first_start > second_start:
					CDS_list = CDS_list[::-1]
			CDS_region = ';'.join(CDS_list)
			CDS_list = []
			coding_flag = 0
			outf.write(last_ENST_ID+'\t'+last_ENST_type+'\t'+last_tag+'\t'+last_ENSG_ID+'\t'+last_gene_name+'\t'+gene_range_dict[last_ENSG_ID]+'\t'+CDS_region+'\n')
		# new transcript
		raw_ENST_ID = re.findall('transcript_id \"(.+?)\";',arr[8])[0]
		if re.findall('_PAR_Y', raw_ENST_ID):
			ENST_ID = raw_ENST_ID.split('.')[0]+'-PAR-Y'
		else:
			ENST_ID = raw_ENST_ID.split('.')[0]
		if re.findall('transcript_type',arr[8]):
			ENST_type = re.findall('transcript_type \"(\w+)\";',arr[8])[0]
		elif re.findall('transcript_biotype',arr[8]):
			ENST_type = re.findall('transcript_biotype \"(\w+)\";',arr[8])[0]
		raw_ENSG_ID = re.findall('gene_id \"(.+?)\";',arr[8])[0]
		if re.findall('_PAR_Y', raw_ENSG_ID):
			ENSG_ID = raw_ENSG_ID.split('.')[0]+'-PAR-Y'
		else:
			ENSG_ID = raw_ENSG_ID.split('.')[0]
		gene_name = re.findall('gene_name \"(.+?)\";',arr[8])[0]
		#tag = 'non_basic'
		#if re.findall('tag \"basic\";',arr[8]):
		#	tag = 'basic'
		tag = ';'.join(re.findall('tag \"(.+?)\";',arr[8]))
		strand = arr[6]
		if ENST_type in protein_coding_group:
			coding_flag = 1
			last_ENST_ID = ENST_ID
			last_ENST_type = ENST_type
			last_ENSG_ID = ENSG_ID
			last_gene_name = gene_name
			last_tag = tag
			last_strand = strand
			continue
		else:
			outf.write(ENST_ID+'\t'+ENST_type+'\t'+tag+'\t'+ENSG_ID+'\t'+gene_name+'\t'+gene_range_dict[ENSG_ID]+'\n')
	elif (coding_flag==1) and (arr[2]=='CDS'):  #only show CDS resgion
	#elif (coding_flag==1) and (arr[2]=='exon'): #show whole transcript
		CDS_exon = str(arr[3])+'#'+str(arr[4])
		CDS_list.append(CDS_exon)
gtf_inf.close()
outf.close()
