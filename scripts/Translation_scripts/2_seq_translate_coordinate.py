import os,sys,re
from collections import defaultdict
from pygr import seqdb
from Bio.Seq import Seq
from Bio.Alphabet import generic_rna
import configargparse

def parse_args():
	parser = configargparse.ArgParser(description='Translate transcripts into proteins')
	parser.add('-i', '--trans_inf', dest='trans_inf', type=str, help='input file of transcript structure', required=True)
	parser.add('-a', '--abundance_inf', dest='abundance_inf', type=str, help='input file of transcript abundance, e.g. samples_N2_R0_abundance.esp', required=True)
	parser.add('-g', '--gtf_inf', dest='gtf_inf', type=str, help='input gtf file of transcript, e.g. samples_N2_R0_updated.gtf', required=True)
	parser.add('-o', '--outf_name', dest='outf_name', type=str, help='name prefix of output files', required=True)
	parser.add('-f', '--genome', dest='genome', type=str, help='human reference genome fasta', required=True)
	parser.add('-c', '--cds', dest='cds', type=str, help='CDS annotation from Gencode file', required=True)
	parser.add('-m', '--mode', dest='mode', type=str, help='short-read input or long-read input', required=True)
	args = parser.parse_args()
	return parser, args

def read_genome(fn):
	#the input of this function should be a FASTA file
	genome = seqdb.SequenceFileDB(fn)
	return genome

def fetch_seq(genome, chr, start, end, strand):
	seq = ""
	try:
		seq = genome[chr][start:end]
		if strand == "-":
			seq = -seq  #ywang: make it reverse complement
	except:
		raise Exception('pygr cannot fetch sequences')
	return str(seq).upper()

def codon_judge(start_list,stop_list):
	pair_codon = []	
	used_end_0 = -1 
	used_end_1 = -1
	used_end_2 = -1
	for each_start in start_list:
		for each_end in stop_list:
			distance = int(each_end) - int(each_start)
			if 0<distance <60 and distance%3==0:
				break # cannot start with this start_coordinate
			if distance >= 60 and distance%3==0:
				if int(each_start)%3==0 and int(each_start)>used_end_0:
					res = str(each_start)+'#'+str(each_end)+'#0'
					pair_codon.append(res)
					used_end_0 = int(each_end)
					break
				elif int(each_start)%3==1 and int(each_start)>used_end_1:
					res = str(each_start)+'#'+str(each_end)+'#1'
					pair_codon.append(res)
					used_end_1 = int(each_end)
					break
				elif int(each_start)%3==2 and int(each_start)>used_end_2:
					res = str(each_start)+'#'+str(each_end)+'#2'
					pair_codon.append(res)
					used_end_2 = int(each_end)
					break
	return pair_codon

def lower_ss(pep_s,position_list):
	if len(position_list)==0:
		return pep_s
	else:
		for each_pos in position_list:
			aa_pos = int(each_pos)/3
			if aa_pos < len(pep_s):
				temp_pep_s = pep_s[0:aa_pos] + pep_s[aa_pos].lower() + pep_s[aa_pos+1:]
				pep_s = temp_pep_s
		return pep_s	

def determine_equal_CDS(anno_CDS_1,find_CDS_2):
	anno_list = anno_CDS_1.split(';')
	find_list = find_CDS_2.split(';')
	ans = False
	if len(anno_list) == len(find_list): 
		if len(anno_list) >= 2:
			find_list[0] = anno_list[0].split('#')[0]+'#'+find_list[0].split('#')[1] #ignore the start and end position of found CDS
			find_list[-1] = find_list[-1].split('#')[0]+'#'+anno_list[-1].split('#')[1]
			new_find_CDS = ';'.join(find_list)
			if anno_CDS_1 == new_find_CDS:
				ans = True
		else:
			ans = True #if only one CDS region, then consider it is equivalent
	return ans
		
def genome_position(s_trans,e_trans,exon_link):
	split_region = ''
	if re.findall(';',exon_link):
		exon_list = exon_link.split(';')
		pre_pos = 0 
		start_record_tag = 0
		end_record_tag = 0
		strand = exon_list[0].split('#')[2]
		if strand == '+':
			for each_exon in exon_list:
				temp_arr = each_exon.split('#')
				temp_dist = int(temp_arr[4]) - int(temp_arr[3]) + 1
				this_s_trans = int(s_trans) - pre_pos
				this_e_trans = int(e_trans) - pre_pos
				#print pre_pos,temp_dist,s_trans,this_s_trans,e_trans,this_e_trans
				if pre_pos <= int(s_trans) < pre_pos+temp_dist:  #in this region
					absolute_s = int(temp_arr[3]) + int(this_s_trans)
					start_record_tag = 1
					#print 'ss', absolute_s
				if pre_pos <= int(e_trans) <= pre_pos+temp_dist:
					absolute_e = int(temp_arr[3]) + int(this_e_trans) - 1
					end_record_tag = 1
					#print 'ee', absolute_e
				pre_pos = pre_pos + temp_dist
				
				# split_record_tag: 1: find start position; 2: find stop position; 0.5: in the middle part
				if (start_record_tag==1) and (end_record_tag==0):  
					split_region_str = str(absolute_s)+'#'+str(temp_arr[4])
					start_record_tag = 2
				elif (start_record_tag==1) and (end_record_tag == 1):
					split_region_str = str(absolute_s)+'#'+str(absolute_e)
					start_record_tag = 2
					end_record_tag = 2
				elif (start_record_tag==2) and (end_record_tag == 0):
					split_region_str = str(temp_arr[3])+'#'+str(temp_arr[4])
				elif (start_record_tag==2) and (end_record_tag == 1):
					split_region_str = str(temp_arr[3])+'#'+str(absolute_e)
					end_record_tag = 2

				if start_record_tag == 2:
					if split_region=='':
						split_region = split_region_str
					else:
						split_region = split_region +';'+ split_region_str
					if end_record_tag==2:
						break
		elif strand == '-':
			first_start = exon_list[0].split('#')[3]
			second_start = exon_list[1].split('#')[3]
			if int(first_start) < int(second_start):
				exon_list_neg = exon_list[::-1]
			else:
				exon_list_neg = exon_list
			for each_exon in exon_list_neg:
				temp_arr = each_exon.split('#')
				temp_dist = int(temp_arr[4]) - int(temp_arr[3]) + 1
				this_s_trans = int(s_trans) - pre_pos  #distance from the end of this exon; negative strand
				this_e_trans = int(e_trans) - pre_pos
				if pre_pos <= int(s_trans) < pre_pos+temp_dist:  #in this region
					absolute_s = int(temp_arr[4]) - int(this_s_trans)
					start_record_tag = 1
				if pre_pos <= int(e_trans) <= pre_pos+temp_dist:
					absolute_e = int(temp_arr[4]) - int(this_e_trans) + 1 #compatible with gtf coordinate
					end_record_tag = 1
				pre_pos = pre_pos + temp_dist
				# split_record_tag: 1: find start position; 2: find stop position; 0.5: in the middle part
				if (start_record_tag==1) and (end_record_tag==0):
					split_region_str = str(temp_arr[3])+'#'+str(absolute_s)
					start_record_tag = 2
				elif (start_record_tag==1) and (end_record_tag == 1):
					split_region_str = str(absolute_e)+'#'+str(absolute_s)
					start_record_tag = 2
					end_record_tag = 2
				elif (start_record_tag==2) and (end_record_tag == 0):
					split_region_str = str(temp_arr[3])+'#'+str(temp_arr[4])
				elif (start_record_tag==2) and (end_record_tag == 1):
					split_region_str = str(absolute_e)+'#'+str(temp_arr[4])
					end_record_tag = 2

				if start_record_tag == 2:
					if split_region=='':
						split_region = split_region_str
					else:
						split_region = split_region +';'+ split_region_str
					if end_record_tag==2:
						split_region = ';'.join(split_region.split(';')[::-1]) #from 3#4;1#2 to 1#2;3#4
						break
	else:
		strand = exon_link.split('#')[2]
		if strand == '+':
			start_position = exon_link.split('#')[3]
			end_position = exon_link.split('#')[4]
			absolute_s = int(start_position) + int(s_trans)
			absolute_e = int(start_position) + int(e_trans) - 1
			split_region = str(absolute_s)+'#'+str(absolute_e)
		elif strand == '-':
			start_position = exon_link.split('#')[3]
			end_position = exon_link.split('#')[4]
			absolute_s = int(end_position) - int(s_trans)
			absolute_e = int(end_position) - int(e_trans) + 1
			split_region = str(absolute_e)+'#'+str(absolute_s)
	return [absolute_s,absolute_e,split_region]

def check_NMD(this_trans,full_trans,strand):
	#print this_trans,full_trans,strand
	exon_num = len(full_trans.split(';'))
	# condition 1
	if exon_num < 2:
		return False
	# condition 2
	if strand == '+':
		this_stop_codon = int(this_trans.split(';')[-1].split('#')[-1])
		ORF_last_exon_start = int(this_trans.split(';')[-1].split('#')[-2])
		dist = 0   ### last_EJC - used_stop_codon
		for each_exon_index, each_exon in enumerate(full_trans.split(';')[::-1]):
			each_exon_start = int(each_exon.split('#')[-2])
			each_exon_end = int(each_exon.split('#')[-1])
			if each_exon_index == 0:  ## the last exon
				if each_exon_start == ORF_last_exon_start: ### stop codon is in the last exon
					return False
			else:
				if each_exon_start == ORF_last_exon_start:
					dist += (each_exon_end - this_stop_codon)
					break
				elif each_exon_start > ORF_last_exon_start:
					dist += (each_exon_end - each_exon_start + 1)
		if dist < 50: #No NMD
			return False
	else: #strand == '-'
		this_stop_codon = int(this_trans.split(';')[0].split('#')[-2])
		ORF_last_exon_start = int(this_trans.split(';')[0].split('#')[-1])
		dist = 0   ### last_EJC - used_stop_codon
		for each_exon_index, each_exon in enumerate(full_trans.split(';')):
			each_exon_end = int(each_exon.split('#')[-2])
			each_exon_start = int(each_exon.split('#')[-1])
			if each_exon_index == 0:  ## the last exon
				if each_exon_start == ORF_last_exon_start: ### stop codon is in the last exon
					return False
			else:
				if each_exon_start == ORF_last_exon_start:
					dist += (this_stop_codon - each_exon_end)
					break
				elif each_exon_start < ORF_last_exon_start:
					dist += (each_exon_start - each_exon_end + 1)
		if dist < 50:
			return False
	# condition 3
	ORF_len = 0
	exon_list = this_trans.split(';')
	for each_exon in exon_list:
		exon_len = int(each_exon.split('#')[-1]) - int(each_exon.split('#')[-2]) + 1
		ORF_len += exon_len
	#print 'exon_len',ORF_len
	if ORF_len < 200: # or 200, too short transcript doesn't have NMD
		return False
	# if pass all 3 tests
	return True


############ main script ##############
parser, args = parse_args()

genome = read_genome(args.genome)
input_name = args.trans_inf
out_name = args.outf_name
mode = args.mode
gtf = args.gtf_inf
inf = open(input_name,'r')
outf_transcript = open(out_name + '_transcript.txt','w')
outf_protein = open(out_name + '_protein.txt','w')
use_annotation = 'yes'

CDS_dict = defaultdict()
ENST2ENSG_dict = defaultdict()
#protein_coding_group = ['protein_coding','nonsense_mediated_decay','non_stop_decay']
transcript_type_dict = defaultdict()
if use_annotation == 'yes':              ##### loading gtf ########
	inf_gtf = open(args.cds)
	for line in inf_gtf:
		arr = line.strip().split('\t')
		transcript_type_dict[arr[0]] = arr[1]
		ENST2ENSG_dict[arr[0]] = arr[3]
		if (arr[1] == 'protein_coding') and (re.findall('basic',arr[2])):
			#if re.findall('PAR',arr[2]): continue ## skip PAR transcripts
			CDS_dict[arr[0]]=arr[-1]
		elif arr[1] == 'nonsense_mediated_decay':
			CDS_dict[arr[0]]=arr[-1]
	inf_gtf.close()

#### long-read #####
if mode == 'long-read':
	inf_exp_data = open(args.abundance_inf,'r')
	for index,line in enumerate(inf_exp_data):
		arr = line.strip().split('\t')
		if index == 0: continue
		if re.findall('_PAR_Y', arr[0]):
			ENST_ID = arr[0].split('.')[0]+'-PAR-Y'
		else:
			ENST_ID = arr[0].split('.')[0]
		if re.findall('_PAR_Y', arr[2]):
			ENSG_ID = arr[2].split('.')[0]+'-PAR-Y'
		else:
			ENSG_ID = arr[2].split('.')[0]
		if ENST_ID.startswith('ESPRESSO'):
			ENST2ENSG_dict[ENST_ID] = ENSG_ID
		else:
			if ENST_ID in ENST2ENSG_dict:
				if ENST2ENSG_dict[ENST_ID] != ENSG_ID:
					print ('Conflict gene id',ENST_ID,ENST2ENSG_dict[ENST_ID],ENSG_ID)
			else:
				ENST2ENSG_dict[ENST_ID] = ENSG_ID
				print ('not found in gtf',ENST_ID,ENSG_ID)
	inf_exp_data.close()

#### short-read ######
elif mode == 'short-read':
	inf_gtf_data = open(gtf,'r')
	for line in inf_gtf_data:
		if line.startswith('#'): continue
		arr = line.strip().split('\t')
		if re.findall('reference_id',arr[8]):
			raw_ENST_ID = re.findall('reference_id \"(.+?)\"',arr[8])[0]
			raw_ENSG_ID = re.findall('ref_gene_id \"(.+?)\"',arr[8])[0]
		else:	
			raw_ENST_ID = re.findall('transcript_id \"(.+?)\"',arr[8])[0]
			raw_ENSG_ID = re.findall('gene_id \"(.+?)\"',arr[8])[0]
		if re.findall('STRG', raw_ENST_ID):
			ENST_ID = raw_ENST_ID
		elif re.findall('_PAR_Y', raw_ENST_ID):
			ENST_ID = raw_ENST_ID.split('.')[0]+'-PAR-Y'
		else:
			ENST_ID = raw_ENST_ID.split('.')[0]
		if re.findall('STRG', raw_ENSG_ID):
			ENSG_ID = raw_ENSG_ID
		elif re.findall('_PAR_Y', raw_ENSG_ID):
			ENSG_ID = raw_ENSG_ID.split('.')[0]+'-PAR-Y'
		else:
			ENSG_ID = raw_ENSG_ID.split('.')[0]
		ENST2ENSG_dict[ENST_ID] = ENSG_ID
	inf_gtf_data.close()


###################################
print 'number of basic protein coding transcripts',len(CDS_dict)


###################################
for line in inf:
	arr = line.strip().split('\t')
	ENST_ID = ''
	transcript_link = arr[0]
	#####	
	if re.findall(';',transcript_link):
		new_arr = transcript_link.split(';')[0]
	else:	## single exon
		new_arr = transcript_link
	####
	this_exon = new_arr.split('#')
	if re.findall('/',this_exon[0]):
		raw_ENST_ID = this_exon[0].split('/')[0]  #ENST00000548656
	else:
		raw_ENST_ID = this_exon[0]
	if re.findall('STRG', raw_ENST_ID):
		ENST_ID = raw_ENST_ID
	elif re.findall('_PAR_Y', raw_ENST_ID):
		ENST_ID = raw_ENST_ID.split('.')[0]+'-PAR-Y'
	else:
		ENST_ID = raw_ENST_ID.split('.')[0]

	if re.findall('_', this_exon[1]): continue
	#if re.findall('_PAR_Y', this_exon[0]): continue ## remove PAR gene in Y chromsome, since we have same sequence in X.
	chrom = re.findall('(chr\w+)',this_exon[1])[0]
	if chrom == 'chrM' or re.findall('_',chrom):
		continue

	strand = this_exon[2]
	ENSG_ID = '-'
	if ENST_ID not in ENST2ENSG_dict: continue
	ENSG_ID = ENST2ENSG_dict[ENST_ID]
	#print chrom,int(this_exon[3])-1,int(this_exon[4]),strand
	mark_weird_isoform = 0

	if (use_annotation == 'yes') and (ENST_ID in CDS_dict):   #if this ENST is annotated, we use its annotated ORF
		# Obtain the full length sequence of transcripts
		final_seq = ''  #full mRNA sequence
		if re.findall(';', transcript_link):  # multiple exons in transcript
			each_exon_list = transcript_link.split(';')
			if strand == '-':
				first_start = each_exon_list[0].split('#')[3]
				second_start = each_exon_list[1].split('#')[3]
				if int(first_start) < int(second_start):
					each_exon_list = each_exon_list[::-1]
			for each_exon in each_exon_list:
				exon_detail = each_exon.split('#')
				start = int(exon_detail[3])
				end = int(exon_detail[4])
				if start -1 >= end:
					mark_weird_isoform = 1
					break
				temp_seq = fetch_seq(genome, chrom, start-1, end, strand)
				seq = temp_seq[0:-1] + temp_seq[-1].lower()
				final_seq += seq
			if mark_weird_isoform == 1:
				continue
			temp_final_seq = final_seq[0:-1] + final_seq[-1].upper()
			final_seq = temp_final_seq
		else: # only one exon
			if int(this_exon[3])-1 >= int(this_exon[4]):
				mark_weird_isoform = 1
				continue
			final_seq = fetch_seq(genome, chrom, int(this_exon[3])-1, int(this_exon[4]), strand)

		## use annotated CDS region
		Full_CDS_seq = ''
		anno_CDS_list = CDS_dict[ENST_ID].split(';')
		CDS_start = anno_CDS_list[0].split('#')[0]
		CDS_end = anno_CDS_list[-1].split('#')[1]
		if strand == '-':
			if len(anno_CDS_list) >= 2:
				first_start = anno_CDS_list[0].split('#')[0]
				second_start = anno_CDS_list[1].split('#')[0]
				if int(first_start) < int(second_start):
					CDS_start = anno_CDS_list[-1].split('#')[1]
					CDS_end = anno_CDS_list[0].split('#')[0]
					anno_CDS_list = anno_CDS_list[::-1]

		for each_CDS_region in anno_CDS_list:
			start = int(each_CDS_region.split('#')[0])
			end = int(each_CDS_region.split('#')[1])
			temp_seq = fetch_seq(genome, chrom, start-1, end, strand)
			seq = temp_seq[0:-1] + temp_seq[-1].lower()
			Full_CDS_seq += seq
		Full_CDS_seq = Full_CDS_seq[0:-1] + Full_CDS_seq[-1].upper()
		# after obtain full length of CDS sequence
		splicing_site_position = [i.start() for i in re.finditer("[a-z]", Full_CDS_seq)]
		peptides=Seq(Full_CDS_seq).translate()
		peptides_with_ss = lower_ss(peptides,splicing_site_position) #remove final stop codon-derived AA
		if peptides_with_ss[-1]=='*':   
			peptides_with_ss = peptides_with_ss[0:-1]

		final_pep_ID = '>'+str(chrom)+'_'+strand+'_'+str(CDS_start)+'_'+str(CDS_end)+'_CDS_'+ENST_ID+'_'+ENSG_ID+'_'+CDS_dict[ENST_ID]+'_withGTF'
		if transcript_type_dict[ENST_ID] == 'protein_coding':
			final_pep_ID = final_pep_ID + ':PC'
		elif transcript_type_dict[ENST_ID] == 'nonsense_mediated_decay':
			final_pep_ID = final_pep_ID + ':NMD'

		outf_protein.write(str(final_pep_ID)+'\n'+str(peptides_with_ss)+'\n')
		outf_transcript.write(str(final_pep_ID)+'\n'+str(final_seq)+'\n')

	else: # either ENST is not annotated or that paramter is not used
		# Obtain the full length sequence of transcripts
		final_seq = ''  #full mRNA sequence
		if re.findall(';',transcript_link):
			each_exon_list = transcript_link.split(';')
			if strand == '-':
				first_start = each_exon_list[0].split('#')[3]
				second_start = each_exon_list[1].split('#')[3]
				if int(first_start) < int(second_start):
					each_exon_list = each_exon_list[::-1]
			for each_exon in each_exon_list:
				exon_detail = each_exon.split('#')
				start = int(exon_detail[3])
				end = int(exon_detail[4])
				if start -1 >= end:
					mark_weird_isoform = 1
					break
				temp_seq = fetch_seq(genome, chrom, start-1, end, strand)
				seq = temp_seq[0:-1] + temp_seq[-1].lower()
				if final_seq=='':
					final_seq = seq
				else:
					final_seq = final_seq + seq
			if mark_weird_isoform == 1:
				continue
			#temp_final_seq = final_seq[0:-1] + final_seq[-1].upper()
			#final_seq = temp_final_seq
			if strand == '+':
				tail_end = int(each_exon_list[-1].split('#')[4])+1
				tail_3_nt = fetch_seq(genome, chrom, tail_end-1, tail_end+2, strand)
			elif strand == '-':
				tail_end = int(each_exon_list[-1].split('#')[3])
				tail_3_nt = fetch_seq(genome, chrom, tail_end-4, tail_end-1, strand)
			final_seq = final_seq + tail_3_nt

		else: # only one exon
			if int(this_exon[3])-1 >= int(this_exon[4]):
				mark_weird_isoform = 1
				continue
			this_exon_seq = fetch_seq(genome, chrom, int(this_exon[3])-1, int(this_exon[4]), strand)
			if strand == '+':
				tail_end = int(this_exon[4]) + 1
				tail_3_nt = fetch_seq(genome, chrom, tail_end-1, tail_end+2, strand)
			elif strand == '-':
				tail_end = int(this_exon[3])
				tail_3_nt = fetch_seq(genome, chrom, tail_end-4, tail_end-1, strand)
			final_seq = this_exon_seq + tail_3_nt

		if len(final_seq) > 100000: continue
		# Searching for right ORF
		start_codon_position = [i.start() for i in re.finditer('ATG', final_seq,flags=re.IGNORECASE)]
		stop_codon_position = [i.start() for i in re.finditer("TAG|TGA|TAA", final_seq,flags=re.IGNORECASE)]
		pair_codon = codon_judge(start_codon_position,stop_codon_position)
		if len(pair_codon) > 0:
			final_pep_ID = ''
			final_peptide_seq = ''
			peptide_length = 0
			final_pep_link = ''
			for each_pair in pair_codon:
				ss = each_pair.split('#')[0]
				ee = each_pair.split('#')[1]
				final_transcript = final_seq[int(ss):int(ee)+3]
				splicing_site_position = [i.start() for i in re.finditer("[a-z]", final_transcript)]
				absolute_position = genome_position(ss,ee,transcript_link)
				#print absolute_position
				peptides=Seq(final_transcript).translate()
				peptides_with_ss = lower_ss(peptides,splicing_site_position)
				#remove final stop codon-derived AA
				if peptides_with_ss[-1]=='*':   
					peptides_with_ss = peptides_with_ss[0:-1]
			
				this_pep_len = len(peptides_with_ss)
				this_pep_link = absolute_position[2]
				CDS_str = 'CDS:'+str(int(ss)+1)+'-'+str(int(ee)+3)
				this_pep_ID = '>'+str(chrom)+'_'+strand+'_'+str(absolute_position[0])+'_'+str(absolute_position[1])+'_'+CDS_str+'_'+ENST_ID+'_'+ENSG_ID+'_'+this_pep_link
				
				#(either use_annotion == 'No' or ENST_ID not in CDS_dict)
				if this_pep_len >= peptide_length:
					if this_pep_len == peptide_length:
						if len(this_pep_link.split(';')) <= len(final_pep_link.split(';')):
							continue
					final_pep_ID = this_pep_ID+'_withoutGTF'
					final_peptide_seq = peptides_with_ss
					final_pep_link = this_pep_link
					peptide_length = this_pep_len
					continue
			# Decide NMD
			NMD_result = check_NMD(final_pep_link,transcript_link,strand)
			if NMD_result:
				final_pep_ID = final_pep_ID + ':NMD'
			else:
				final_pep_ID = final_pep_ID + ':PC'
			# After loop, select the longest ORF
			outf_protein.write(str(final_pep_ID)+'\n'+str(final_peptide_seq)+'\n')
		else:  # No paired start/stop codon is found
			final_pep_ID = '>'+str(chrom)+'_'+strand+'_0_1_CDS:none_'+ENST_ID+'_'+ENSG_ID+'_0#1_withoutGTF'

		final_seq_original = final_seq[0:-3]
		outf_transcript.write(str(final_pep_ID)+'\n'+str(final_seq_original)+'\n')	

inf.close()
outf_transcript.close()
outf_protein.close()	
