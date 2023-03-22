import os,re,sys
from collections import defaultdict

inf_name = sys.argv[1]
outf_name = sys.argv[2]
outf = open(outf_name,'w')
inf = open(inf_name,'r')
key_list = []
last_transcript_annotation = ''
#print 'length of file',len(list(inf))

file_len_str = os.popen('wc -l %s' % inf_name).readlines()[0]
file_len = int(file_len_str.split(' ')[0])
#print 'file_len',file_len


for i,line in enumerate(inf):
	if line.startswith('#'):
		continue
	arr = line.strip().split('\t')
	if arr[2] == 'transcript':
		if len(key_list)>0:
			transcript_link = ';'.join(key_list)
			outf.write(transcript_link+'\t'+last_transcript_annotation+'\n')
		key_list = []
		continue
	elif arr[2] == 'exon':
		if re.findall('reference_id', arr[8]):
			raw_ENST_ID = re.findall('reference_id \"(.+?)\";',arr[8])[0]
		else:
			raw_ENST_ID = re.findall('transcript_id \"(.+?)\";',arr[8])[0]
		if re.findall('STRG', raw_ENST_ID):
			ENST_ID = raw_ENST_ID
		elif re.findall('_PAR_Y', raw_ENST_ID):
			ENST_ID = raw_ENST_ID.split('.')[0]+'-PAR-Y'
		else:
			ENST_ID = raw_ENST_ID.split('.')[0]
		if arr[0].startswith('chr'):
			pass
		else:
			if arr[0] in (list(map(str,range(1,23)))+['X','Y']):
				arr[0] = 'chr'+arr[0]
			else:
				continue
		key = ENST_ID+'#'+arr[0]+'#'+arr[6]+'#'+str(int(arr[3]))+'#'+str(arr[4])
		key_list.append(key)
		last_transcript_annotation = arr[1]
		if i == file_len-1: # the end of entire file
			transcript_link = ';'.join(key_list)
			outf.write(transcript_link+'\t'+last_transcript_annotation+'\n')
inf.close()
outf.close()
