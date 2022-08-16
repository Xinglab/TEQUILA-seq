import re,os,sys
from collections import defaultdict
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.path import Path as mPath
import matplotlib.patches as mpatches
import numpy as np


##### pre-set parameters ##############
fig = plt.figure(figsize=(12, 4),dpi=300)
ax = fig.add_subplot(111)

exon_height = 6
line_space = 13
font_size = 12
# parameter: plot size: 1000 * 600
####### input parameters ##################
CPM_dict = defaultdict()
exon_dict = defaultdict()
corrected_exon_dict = defaultdict()
corrected_intron_dict = defaultdict()

min_pos = 10000000000
max_pos = 0
num_trans = 0
trans_ID_list = []
top_rank_isoform = []
bed_list_name = sys.argv[1]
gene_name = sys.argv[2]
gene_column = int(sys.argv[3])

tool = 'ESPRESSO'
strand = '+'
if len(sys.argv)>4:
	tool = sys.argv[4]
#if len(sys.argv)>5:
#	strand = sys.argv[5]
max_total_stru_len = 0
text_x = ''

def decide_strand(input_strand):
	if input_strand == '+':
		left_end = "5'"
		right_end = "3'"
	elif input_strand == '-':
		left_end = "3'"
		right_end = "5'"
	return [left_end, right_end]

strand_dict = defaultdict()
tag_dict = defaultdict()
abundance_index_dict = defaultdict()
trans_ID_abundance_list = []
with open(bed_list_name,'r') as list_inf:
	for list_index,list_line in enumerate(list_inf):
		num_trans += 1
		first_arr = list_line.strip().split('\t')
		file_name = first_arr[1].split('/')[-1]
		this_strand = first_arr[2]
		file_name_list = file_name.split('.')[0].split('_')
		group = '_'.join(file_name_list[0:gene_column])
		#gene_name = file_name_list[2]
		if tool == 'ESPRESSO':
			trans_ID = '_'.join(file_name_list[gene_column+1:len(file_name_list)])
		elif tool == 'FLAMES':
			trans_ID = re.findall('Target.+ENS.+(ENS.+)',file_name.split('.')[0])[0]
		trans_ID_list.append(trans_ID)
		strand_dict[trans_ID] = this_strand
		abundance_index_dict[trans_ID] = int(first_arr[0])
		trans_ID_abundance_list.append(int(first_arr[0]))
		if len(first_arr)>3:
			tag_dict[trans_ID] = first_arr[3]
		with open(first_arr[1],'r') as inf:
			for index,line in enumerate(inf):
				arr = line.strip().split('\t')
				if index == 0: continue
				CPM_dict[trans_ID] = float(arr[3])
				if trans_ID not in exon_dict:
					exon_dict[trans_ID] = []
				exon_dict[trans_ID].append('#'.join(arr[0:3]))
				min_pos = min(min_pos, float(arr[1]))
				max_pos = max(max_pos, float(arr[2]))
corrected_max_pos = max_pos - min_pos


min_i2e_len_ratio = 10000000000
for each_trans in trans_ID_list:
	if len(exon_dict[each_trans]) >= 2:
		for i in range(len(exon_dict[each_trans])-1):
			exon_len = float(exon_dict[each_trans][i].split('#')[2]) - float(exon_dict[each_trans][i].split('#')[1])
			intron_len = (float(exon_dict[each_trans][i+1].split('#')[1]) - float(exon_dict[each_trans][i].split('#')[2]))
			next_exon_len = float(exon_dict[each_trans][i+1].split('#')[2]) - float(exon_dict[each_trans][i+1].split('#')[1])
			i2e_len_ratio = min(intron_len/exon_len, intron_len/next_exon_len)
			#print exon_len, intron_len, next_exon_len, i2e_len_ratio
			min_i2e_len_ratio = min(min_i2e_len_ratio, i2e_len_ratio)

### sort color list ###
#index2color_list_original = ['#7D4655','#B8604D','#E0BBA4','#86a65d','#438B80','#413E60']
#index2color_list = []
#for each_abun_index in trans_ID_abundance_list:
#	index2color_list.append(index2color_list_original[int(each_abun_index)-1])

#index2color_list_original = ['#a90f0f','#458baf','#cde1eb','#87b6cf','#295369','#94a3aa']
#index2color_list_original = ['#a90f0f','#458baf','#cde1eb','#87b6cf','#295369']
#index2color_list_original = ['#FC402D','#08589e','#2b8cbe','#4eb3d3','#7bccc4']
#index2color_list_original = ['#FC402D','#2171b5','#4292C6','#6baed6','#9ecae1']
index2color_list_original = ['#F31E08','#2171b5','#4292C6','#6baed6','#9ecae1']
index2color_list = index2color_list_original[::-1]


with open(bed_list_name,'r') as list_inf_2:
	for list_index,list_line in enumerate(list_inf_2):
		file_name = list_line.strip().split('\t')[1].split('/')[-1]
		file_name_list = file_name.split('.')[0].split('_')
		group = '_'.join(file_name_list[0:gene_column])
		#gene_name = file_name_list[2]
		if tool == 'ESPRESSO':
			trans_ID = '_'.join(file_name_list[gene_column+1:len(file_name_list)])
		elif tool == 'FLAMES':
			trans_ID = re.findall('Target.+ENS.+(ENS.+)',file_name.split('.')[0])[0]
		### adjust initial current_pos
		initial_pos = (float(exon_dict[trans_ID][0].split('#')[1]) - min_pos)
		current_pos = initial_pos
		if len(exon_dict[trans_ID]) >= 2:
			for i in range(len(exon_dict[trans_ID])-1):
				exon_len = (float(exon_dict[trans_ID][i].split('#')[2]) - float(exon_dict[trans_ID][i].split('#')[1]))
				intron_len = (float(exon_dict[trans_ID][i+1].split('#')[1]) - float(exon_dict[trans_ID][i].split('#')[2]))
				if trans_ID not in corrected_exon_dict:
					corrected_exon_dict[trans_ID] = []
				if trans_ID not in corrected_intron_dict:
					corrected_intron_dict[trans_ID] = []
				current_exon_pos = str(current_pos)+':'+str(current_pos+exon_len)
				current_pos = current_pos + exon_len
				current_intron_pos = str(current_pos)+':'+str(current_pos+intron_len)
				current_pos = current_pos + intron_len
				corrected_exon_dict[trans_ID].append(current_exon_pos)
				corrected_intron_dict[trans_ID].append(current_intron_pos)
			last_exon_len = float(exon_dict[trans_ID][-1].split('#')[2]) - float(exon_dict[trans_ID][-1].split('#')[1])
			current_exon_pos = str(current_pos)+':'+str(current_pos+last_exon_len)
			corrected_exon_dict[trans_ID].append(current_exon_pos)
			current_pos = current_pos+last_exon_len

		#### single exon #####
		else:
			exon_len = (float(exon_dict[trans_ID][0].split('#')[2]) - float(exon_dict[trans_ID][0].split('#')[1]))
			if trans_ID not in corrected_exon_dict:
				corrected_exon_dict[trans_ID] = []
			current_exon_pos = str(current_pos)+':'+str(current_pos+exon_len)
			current_pos = current_pos + exon_len
			corrected_exon_dict[trans_ID].append(current_exon_pos)
                
		####
		final_pos = current_pos
		max_total_stru_len = max(max_total_stru_len, current_pos)
		if text_x == '':
			text_x = -0.6*max_total_stru_len
		#print trans_ID, corrected_exon_dict[trans_ID], corrected_intron_dict[trans_ID], current_pos

		######## plot ######
		selected_color = index2color_list[list_index]
		#if (trans_ID.startswith('ESPRESSO') or (re.findall('_',trans_ID))):
		#	selected_color = '#E04036'
		#else:
		#	selected_color = '#2893D4'
		######## draw intron line ########	
		if len(corrected_exon_dict[trans_ID]) >= 2:
			intron_left = float(corrected_exon_dict[trans_ID][0].split(':')[1])
			intron_right = float(corrected_exon_dict[trans_ID][-1].split(':')[0])
			intron_x_region = np.linspace(intron_left, intron_right, round(max(10,0.1*(intron_right-intron_left))))
			intron_y_region = np.array([line_space*(list_index+1)]*len(intron_x_region))
			#print intron_x_region
			if strand_dict[trans_ID] == '-':
				plt.plot(-intron_x_region,intron_y_region,c=selected_color,lw=1)
			else:
				plt.plot(intron_x_region,intron_y_region,c=selected_color,lw=1)

		######## draw exon box ########
		for each_exon in corrected_exon_dict[trans_ID]:
			exon_left = float(each_exon.split(':')[0])
			exon_right = float(each_exon.split(':')[1])
			exon_length = exon_right - exon_left
			#left_top_point_y = line_space*(list_index+1) + exon_height/2
			#exon_box = mpatches.Rectangle((exon_left,left_top_point_y),exon_length,exon_height,linewidth=2,edgecolor='r',facecolor='none')
			#ax.add_patch(exon_box)
			exon_x_region = np.linspace(exon_left, exon_right, round(max(10,0.1*exon_length)))
			exon_y_region_top = np.array([line_space*(list_index+1)+exon_height/2]*len(exon_x_region))
			exon_y_region_bottom = np.array([line_space*(list_index+1)-exon_height/2]*len(exon_x_region))
			#print exon_left,exon_right,exon_length,exon_x_region
			if strand_dict[trans_ID] == '-':
				ax.fill_between(-exon_x_region, exon_y_region_top, exon_y_region_bottom, facecolor=selected_color, alpha=1)
			else:
				ax.fill_between(exon_x_region, exon_y_region_top, exon_y_region_bottom, facecolor=selected_color, alpha=1)

		displayed_trans_ID = re.sub('#Target_RBP_SDA|#Target_IDT|#Target_Ctrl|#Direct_RNA|#Target_SDA|_Sendai_03451_hICset4_D7','',trans_ID)
		if re.findall('ENST', displayed_trans_ID):
			displayed_trans_ID = displayed_trans_ID +' (%s)'%tag_dict[trans_ID]
		############# draw 5' and 3' end ###########
		[left_end, right_end] = decide_strand(strand_dict[trans_ID])
		if strand_dict[trans_ID] == '-':
			plt.text(-1*corrected_max_pos, line_space*(list_index+0.9)+exon_height, displayed_trans_ID, fontsize=font_size+1)
			plt.text(-1.03*corrected_max_pos, line_space*(list_index+1)-2, '5\'', fontsize=font_size+2)
		else:
			plt.text(0, line_space*(list_index+0.9)+exon_height, displayed_trans_ID, fontsize=font_size+1)
			plt.text(-0.03*corrected_max_pos, line_space*(list_index+1)-2, '5\'', fontsize=font_size+2)
		#plt.text(1.1*corrected_max_pos, line_space*(list_index+1)-3, right_end, fontsize=11)



#print ('corrected_max_pos', corrected_max_pos)
#print ('max_total_stru_len', max_total_stru_len)

if strand_dict[trans_ID] == '-':
	ax.set_xlim(-1.1*max_total_stru_len, 0.1*max_total_stru_len)
	plt.text(-0.55*max_total_stru_len, line_space*(num_trans+1), gene_name, fontsize=font_size+2)
else:
	plt.text(0.45*max_total_stru_len, line_space*(num_trans+1), gene_name, fontsize=font_size+2)
	ax.set_xlim(-0.1*max_total_stru_len, 1.1*max_total_stru_len)
ax.set_ylim(0, line_space*(num_trans+1.5))
#plt.axis('equal')
plt.axis('off')
plt.subplots_adjust(top=1, bottom=0, right=1, left=0, hspace=0, wspace=0)
plt.margins(0,0)
plt.show()
outf_name = re.sub('.txt','_single_color.pdf',bed_list_name)
plt.savefig(outf_name,dpi=300, pad_inches=0)

outf_name_2 = re.sub('.txt','_single_color.png',bed_list_name)
plt.savefig(outf_name_2, dpi=300, pad_inches=0)
















##### isoform ######
'''
Gene_name = 'SCAMP5'
Transcript_ID = 'ENST00000565989'
Uni_ID = 'H3BTD1'
topology = '0i58-78o86-86i106-116' #extracellular change
topology = '0i374-400o414m447o461-482i500' #extracellular change
ss_list = [1, 20, 63]
exon_list = ['Exon1','Exon2','Exon3','Exon4']


if re.findall('-',topology)==False:             # It must have transmembrane domain
	print 'Wrong: no TM'
protein_length = map(int,re.split('o|-|i|m',topology))[-1]

max_extra_len = 0
max_intra_len = 0
for each_domain in re.split('-|m',topology):
	if re.findall('o',each_domain):
		extra_length = float(each_domain.split('o')[1]) - float(each_domain.split('o')[0]) 
		if extra_length > max_extra_len:
			max_extra_len = extra_length
	elif re.findall('i',each_domain):
		intra_length = float(each_domain.split('i')[1]) - float(each_domain.split('i')[0]) 
		if intra_length > max_intra_len:
			max_intra_len = intra_length
print 'maximum extracellular domain length:', max_extra_len
print 'maximum intracellular domain length:', max_intra_len

TM_region_list = []
for each_region_temp in re.split('o|i',topology):
	if re.findall('-',each_region_temp):
		TM_region_list.append(each_region_temp)

expected_figure_width = 800 #height is 20, set length is 40
shrink_factor = float(expected_figure_width) / float(protein_length)
shrink_factor_trans = float(expected_figure_width) / (float(protein_length) * 2)

########## protein curve  ##############
print "calculating protein curves"
coor_list = map(int,re.split('o|-|i|m',topology))   #[0,10,]
split_list = re.findall('o|-|i|m',topology)  #['o','-',]
current_height = 0
current_width = 0
start_point = 0
ctl_point = 0
end_point = 0
TM_distance = 0
curve_x_list = []
curve_y_list = []
intra_membrane_list = []



########## give each exon different color ######
# exon_region_list: [0, 30, 60, 90, 120, 150]
# TM_region_list: [20-40,95-115]
exon_region_list = [0] + ss_list + [protein_length]
curve_x_list = np.array(curve_x_list)
curve_y_list = np.array(curve_y_list)
color_sequential = ['#00429d', '#2e59a8', '#4771b2', '#5d8abd', '#73a2c6', '#8abccf', '#a5d5d8', '#c5eddf', '#ffffe0']
color_diverging = ['#00429d', '#4771b2', '#73a2c6', '#a5d5d8', '#ffbcaf', '#f4777f', '#cf3759', '#93003a'] * 10 
color_double = ['#99b898','#feceab','#99b898','#feceab','#99b898','#feceab','#99b898','#feceab','#99b898','#feceab',]
#curve_plot_list = []
########## draw protein curve ############
print "plotting protein curves"
for j in range(len(exon_region_list)-1):
	[adjust_exon_region_left, adjust_exon_region_left_TM_percent] = convert_x_coordinate(exon_region_list[j],topology,shrink_factor)
	[adjust_exon_region_right, adjust_exon_region_right_TM_percent] = convert_x_coordinate(exon_region_list[j+1],topology,shrink_factor)
	this_exon_x = curve_x_list[np.where((curve_x_list>=adjust_exon_region_left) & (curve_x_list<=adjust_exon_region_right))]
	this_exon_y = curve_y_list[np.where((curve_x_list>=adjust_exon_region_left) & (curve_x_list<=adjust_exon_region_right))]
	# part of exon is in left TM domain:
	if adjust_exon_region_left_TM_percent > 0:
		point_num_in_TM = len(np.where(np.around(this_exon_x,decimals=2) == round(adjust_exon_region_left,2))[0])
		exon_proportion_in_TM = float(point_num_in_TM * (1-adjust_exon_region_left_TM_percent))
		delete_exon_x_in_TM_index = np.where(np.around(this_exon_x,decimals=2) == round(adjust_exon_region_left,2))[0][:-exon_proportion_in_TM]
		this_exon_x = np.delete(this_exon_x,delete_exon_x_in_TM_index)
		this_exon_y = np.delete(this_exon_y,delete_exon_x_in_TM_index)
	# part of exon is in right TM domain:
	if adjust_exon_region_right_TM_percent > 0:
		point_num_in_TM = len(np.where(np.around(this_exon_x,decimals=2) == round(adjust_exon_region_right,2))[0])
		exon_proportion_in_TM = float(point_num_in_TM * adjust_exon_region_right_TM_percent)
		delete_exon_x_in_TM_index = np.where(np.around(this_exon_x,decimals=2) == round(adjust_exon_region_right,2))[0][exon_proportion_in_TM:]
		this_exon_x = np.delete(this_exon_x,delete_exon_x_in_TM_index)
		this_exon_y = np.delete(this_exon_y,delete_exon_x_in_TM_index)
	#print this_exon_x,this_exon_y
	#selected_color = color_double[j]
	selected_color = color_diverging[exon_color_map_dict[exon_list[j]]]
	plt.plot(this_exon_x,this_exon_y,c=selected_color,lw=3)

	this_exon_y_TM = this_exon_y[np.where((this_exon_y>membrane_bottom) & (this_exon_y<membrane_top))]
	this_exon_x_TM = this_exon_x[np.where((this_exon_y>membrane_bottom) & (this_exon_y<membrane_top))]
	this_exon_x_TM_non_redundancy = list(set(this_exon_x_TM))
	#print this_exon_x_TM_non_redundancy
	for each_x_value in this_exon_x_TM_non_redundancy:
		this_exon_x_TM_splitted = this_exon_x_TM[np.where(this_exon_x_TM == each_x_value)]
		this_exon_y_TM_splitted = this_exon_y_TM[np.where(this_exon_x_TM == each_x_value)]
		if whether_x_in_intramembrane(each_x_value,intra_membrane_list):
			continue
		else:
			plt.plot(this_exon_x_TM_splitted,this_exon_y_TM_splitted,c=selected_color,lw=8)
			pass

########## draw transcript box ############
print "plotting transcript box"
horizontal_c = 50
for j in range(len(exon_region_list)-1):
	width = exon_region_list[j+1] - exon_region_list[j] 
	trans_x_region = np.linspace(exon_region_list[j]*shrink_factor_trans*1.2 + horizontal_c, exon_region_list[j+1]*shrink_factor_trans*1.2 + horizontal_c, 1*width)
	trans_y_region_top = np.array([trans_box_height+10]*len(trans_x_region))
	trans_y_region_bottom = np.array([trans_box_height-10]*len(trans_x_region))
	#selected_color = color_double[j]
        selected_color = color_diverging[exon_color_map_dict[exon_list[j]]]
	ax.fill_between(trans_x_region, trans_y_region_top, trans_y_region_bottom, facecolor=selected_color, alpha=0.1)	
	if len(exon_region_list)<=25:
		plt.text(np.mean(trans_x_region)-5, trans_box_height-20, str(exon_list[j]), fontsize=4)

########## draw transcript structure ############
print "plotting transcript structure"
transcript_inf = open(sample_processed_gtf_inf_name,'r')
transcript_structure_list = []
exon_structure_tag_list = []
transcript_strand = ''
for line in transcript_inf:
	line = line.strip()
	arr = line.split('\t')
	ID = arr[0].split('#')[0].split('.')[0]
	if ID == Transcript_ID:
		transcript_strand = arr[0].split(';')[0].split('#')[2]
		for each_segment in arr[0].split(';'):
			exon_pair = each_segment.split('#')[3]+':'+each_segment.split('#')[4]
			transcript_structure_list.append(exon_pair)
		break
transcript_inf.close()


transcript_total_length = float(transcript_structure_list[-1].split(':')[1]) - float(transcript_structure_list[0].split(':')[0])
transcript_shrink_factor = float(right_most)/transcript_total_length

exon_structure_tag_list = exon_trans_list
if transcript_strand == '-':
	transcript_structure_list = transcript_structure_list[::-1]


transcript_structure_height = -180
current_left = 0

for k in range(len(transcript_structure_list)-1):
	exon_width = float(transcript_structure_list[k].split(':')[1]) - float(transcript_structure_list[k].split(':')[0])
	exon_x_region = np.linspace(current_left*transcript_shrink_factor + horizontal_c, (current_left+exon_width)*transcript_shrink_factor + horizontal_c, 0.1*exon_width)
	exon_y_region_top = np.array([transcript_structure_height+10]*len(exon_x_region))
	exon_y_region_bottom = np.array([transcript_structure_height-10]*len(exon_x_region))
	#selected_color = color_double[j]
	selected_color = '#AFABAB'
	#plt.plot(exon_x_region,exon_y_region,c=selected_color,lw=4)
	ax.fill_between(exon_x_region, exon_y_region_top, exon_y_region_bottom, facecolor=selected_color, alpha=0.1)
	if len(exon_structure_tag_list) == len(transcript_structure_list):
		plt.text(np.mean(exon_x_region)-3, transcript_structure_height-20, str(exon_structure_tag_list[k]), fontsize=4)
	current_left = current_left + exon_width

	if transcript_strand == '+':
		intron_width = float(transcript_structure_list[k+1].split(':')[0]) - float(transcript_structure_list[k].split(':')[1])
	elif transcript_strand == '-':
		intron_width = abs(float(transcript_structure_list[k].split(':')[0]) - float(transcript_structure_list[k+1].split(':')[1]))
	intron_x_region = np.linspace(current_left*transcript_shrink_factor + horizontal_c, (current_left+intron_width)*transcript_shrink_factor + horizontal_c, 0.1*intron_width)
	intron_y_region = np.array([transcript_structure_height]*len(intron_x_region))
	selected_color = 'black'
	plt.plot(intron_x_region,intron_y_region,c=selected_color,lw=0.8)
	#print k,intron_x_region
	current_left = current_left + intron_width
	#if len(exon_region_list)<=25:
	#	plt.text((exon_region_list[j]+exon_region_list[j+1])*0.6*shrink_factor_trans + 3, trans_box_height-2, str(exon_list[j]), fontsize=4)

exon_width = float(transcript_structure_list[-1].split(':')[1]) - float(transcript_structure_list[-1].split(':')[0])
exon_x_region = np.linspace(current_left*transcript_shrink_factor + horizontal_c, (current_left+exon_width)*transcript_shrink_factor + horizontal_c, 0.1*exon_width)
exon_y_region_top = np.array([transcript_structure_height+10]*len(exon_x_region))
exon_y_region_bottom = np.array([transcript_structure_height-10]*len(exon_x_region))
#selected_color = color_double[j]
selected_color = '#AFABAB'
ax.fill_between(exon_x_region, exon_y_region_top, exon_y_region_bottom, facecolor=selected_color, alpha=0.1)
if len(exon_structure_tag_list) == len(transcript_structure_list):
	plt.text(np.mean(exon_x_region)-3, transcript_structure_height-20, str(exon_structure_tag_list[-1]), fontsize=4)
#plt.plot(exon_x_region,exon_y_region,c=selected_color,lw=4)


########## membrane box #################
x = np.arange(left_most-120, right_most+20, 1)
mem_cur_1_verts = np.array([membrane_top] * len(x))
mem_cur_2_verts = np.array([membrane_bottom] * len(x))

ax.fill_between(x, mem_cur_1_verts, mem_cur_2_verts, facecolor='#DEEBF7', alpha=0.1)
plt.text(-100,trans_box_height+50, Gene_name+':'+Transcript_ID, fontsize=12)
plt.text(-110,trans_box_height, 'ORF in transcript', fontsize=8)
plt.text(-130,transcript_structure_height, 'Transcript structure', fontsize=8)
plt.text(-110,transcript_structure_height-20, 'Strand:'+transcript_strand, fontsize=8)
plt.text(-100,0,'Cell membrane',fontsize=8)
plt.text(-100,extra_height-30,'Extracellular',fontsize=8)
plt.text(-100,intra_height+30,'Cytoplasmic',fontsize=8)

############ coordinates of protein ############


########## figure parameters ########
ax.set_xlim(-120,right_most+20)
ax.set_ylim(-120,120)
plt.axis('equal')
plt.axis('off')

#plt.plot(curve_x_list,curve_y_list)
#plt.show()

#outf_name = '/Users/violin/Desktop/'+Gene_name+'_'+Transcript_ID+'_TMHMM.png'
outf_name = './'+Gene_name+'_'+Transcript_ID+'_TMHMM.eps'
plt.savefig(outf_name,dpi=300)

'''
