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
#index2color_list_original = ['#F31E08','#2171b5','#4292C6','#6baed6','#9ecae1']
index2color_list_original = ['#FF5F42','#003E7F','#0068AF','#5495E1','#A1C2E8']
index2color_list_original = index2color_list_original[0:num_trans]
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
				plt.plot(-intron_x_region,intron_y_region,c="black",lw=0.5)
			else:
				plt.plot(intron_x_region,intron_y_region,c="black",lw=0.5)

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
				ax.fill_between(-exon_x_region, exon_y_region_top, exon_y_region_bottom, facecolor=selected_color, alpha=1, edgecolor="black", linewidth=0.5, zorder=100)
			else:
				ax.fill_between(exon_x_region, exon_y_region_top, exon_y_region_bottom, facecolor=selected_color, alpha=1, edgecolor="black", linewidth=0.5, zorder=100)

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

