import sys
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt


def heatmap(input, ax=None,cbar_kw=None, cbarlabel="", **kwargs):
	#Parse data
	dp=pd.read_table(input,header=None)
	mat=dp.to_numpy()
	nrow,ncol=mat.shape
	chr=mat[0,0]
	group_vec=mat[0,6:ncol]
	loci_vec=mat[1:nrow,1]
	freq_mat=mat[1:nrow,6:ncol].T

	# Read dir name
	dir=""
	for i in range(0,len(input.split('/'))-1):
		dir+=input.split('/')[i]+"/"
		
	# Parse file name
	file_name=input.split('/')[-1]
	tar_group=file_name.split('_')[0]
	seg_start=file_name.split('_')[-2]
	seg_end=file_name.split('_')[-1].split('.')[0]
	chr=""
	for i in range(1,len(file_name.split('_'))-2):
		if i!=len(file_name.split('_'))-3:
			chr+=file_name.split('_')[i]+"_"
		else:
			chr+=file_name.split('_')[i]
	
	for i in range(0,ncol-6):
		for j in range(0,nrow-1):
			num=int(freq_mat[i][j].split('/')[0])
			denom=int(freq_mat[i][j].split('/')[1])
			if denom!=0:
				freq=int(freq_mat[i][j].split('/')[0])/int(freq_mat[i][j].split('/')[1])
			else:
				freq=np.nan
			freq_mat[i][j]=freq

	freq_mat_float=np.vstack((freq_mat[:, :])).astype(float)

	#Initialize plotting parameters
	if ax is None:
		ax = plt.gca()

	if cbar_kw is None:
		cbar_kw = {}	

	#Making heatmap
	im=ax.imshow(freq_mat_float,aspect='auto',**kwargs)

 	# Create colorbar
	cbar = ax.figure.colorbar(im,ax=ax,shrink=0.5,**cbar_kw)
	cbar.ax.set_ylabel(cbarlabel,size=5,rotation=-90, va="bottom")
	cbar.ax.tick_params(labelsize=4)

	# Show all ticks and label them with the respective list entries
	ax.set_xticks(np.arange(len(loci_vec)), labels=loci_vec)
	ax.set_yticks(np.arange(len(group_vec)), labels=group_vec)

	# Let the horizontal axes labeling appear on top.
	ax.tick_params(top=False, bottom=True,labeltop=False, labelbottom=True)

	# Rotate the tick labels and set their alignment.
	plt.setp(ax.get_xticklabels(), size=2.5,rotation=60, ha="right",rotation_mode="anchor")
	plt.setp(ax.get_yticklabels(), size=3)

	# Set plot title
	plt.xlabel("Variant genomic location",fontsize=6)
	plt.ylabel("Sample group",fontsize=6)
	plt.title("Enrichment pattern over "+tar_group+"-specific signature segment on chr "+chr,fontsize=8)

    # Turn spines off and create white grid.
	ax.spines[:].set_visible(False)
	ax.grid(which="minor", color="white", linestyle='-', linewidth=0.5)
	ax.tick_params(which="minor", bottom=False, left=False)
	ax.set_xticks(np.arange(len(loci_vec)+1)-0.5, minor=True)
	ax.set_yticks(np.arange(len(group_vec)+1)-0.5, minor=True)
	plt.rcParams["figure.autolayout"] = True

	op_fig_name=dir+tar_group+"_"+chr+"_"+seg_start+"_"+seg_end+".png"

	return im,cbar,op_fig_name

file=sys.argv[1]
fig,ax=plt.subplots()
im,cbar,op_fig=heatmap(file,ax=ax,cmap="YlGn",cbarlabel="Signature-carrying sample frequency")
fig.tight_layout()
plt.savefig(op_fig,format="png",dpi=800)


