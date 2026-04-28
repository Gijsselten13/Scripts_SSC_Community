## Imports
import glob
import numpy as np
import pandas as pd
import seaborn as sns
import scipy.stats as stats
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans 
from sklearn.decomposition import PCA
from sklearn.metrics import accuracy_score
from matplotlib.colors import LinearSegmentedColormap, BoundaryNorm, TwoSlopeNorm

path <- ""

## Functions (v1)

# Divide and organize data
def divide_and_organize_data(kos, ko_list, isos, taxonomy):
	# Define other side of isolates
	nontop_isos = [i for i in kos.columns.to_list() if i not in isos]
	# Split data
	top = kos.loc[ko_list,isos]
	nontop = kos.loc[ko_list,nontop_isos]
	# Binarize and merge taxonomy
	top = pd.merge((top.T > 0).astype(int), taxonomy[['isolate','family']],
               left_index=True, right_on='isolate').set_index('isolate')
	nontop = pd.merge((nontop.T > 0).astype(int), taxonomy[['isolate','family']],
                  left_index=True, right_on='isolate').set_index('isolate')	# Count families
	fam_top = top[['family']].reset_index().groupby('family').count()
	fam_nontop = nontop[['family']].reset_index().groupby('family').count()
	# KO count
	ko_top = top.groupby('family').sum()
	ko_nontop = nontop.groupby('family').sum()
	# Define order rows
	orderrow_top = fam_top.sort_values('isolate', ascending=False).index.to_list()
	orderrow_nontop = fam_nontop.sort_values('isolate', ascending=False).index.to_list()
	# Define order columns
	ordercol = top.drop('family', axis=1).sum(axis=0).sort_values().index[::-1].to_list()
	# Define rownames
	rownames_top = [f"{i} (n={j})" for i,j in fam_top.sort_values('isolate', ascending=False).to_dict()['isolate'].items()]
	rownames_nontop = [f"{i} (n={j})" for i,j in fam_nontop.sort_values('isolate', ascending=False).to_dict()['isolate'].items()]
	# Order tables
	ko_top = ko_top.loc[orderrow_top,ordercol]
	ko_nontop = ko_nontop.loc[orderrow_nontop,ordercol]
	# Compute distribution
	distribution_top = ko_top.sum(axis=0)
	distribution_nontop = ko_nontop.sum(axis=0)
	# Prevalence
	prevalence_top = ko_top.divide(fam_top['isolate'], axis=0).multiply(100).loc[orderrow_top,ordercol]
	prevalence_nontop = ko_nontop.divide(fam_nontop['isolate'], axis=0).multiply(100).loc[orderrow_nontop,ordercol]
	# Organize data
	res = {'top': {'dist': distribution_top, 'prev': prevalence_top, 'rownames': rownames_top, 'orig': top, 'fam_nisos': fam_top}, 'nontop': {'dist': distribution_nontop, 'prev': prevalence_nontop, 'rownames': rownames_nontop, 'orig': nontop, 'fam_nisos': fam_nontop}}
	# Return
	return(res)

# Plot data
def plot_upset(dist, prev, rownames, cmap, norm, figsize, outfile):
	# Melt data and remove zeros
	prev = prev.melt(ignore_index=False).reset_index()
	prev['value'] = prev['value'].replace({0:np.nan})
	
	# Fomat
	fig,ax = plt.subplots(ncols=1, nrows=2, figsize=figsize)

	### Distribution
	sns.barplot(dist, color="grey", edgecolor="white", linewidth=0.5, ax=ax[0])
	ax[0].set_xticklabels([])
	ax[0].tick_params(axis='y', labelsize=22)
	
	### Prevalence
	# Background
	ax[1].set_facecolor("black")
	# Scatterplot/Heatmap
	sns.scatterplot(prev, x="variable", y = "family", size="value", hue="value", palette=cmap, edgecolor=None, ax=ax[1], legend=False, zorder=1)
	plt.margins(x=0.005, y=0.025)
	plt.xlabel("")
	plt.ylabel("")
	plt.yticks(size=22)
	plt.xticks(size=12, rotation=90)
	# Add vertical lines
	for xtick in ax[1].get_xticks():
		ax[1].axvline(x=xtick, color='lightgray', linestyle='-', alpha=0.5, zorder=0)
	# Rownames
	ax[1].set_yticklabels(rownames)
	# Color scale
	sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
	sm.set_array([])
	cax = fig.add_axes([ax[1].get_position().x1+0.01, ax[1].get_position().y0+0.05, 0.01, ax[1].get_position().height/1.5])
	cbar = ax[1].figure.colorbar(sm, cax=cax)
	cbar.ax.tick_params(labelsize=26)
	## Save and show plot
	plt.savefig(outfile, format="pdf", bbox_inches="tight", dpi=300)
	plt.tight_layout()
	plt.show()


## Load data

file_name_SSC = "SSC_taxonomy_GTDB.tsv"
file_name_isos_90_266 = "Barnard_stat/isos_90perc_266kos.tsv"
file_name_kos = "KO_genome/KO_SSC.tsv"
file_name_266 = "Barnard_stat/266_KOs.txt"

full_path_SSC = os.path.join(path, file_name_SSC)
full_path_isos = os.path.join(path, file_name_isos_90_266)
full_path_KOs = os.path.join(path, file_name_kos)
full_path_266 = os.path.join(path, file_name_266)

# GTDB taxonomy
taxonomy = pd.read_table(full_path_SSC)

# KO Data
kos = pd.read_table(full_path_KOs).set_index('sequence')
kos266 = pd.read_table(full_path_266, header=None)[0].to_list()

# Isolates lists: 266
isos_90_266 = pd.read_table(file_name_isos_90_266, header=None)[0].to_list()

# Compute distribution and prevalence: 266 KOs
res_90_266 = divide_and_organize_data(kos, kos266, isos_90_266, taxonomy)

## Color map

# Define section colors
boundaries = [0, 30, 60, 90, 500]
colors = ["lightgrey", "salmon", "gold", "#66bd63"]
cmap = LinearSegmentedColormap.from_list("custom", colors, N=len(colors))
norm = BoundaryNorm(boundaries, cmap.N, clip=True)
cmap

## Visualize
# Top 90 266

# Save results
file_name_out = "results/Figure_S30b_upset.top.90_266.pdf"
file_name_out_full = os.path.join(path, file_name_out)

plot_upset(dist=res_90_266['top']['dist'], prev=res_90_266['top']['prev'], rownames=res_90_266['top']['rownames'], cmap=cmap, norm=norm, figsize=(50,16), outfile=file_name_out_full)
