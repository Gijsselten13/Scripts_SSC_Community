## Imports
import glob
import numpy as np
import pandas as pd
import scipy.stats as stats
import os

path = ""

## Functions

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
                  left_index=True, right_on='isolate').set_index('isolate')
	# Count families
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
	res = {'top': {'dist': distribution_top, 'prev': prevalence_top, 'rownames': rownames_top, 'orig': top, 'fam_nisos': fam_top, 'ko_top': ko_top}, 'nontop': {'dist': distribution_nontop, 'prev': prevalence_nontop, 'rownames': rownames_nontop, 'orig': nontop, 'fam_nisos': fam_nontop, 'ko_nontop': ko_nontop}}
	# Return
	return(res)

# Barnard Test
def barnard_test(top_present, nontop_present, top_isos, nontop_isos):
	# Obtain values for contingency table
	top_absent = top_isos - top_present
	nontop_absent = nontop_isos - nontop_present
	# Create table
	m = np.array([[top_present, nontop_present], [top_absent, nontop_absent]])
	# Compute Barnard test
	exact_test = stats.barnard_exact(m, alternative='two-sided', pooled=False)
	return [exact_test.statistic, exact_test.pvalue, top_present, top_absent, nontop_present, nontop_absent]

# Match top to non-top comparison
def filter_families(top, nontop, minval):
	# Define families to keep
	topfams = top['fam_nisos']
	keep = topfams[topfams['isolate'] >= minval].index.to_list()
	# Filter tables
	filt_top = top['ko_top'].loc[keep,:]
	filt_nontop = nontop['ko_nontop'].loc[keep,:]
	# Return
	return [filt_top, filt_nontop]

# Function to run Barnard's Test per KO for each family between Top and Non-Top
def test_group_proportions(top, nontop, minval):
	# Dictionary to store results
	res = []
	# Filter families given minimum number of isolates
	ko_top, ko_nontop = filter_families(top, nontop, minval)
	# Define families to test
	fams = ko_top.index.to_list()
	# Define KOs
	kos = ko_top.columns.to_list()
	# Loop KOs and families
	for ko in kos:
		for fam in fams:
			# Define KO values
			top_present = ko_top.loc[fam,ko]
			nontop_present = ko_nontop.loc[fam,ko]
			# Define number of isolates
			top_isos = top['fam_nisos'].loc[fam,:].to_list()[0]
			nontop_isos = nontop['fam_nisos'].loc[fam,:].to_list()[0]
			# Compute exact test
			statistic, pvalue, top_present, top_absent, nontop_present, nontop_absent = barnard_test(top_present, nontop_present, top_isos, nontop_isos)
			# Append results
			res.append([ko, fam, top_present, top_absent, nontop_present, nontop_absent, statistic, pvalue])
	# Results to dataframe
	res = pd.DataFrame(res)
	# Add columns
	res.columns = ['ko','family','TP','TA','NP','NA','statistic','pvalue']
	# Return
	return res


## Execute

#set files
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

# Isolates lists
isos_90_266 = pd.read_table(full_path_isos, header=None)[0].to_list()

# Compute distribution and prevalence: 266 KOs
res_90_266 = divide_and_organize_data(kos, kos266, isos_90_266, taxonomy)

# Exact test
barnard_90_266 = test_group_proportions(res_90_266['top'], res_90_266['nontop'], 5)

# Save results
file_name_out = "Barnard_stat/barnard_90_266.tsv"
file_name_out_full = os.path.join(path, file_name_out)

barnard_90_266.to_csv(file_name_out_full, sep="\t", index=None)
