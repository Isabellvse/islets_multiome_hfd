#!/usr/bin/bash

eval "$(conda shell.bash hook)"
conda activate /work/islet_scenic_plus/myenv

exec > >(tee -i /work/Home/islets_multiome_hfd/data/scenicplus/beta/output_export.log)
exec 2>&1

# Run the Python script
python -c "
print('Importing necessary packages...');
import scenicplus
from pycisTopic.cistopic_class import *
from pycisTopic.clust_vis import *
from pycisTopic.lda_models import *
from pycisTopic.topic_binarization import *
from pycisTopic.diff_features import *
from pycistarget.utils import region_names_to_coordinates
from scenicplus.scenicplus_class import create_SCENICPLUS_object
import pandas as pd
import os
import pickle
import scanpy as sc
import pyranges as pr
import dill
import warnings
import numpy as np
from scenicplus.wrappers.run_pycistarget import run_pycistarget
from scenicplus.preprocessing.filtering import apply_std_filtering_to_eRegulons
from scenicplus.eregulon_enrichment import score_eRegulons
from scenicplus.RSS import *
from pycistarget.utils import get_motifs_per_TF
import csv
from scenicplus.utils import *

# Suppress warnings
warnings.filterwarnings('ignore')

# Set stderr to null to avoid strange messages from ray
import sys
_stderr = sys.stderr
null = open(os.devnull,'wb')

print('setting up analysis directories and common variables');
baseDir = '/work/islet_scenic_plus/';
	
# project directory
projDir = '/work/Home/islets_multiome_hfd/data/scenicplus/beta/';

# output directory
outDir = projDir + 'results/';
if not os.path.exists(outDir):
    os.makedirs(outDir);
    
motifsDir = outDir + 'motifs/';
if not os.path.exists(motifsDir):
    os.makedirs(motifsDir);

scatacDir = outDir + 'scatac/';
if not os.path.exists(scatacDir):
    os.makedirs(scatacDir);

scenicplusDir = outDir + 'scenicplus/';
if not os.path.exists(scenicplusDir):
    os.makedirs(scenicplusDir);

expDir = scenicplusDir + 'export/';
if not os.path.exists(expDir):
    os.makedirs(expDir);

# Load the objects
print('Loading objects...');
try:
    scplus_obj = dill.load(open(os.path.join(scenicplusDir, 'scplus_obj.pkl'), 'rb'))
    region_ranking = dill.load(open(os.path.join(scenicplusDir, 'region_ranking.pkl'), 'rb'))
    gene_ranking = dill.load(open(os.path.join(scenicplusDir, 'gene_ranking.pkl'), 'rb'))
except Exception as e:
    print(f'Error loading objects: {e}')
    raise e

# Save raw metadata
# print('Saving raw metadata...');
eRegulon_meta = scplus_obj.uns['eRegulon_metadata']
eRegulon_meta.to_csv(expDir + 'eRegulon_meta_raw.csv')

# Simplify object
print('Simplifying object...');
apply_std_filtering_to_eRegulons(scplus_obj)

# Save simplified metadata
print('Saving simplified metadata...');
eRegulon_meta = scplus_obj.uns['eRegulon_metadata_filtered']
eRegulon_meta.to_csv(expDir + 'eRegulon_metadata_filtered.csv')

# Calculate enrichment scores (AUC) for target genes and regions
print('Calculating eRegulon enrichment scores...');
score_eRegulons(scplus_obj=scplus_obj, ranking=region_ranking, eRegulon_signatures_key='eRegulon_signatures_filtered', key_added='eRegulon_AUC_filtered', enrichment_type='region', auc_threshold=0.05, normalize=False, n_cpu=64)
score_eRegulons(scplus_obj=scplus_obj, ranking=gene_ranking, eRegulon_signatures_key='eRegulon_signatures_filtered', key_added='eRegulon_AUC_filtered', enrichment_type='gene', auc_threshold=0.05, normalize=False, n_cpu=64)

auc_region = scplus_obj.uns['eRegulon_AUC_filtered']['Region_based']
auc_gene = scplus_obj.uns['eRegulon_AUC_filtered']['Gene_based']

auc_region.to_csv(expDir + 'auc_region.csv')
auc_gene.to_csv(expDir + 'auc_gene.csv')

# Save the simplified scenic plus object
print('Saving the simplified SCENICPlus object...');
dill.dump(scplus_obj, open(os.path.join(scenicplusDir, 'scplus_obj_2.pkl'), 'wb'), protocol=-1)

# Load again
scplus_obj = dill.load(open(os.path.join(scenicplusDir, 'scplus_obj_2.pkl'), 'rb'))

# Save motifs associated with each TF
print('Saving motifs associated with each TF...');
def get_motifs_for_TF(scplus_obj, TF):
    annotations = []
    for k in scplus_obj.menr.keys():
        if isinstance(scplus_obj.menr[k], dict):
            for subkey in scplus_obj.menr[k].keys():
                annotations.extend(flatten_list([get_motifs_per_TF(scplus_obj.menr[k][subkey].motif_enrichment, TF, 'Index')]))
        else:
            for subkey in scplus_obj.menr[k].motif_enrichment.keys():
                annotations.extend(flatten_list([get_motifs_per_TF(scplus_obj.menr[k].motif_enrichment[subkey], TF, 'Index')]))
    return set(flatten_list([i.split(', ') if isinstance(i, str) else [i] for i in annotations]))

inputs = scplus_obj.uns['eRegulon_metadata_filtered']['TF'].unique()

# Write to CSV
with open(expDir + 'motifs_associated_with_eregulon_TFs.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    for input in inputs:
        output = get_motifs_for_TF(scplus_obj, TF=input)
        writer.writerow([input, output])

print('Done!')
"

