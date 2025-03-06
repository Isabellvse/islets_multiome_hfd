#!/usr/bin/bash

eval "$(conda shell.bash hook)"
conda activate /work/islet_scenic_plus/myenv

exec > >(tee -i /work/Home/islets_multiome_hfd/data/scenicplus/beta/output_last.log)
exec 2>&1

python -c "
print('Import packages');
import scenicplus;
from pycisTopic.cistopic_class import *;
from pycisTopic.clust_vis import *;
from pycisTopic.lda_models import *;
from pycisTopic.topic_binarization import *;
from pycisTopic.diff_features import *;
from pycistarget.utils import region_names_to_coordinates;
from pycistarget.utils import region_names_to_coordinates;
from scenicplus.scenicplus_class import create_SCENICPLUS_object;
import pandas as pd;
import os;
import pickle;
import scanpy as sc;
import pyranges as pr;
import dill;
import warnings;
import numpy as np;
from scenicplus.wrappers.run_pycistarget import run_pycistarget;
warnings.filterwarnings('ignore');
# Set stderr to null to avoid strange messages from ray
import sys;
_stderr = sys.stderr;
null = open(os.devnull,'wb');

print('setting up analysis directories and common variables');
baseDir = '/work/islet_scenic_plus/';
	
# project directory
projDir = '/work/Home/islets_multiome_hfd/data/scenicplus/beta/';


# output directory
outDir = projDir + 'results/';
if not os.path.exists(outDir):
    os.makedirs(outDir);
    
# temporary directory
tempDir = projDir + 'temp/';
if not os.path.exists(tempDir):
    os.makedirs(tempDir);

malletDir = tempDir + 'mallet/';
if not os.path.exists(malletDir):
    os.makedirs(malletDir);
    
# figure
figureDir = outDir + 'figures/';
if not os.path.exists(figureDir):
    os.makedirs(figureDir);

motifsDir = outDir + 'motifs/';
if not os.path.exists(motifsDir):
    os.makedirs(motifsDir);

scatacDir = outDir + 'scatac/';
if not os.path.exists(scatacDir):
    os.makedirs(scatacDir);

scenicplusDir = outDir + 'scenicplus/';
if not os.path.exists(scenicplusDir):
    os.makedirs(scenicplusDir);

print('the paths are');
print('Base directory:', baseDir);
print('Project directory:', projDir);
print('Temprorary directory:', tempDir);
print('Output directory:', outDir);
print('Figure directory:', figureDir);
print('Motif directory:', motifsDir);
print('Scatac directory:', scatacDir);


print('define colors for plotting');
color_dict = {
    'Beta':'#A83708',
    'Alpha':'#2e8b57',
    'Delta':'#004B7A',
    'Gamma':'#F8AD4B',
    'Acinar':'#1e90ff',
    'Stellate':'#FADEB8',
    'Proliferating':'#ff4500',
    'Endothelial':'#ba55d3',
    'Immune':'#e6c2dc'
}

sample_color = {
    'LFD_R1':'#FABA3E',
    'LFD_R2':'#1B8235',
    'LFD_R3':'#1977FA',
    'HFD_1_R1':'#FA8231',
    'HFD_1_R2':'#004B7A',
    'HFD_3_R1':'#8F426D',
    'HFD_3_R2':'#FA3C25'
}

condition_color = {
    'LFD':'#FA8231',
    'HFD_1':'#46657A',
    'HFD_3':'#A86041'
}


print('import data');

# count matrix
count_matrix = pd.read_csv(projDir+'files/counts_atac.csv', sep=',', index_col=0, engine='pyarrow');

# cell data
cell_data = pd.read_csv(projDir+'files/atac_meta.csv', sep='\t', index_col=0);

# path to blacklist
path_to_blacklist=projDir+'files/mm10-blacklist.v2.bed';

print('create cistopic object');

# create Cistopic object
cistopic_obj = create_cistopic_object(fragment_matrix=count_matrix, path_to_blacklist=path_to_blacklist);

# Adding cell information
cistopic_obj.add_cell_data(cell_data);

pickle.dump(cistopic_obj,
            open(os.path.join(outDir, 'cistopic_obj.pkl'), 'wb'));
print(cistopic_obj);

print('topic modelling');
from pycisTopic.lda_models import *
path_to_mallet_binary='/work/islet_scenic_plus/Mallet/bin/mallet'

models=run_cgs_models_mallet(cistopic_obj = cistopic_obj,
                    n_topics=[2,5,10,15,20,25,30,35,40,45,50],
                    path_to_mallet_binary = path_to_mallet_binary,
                    n_cpu=1,
                    n_iter=500, 
                    random_state=555,
                    alpha=50,
                    alpha_by_topic=True,
                    eta=0.1,
                    eta_by_topic=False,
                   tmp_path = malletDir,
                            save_path = malletDir)
print('save');
pickle.dump(models,
            open(os.path.join(outDir, 'models'), 'wb'));


print('evaluate model');
model = evaluate_models(models,
                       return_model=True,
                       metrics=['Arun_2010','Cao_Juan_2009', 'Minmo_2011', 'loglikelihood'],
                       plot_metrics=False,
                       save= figureDir + 'model_selection.pdf');

print('add model to cistopic object');

# add model
cistopic_obj.add_LDA_model(model);

print('export topic contributions');
topic_contributions = cistopic_obj.selected_model.cell_topic;
topic_contributions.to_csv(outDir + 'cistopic_contributions.csv');

# save
pickle.dump(cistopic_obj,
            open(os.path.join(outDir, 'cistopic_obj.pkl'), 'wb'));
			
print('script 2');
print('run umap');

# run umap
run_umap(cistopic_obj, target  = 'cell', scale = True);

# plot umap
plot_metadata(cistopic_obj, 
              reduction_name = 'UMAP', 
              variables = ['manual_anno', 'sample_id', 'condition'],
              num_columns=3,
                 text_size=10,
                 dot_size=3,
                 figsize=(10,5),
                 color_dictionary={'sample_id': sample_color,
                                'manual_anno': color_dict,
                               'condition': condition_color},
                 save= figureDir + 'dimensionality_reduction_label_uncorrected.pdf');

print('run harmony');
# Harmony
harmony(cistopic_obj, 'sample_id', random_state=555);

# UMAP
run_umap(cistopic_obj, 
         reduction_name='harmony_UMAP',
         target  = 'cell', 
         harmony=True);

plot_metadata(cistopic_obj, 
              reduction_name = 'harmony_UMAP', 
              variables = ['manual_anno', 'sample_id', 'condition'],
              num_columns=3,
                 text_size=10,
                 dot_size=1,
                 figsize=(10,5),
                 color_dictionary={'sample_id': sample_color,
                                'manual_anno': color_dict,
                               'condition': condition_color},
                 save= figureDir + 'dimensionality_reduction_label_corrected.pdf');

print('binirization');
region_bin_topics_otsu = binarize_topics(cistopic_obj, method='otsu');
region_bin_topics_top3k = binarize_topics(cistopic_obj, method='ntop', ntop = 3000);

print('save');

if not os.path.exists(os.path.join(scatacDir, 'candidate_enhancers')):
    os.makedirs(os.path.join(scatacDir, 'candidate_enhancers'));

pickle.dump(region_bin_topics_otsu, open(os.path.join(scatacDir, 'candidate_enhancers/region_bin_topics_otsu.pkl'), 'wb'));
pickle.dump(region_bin_topics_top3k, open(os.path.join(scatacDir, 'candidate_enhancers/region_bin_topics_top3k.pkl'), 'wb'));

print('DARs across cell type');
# impute the region accessibility exploting the cell-topic and topic-region probabilities. 
# To shrink very low probability values to 0, we use a scale factor (by default: 10^6).
imputed_acc_obj = impute_accessibility(cistopic_obj, 
                                       selected_cells=None, 
                                       selected_regions=None, 
                                       scale_factor=10**6);
# log-normalize the imputed data.
normalized_imputed_acc_obj = normalize_scores(imputed_acc_obj, 
                                              scale_factor=10**4);
                                              
# Identify highly variable regions. This is not mandatory, but will speed up the hypothesis testing step for identifying DARs.
variable_regions = find_highly_variable_features(normalized_imputed_acc_obj, 
                                                 plot = False);
# Identify differentially accessible regions between groups. 
# By default, this function will perform a Wilcoxon rank-sum test between each group in the specified variable and the rest. 
markers_dict = find_diff_features(cistopic_obj, 
                                  imputed_acc_obj, 
                                  variable='condition', 
                                  var_features=variable_regions, 
                                  contrasts= [[['LFD'], ['HFD_1']], 
                                   [['LFD'], ['HFD_3']], 
                                   [['LFD'], ['HFD_1', 'HFD_3']],
                                   [['HFD_1'], ['LFD']],
                                   [['HFD_1'], ['HFD_3']],
                                   [['HFD_1'], ['LFD', 'HFD_3']],
                                   [['HFD_3'], ['LFD']],
                                   [['HFD_3'], ['HFD_1']],
                                   [['HFD_3'], ['LFD', 'HFD_1']]],
                                  split_pattern = '-');
print('save');
pickle.dump(markers_dict, open(os.path.join(scatacDir, 'candidate_enhancers/markers_dict.pkl'), 'wb'));

print('script 3');
print('import data');
region_bin_topics_otsu = pickle.load(open(os.path.join(scatacDir, 'candidate_enhancers/region_bin_topics_otsu.pkl'), 'rb'));
region_bin_topics_top3k = pickle.load(open(os.path.join(scatacDir, 'candidate_enhancers/region_bin_topics_top3k.pkl'), 'rb'));
markers_dict = pickle.load(open(os.path.join(scatacDir, 'candidate_enhancers/markers_dict.pkl'), 'rb'));

print('Convert to dictionary of pyranges objects.');

region_sets = {};
region_sets['topics_otsu'] = {};
region_sets['topics_top_3'] = {};
region_sets['DARs'] = {};
for topic in region_bin_topics_otsu.keys():
    regions = region_bin_topics_otsu[topic].index[region_bin_topics_otsu[topic].index.str.startswith('chr')] #only keep regions on known chromosomes
    region_sets['topics_otsu'][topic] = pr.PyRanges(region_names_to_coordinates(regions));
for topic in region_bin_topics_top3k.keys():
    regions = region_bin_topics_top3k[topic].index[region_bin_topics_top3k[topic].index.str.startswith('chr')] #only keep regions on known chromosomes
    region_sets['topics_top_3'][topic] = pr.PyRanges(region_names_to_coordinates(regions));
for DAR in markers_dict.keys():
    regions = markers_dict[DAR].index[markers_dict[DAR].index.str.startswith('chr')] #only keep regions on known chromosomes
    region_sets['DARs'][DAR] = pr.PyRanges(region_names_to_coordinates(regions));

for key in region_sets.keys():
    print(f'{key}: {region_sets[key].keys()}');

print('Define rankings, score and motif annotation database.');

rankings_db = os.path.join(projDir, 'database/mm10_screen_v10_clust.regions_vs_motifs.rankings.feather');
scores_db = os.path.join(projDir, 'database/mm10_screen_v10_clust.regions_vs_motifs.scores.feather');
motif_annotation = os.path.join(projDir, 'database/motifs-v10-nr.mgi-m0.00001-o0.0.tbl');
biomart_host = 'http://nov2020.archive.ensembl.org/';

run_pycistarget(
    region_sets = region_sets,
    species = 'mus_musculus',
    save_path = motifsDir,
    ctx_db_path = rankings_db,
    dem_db_path = scores_db,
    biomart_host = biomart_host,
    path_to_motif_annotations = motif_annotation,
    run_without_promoters = False,
    n_cpu = 15,
    annotation_version = 'v10nr_clust',
    );

print('script 3');
print('import data');
adata = sc.read_h5ad(os.path.join(projDir, 'data/rna.h5ad')); 
print(adata)
cistopic_obj = dill.load(open(os.path.join(outDir, 'cistopic_obj.pkl'), 'rb'));
print(cistopic_obj)
menr = dill.load(open(os.path.join(motifsDir, 'menr.pkl'), 'rb'));

print('create object');
scplus_obj = create_SCENICPLUS_object(
    GEX_anndata = adata,
    cisTopic_obj = cistopic_obj,
    menr = menr, 
    bc_transform_func = lambda x: f'{x}___cisTopic');
scplus_obj.X_EXP = np.array(scplus_obj.X_EXP.todense());
scplus_obj;

biomart_host = 'http://nov2020.archive.ensembl.org/';
list_tf_path = baseDir+'scenicplus/resources/allTFs_mm.txt';

from scenicplus.wrappers.run_scenicplus import run_scenicplus;

try:
    run_scenicplus(
        scplus_obj = scplus_obj,
        variable = ['GEX_condition'],
        species = 'mmusculus',
        assembly = 'mm10',
        tf_file = list_tf_path,
        save_path = os.path.join(scenicplusDir),
        biomart_host = biomart_host,
        upstream = [1000, 150000],
        downstream = [1000, 150000],
        calculate_TF_eGRN_correlation = True,
        calculate_DEGs_DARs = False,
        export_to_loom_file = False,
        export_to_UCSC_file = False,
        n_cpu = 15,
    _temp_dir = None);
except Exception as e:
    dill.dump(scplus_obj, open(os.path.join(scenicplusDir, 'scplus_obj.pkl'), 'wb'), protocol=-1);
    raise(e)"