# Description -------------------------------------------------------------
# Investigating cell-to-cell communication with multi nichenet
# https://github.com/saeyslab/multinichenetr/blob/main/vignettes/detailed_analysis_steps_empirical_pvalues.md


# Overview of steps -------------------------------------------------------

# 1. Extract cell type abundance and expression information from receiver and sender cell types,
#    and link this expression information for ligands of the sender cell types to the corresponding receptors of the receiver cell types

# 2. Perform genome-wide differential expression analysis of receiver and sender cell types to define DE genes between the conditions of interest.
#    Based on this analysis, we can define the logFC/p-value of ligands in senders and receptors in receivers, and define the set of affected target genes in the receiver.

# 3. Predict NicheNet ligand activities and NicheNet ligand-target links based on these differential expression results

# 4. Use the information collected above to prioritize all sender-ligand—receiver-receptor pairs.

# 5. Calculate correlation in expression between ligand-receptor pairs and their predicted target genes

# Setup -------------------------------------------------------------------
base::source(here::here("R/set_up_nichenet.R"))
set.seed(1000)
create_directories(c(here::here("data/nichenet/beta/files/"),
                     here::here("data/nichenet/beta/figures/")))
# Load --------------------------------------------------------------------
# seurat object
seurat_5 <- base::readRDS(here::here("data/seurat_objects/seurat_5.rds"))

# rename conditions
seurat_5@meta.data <-
  seurat_5@meta.data %>%
  BiocGenerics::as.data.frame() %>%
  dplyr::mutate(condition = stringr::str_replace(condition, "_", ""))


# NicheNet ligand-receptor network and ligand-target matrix
organism = "mouse"
if(organism == "human"){
  lr_network = readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
  lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% dplyr::distinct(ligand, receptor) %>% dplyr::mutate(ligand = base::make.names(ligand), receptor = base::make.names(receptor))
  ligand_target_matrix = readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"))
  BiocGenerics::colnames(ligand_target_matrix) = BiocGenerics::colnames(ligand_target_matrix) %>% base::make.names()
  BiocGenerics::rownames(ligand_target_matrix) = BiocGenerics::rownames(ligand_target_matrix) %>% base::make.names()
} else if(organism == "mouse"){
  lr_network = readRDS(url("https://zenodo.org/record/7074291/files/lr_network_mouse_21122021.rds"))
  lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% distinct(ligand, receptor) %>% mutate(ligand = base::make.names(ligand), receptor = base::make.names(receptor))
  ligand_target_matrix = readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds"))
  BiocGenerics::colnames(ligand_target_matrix) = BiocGenerics::colnames(ligand_target_matrix) %>% base::make.names()
  BiocGenerics::rownames(ligand_target_matrix) = BiocGenerics::rownames(ligand_target_matrix) %>% base::make.names()
}

# Prepare data ------------------------------------------------------------
## Convert seurat object to single-cell object ----
sce = Seurat::as.SingleCellExperiment(seurat_5, assay = "RNA")

## Covnert gene symbols to the most recent version of gene symbol ----
sce = multinichenetr::alias_to_symbol_SCE(sce, "mouse") %>%
  multinichenetr::makenames_SCE()

# Check data --------------------------------------------------------------

## Number of cell for each replicate + condition ----
base::table(SummarizedExperiment::colData(sce)$manual_anno, SummarizedExperiment::colData(sce)$'orig.ident')

# Setup analysis ----------------------------------------------------------

## Define variables ----
sample_id = "orig.ident"
group_id = "condition"
celltype_id = "manual_anno"
covariates = NA
batches = NA

## Sender and reciver cells ----
senders_oi = manual_anno_levels
receivers_oi = "Beta"
# 
# # ## Remove Acinar because there is not enough cells
# senders_oi <- senders_oi[ !senders_oi %in% c("Acinar")]
# receivers_oi <- "Beta"

## Subset object
sce = sce[, SummarizedExperiment::colData(sce)[,celltype_id] %in% senders_oi]

## Subset sce object

## Set minimum number of cells ----
min_cells = 0

## Define contrast ----
contrasts_oi = c("'LFD-(HFD1+HFD3)/2','HFD1-(LFD)/1', 'HFD3-(LFD)/1'")
contrasts_oi

## Contrast table ----
contrast_tbl = tibble(contrast = c('LFD-(HFD1+HFD3)/2','HFD1-(LFD)/1', 'HFD3-(LFD)/1'),
                      group = c("LFD","HFD1", "HFD3"))
contrast_tbl


## define parameters ----
logFC_threshold = 1
p_val_threshold = 0.05
fraction_cutoff = 0.05

## We will set p-value cutoff at the adjusted p-value ----
p_val_adj = TRUE
empirical_pval = FALSE

## Set top number of targets per ligand ----
top_n_target = 250

## Set up the system to run in parallel
cores_system = 64
n.cores = base::min(cores_system, base::union(senders_oi, receivers_oi) %>% base::length()) # use one core per receiver cell type

verbose = TRUE
# Step 1 ------------------------------------------------------------------

## Calculate abundance and expression information of each celltype ----
## And link this expression information for ligands of the sender cell types to the corresponding receptors of the receiver cell types

abundance_expression_info = multinichenetr::get_abundance_expression_info(
  sce = sce,
  sample_id = sample_id,
  group_id = group_id,
  celltype_id = celltype_id,
  min_cells = min_cells,
  senders_oi = senders_oi,
  receivers_oi = receivers_oi,
  lr_network = lr_network,
  batches = batches
)

## Plot cell abundance ----
# The first plot visualizes the number of cells per celltype-sample combination,
# and indicates which combinations are removed during the DE analysis
# because there are less than min_cells in the celltype-sample combination.

# The red dotted line indicates the required minimum of cells as defined above in min_cells
abundance_expression_info$abund_plot_sample

## plot cell abundance per group ----
# Look at Differential abundance between the conditions.
# This because the pseudobulking approach behind Muscat could potentially suffer
# from some biases if there would be huge differences in abundances of a cell type between different groups.
# Downstream results of these cell types should then be considered with some caution.
abundance_expression_info$abund_plot_group

## The numbers behind the plots ----
abundance_expression_info$abundance_data_sender

# prod is product
## Looking at genes per sample ----
# Average normalized expression value
abundance_expression_info$celltype_info$avg_df %>% head()
# Average fraction of expression cells with non-zero counts
abundance_expression_info$celltype_info$frq_df %>% head()
# logCPM-pseudocounts
abundance_expression_info$celltype_info$pb_df %>% head()

## per group ----
# Average normalized expression value
abundance_expression_info$celltype_info$avg_df_group %>% head()
# Average fraction of expression cells with non-zero counts
abundance_expression_info$celltype_info$frq_df_group %>% head()
# logCPM-pseudocounts
abundance_expression_info$celltype_info$pb_df_group %>% head()

## Ligand -receptor pair combinations - sample based ----
# Average normalized expression value
abundance_expression_info$sender_receiver_info$avg_df %>% head()
# Average fraction of expression cells with non-zero counts
abundance_expression_info$sender_receiver_info$frq_df %>% head()
# logCPM-pseudocounts
abundance_expression_info$sender_receiver_info$pb_df %>% head()

## Ligand -receptor pair combinations - group based ----
# Average normalized expression value
abundance_expression_info$sender_receiver_info$avg_df_group %>% head()
# Average fraction of expression cells with non-zero counts
abundance_expression_info$sender_receiver_info$frq_df_group %>% head()
# logCPM-pseudocounts
abundance_expression_info$sender_receiver_info$pb_df_group %>% head()



# Step2: ------------------------------------------------------------------

## Perform DE analysis for each celltype ----
DE_info = multinichenetr::get_DE_info(
  sce = sce,
  sample_id = sample_id,
  group_id = group_id,
  celltype_id = celltype_id,
  batches = batches,
  covariates = covariates,
  contrasts_oi = contrasts_oi,
  min_cells = min_cells,
  de_method_oi = "DESeq2")

## Check results ----
DE_info$celltype_de$de_output_tidy %>%
  dplyr::arrange(p_val) %>%
  utils::head()

## Distribution of p-values ----
# It is always a good idea to check distribution of the p-values resulting from this DE expression analysis.
# In order to trust the p-values, the p-value distributions should be uniform distributions,
# with a peak allowed between 0 and 0.05 if there would be a clear biological effect in the data.

DE_info$hist_pvals

# continue with this https://github.com/saeyslab/multinichenetr/blob/main/vignettes/detailed_analysis_steps_MISC.md
if(empirical_pval == FALSE){
  celltype_de = DE_info$celltype_de$de_output_tidy
} else {
  celltype_de = DE_info_emp$de_output_tidy_emp %>% dplyr::select(-p_val, -p_adj) %>% dplyr::rename(p_val = p_emp, p_adj = p_adj_emp)
}

## Combine DE information for ligand-senders and receptors-receivers
sender_receiver_de = multinichenetr::combine_sender_receiver_de(
  sender_de = celltype_de,
  receiver_de = celltype_de,
  senders_oi = senders_oi,
  receivers_oi = receivers_oi,
  lr_network = lr_network
)


# Step 3. -----------------------------------------------------------------

## Run ligand activity analysis
ligand_activities_targets_DEgenes = multinichenetr::get_ligand_activities_targets_DEgenes(
  receiver_de = celltype_de,
  receivers_oi = receivers_oi,
  ligand_target_matrix = ligand_target_matrix,
  logFC_threshold = logFC_threshold,
  p_val_threshold = p_val_threshold,
  p_val_adj = p_val_adj,
  top_n_target = top_n_target,
  verbose = verbose,
  n.cores = n.cores
)

## Check DE genes for analysis
ligand_activities_targets_DEgenes$de_genes_df %>% head(20)

## Check output of activity analysis
ligand_activities_targets_DEgenes$ligand_activities %>% head(20)


# Step 4 ------------------------------------------------------------------
# In the 3 previous steps, we calculated expression, differential expression and NicheNet activity information.
# Now we will combine these different types of information in one prioritization scheme.

# MultiNicheNet allows the user to define the weights of the following criteria to prioritize ligand-receptor interactions:

# 1. Upregulation of the ligand in a sender cell type and/or upregulation of the receptor in a receiver cell type - in the condition of interest. : de_ligand and de_receptor
# 2. Sufficiently high expression levels of ligand and receptor in many samples of the same group (to mitigate the influence of outlier samples). : frac_exprs_ligand_receptor
# 3. Cell-type and condition specific expression of the ligand in the sender cell type and receptor in the receiver cell type (to mitigate the influence of upregulated but still relatively weakly expressed ligands/receptors) : exprs_ligand and exprs_receptor
# 4. High NicheNet ligand activity, to further prioritize ligand-receptor pairs based on their predicted effect of the ligand-receptor interaction on the gene expression in the receiver cell type : activity_scaled
# 5. High relative abundance of sender and/or receiver in the condition of interest: abund_sender and abund_receiver (experimental feature - not recommended to give non-zero weights for default analyses)
# The different properties of the sender-ligand—receiver-receptor pairs can be weighted according to the user’s preference and insight in the dataset at hand.


## Define the prioritization weights, and prepare grouping objects ----
prioritizing_weights_DE = c("de_ligand" = 1,
                            "de_receptor" = 1)
prioritizing_weights_activity = c("activity_scaled" = 2)

prioritizing_weights_expression_specificity = c("exprs_ligand" = 2,
                                                "exprs_receptor" = 2)

prioritizing_weights_expression_sufficiency = c("frac_exprs_ligand_receptor" = 1)

prioritizing_weights_relative_abundance = c( "abund_sender" = 0,
                                             "abund_receiver" = 0)

prioritizing_weights = c(prioritizing_weights_DE,
                         prioritizing_weights_activity,
                         prioritizing_weights_expression_specificity,
                         prioritizing_weights_expression_sufficiency,
                         prioritizing_weights_relative_abundance)


sum_prioritization_weights = 2*prioritizing_weights["de_ligand"] +
  2*prioritizing_weights["de_receptor"] +
  prioritizing_weights["activity_scaled"] +
  prioritizing_weights["exprs_ligand"] +
  prioritizing_weights["exprs_receptor"] +
  prioritizing_weights["frac_exprs_ligand_receptor"] +
  prioritizing_weights["abund_sender"] +
  prioritizing_weights["abund_receiver"]


## Make a grouping tbl
sender_receiver_tbl = sender_receiver_de %>%
  dplyr::distinct(sender, receiver)

metadata_combined = SummarizedExperiment::colData(sce) %>%
  tibble::as_tibble()

if(!is.na(batches)){
  grouping_tbl = metadata_combined[,c(sample_id, group_id, batches)] %>% tibble::as_tibble() %>% dplyr::distinct()
  colnames(grouping_tbl) = c("sample","group",batches)
} else {
  grouping_tbl = metadata_combined[,c(sample_id, group_id)] %>% tibble::as_tibble() %>% dplyr::distinct()
  colnames(grouping_tbl) = c("sample","group")
}

## Run the prioritization ----
prioritization_tables = multinichenetr::generate_prioritization_tables(
  sender_receiver_info = abundance_expression_info$sender_receiver_info,
  sender_receiver_de = sender_receiver_de,
  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
  contrast_tbl = contrast_tbl,
  sender_receiver_tbl = sender_receiver_tbl,
  grouping_tbl = grouping_tbl,
  prioritizing_weights = prioritizing_weights,
  fraction_cutoff = fraction_cutoff,
  abundance_data_receiver = abundance_expression_info$abundance_data_receiver,
  abundance_data_sender = abundance_expression_info$abundance_data_sender
)

## Look at table ----

# group based summary
prioritization_tables$group_prioritization_tbl %>% head(20)

# sample-based summary table: contains expression information of each LR pair per sample
prioritization_tables$sample_prioritization_tbl %>% head(20)


# Step 5 ------------------------------------------------------------------
# Add information on prior knowledge and expression correlation between LR and target expression.
# In multi-sample datasets, we have the opportunity to look whether expression of ligand-receptor
# across all samples is correlated with the expression of their by NicheNet predicted target genes.
# This is what we will do with the following line of code:

# Calculate the pearson and spearman expression correlation between ligand-receptor pair pseudobulk expression products
# and DE gene expression product.
# Add the NicheNet ligand-target regulatory potential score as well as prior knowledge support for the LR-->Target link.
lr_target_prior_cor = multinichenetr::lr_target_prior_cor_inference(
  prioritization_tables$group_prioritization_tbl$receiver %>% unique(),
  abundance_expression_info,
  celltype_de,
  grouping_tbl,
  prioritization_tables,
  ligand_target_matrix,
  logFC_threshold = logFC_threshold,
  p_val_threshold = p_val_threshold,
  p_val_adj = p_val_adj
)


# Combine all data --------------------------------------------------------

multinichenet_output = base::list(
  celltype_info = abundance_expression_info$celltype_info,
  celltype_de = celltype_de,
  sender_receiver_info = abundance_expression_info$sender_receiver_info,
  sender_receiver_de =  sender_receiver_de,
  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
  prioritization_tables = prioritization_tables,
  grouping_tbl = grouping_tbl,
  lr_target_prior_cor = lr_target_prior_cor
)

multinichenet_output = multinichenetr::make_lite_output(multinichenet_output)


# Save --------------------------------------------------------------------
base::saveRDS(multinichenet_output, here::here("data/nichenet/beta/files/beta_vs_all_multinichenet.rds"))

multinichenet_output <- base::readRDS(here::here("data/nichenet/beta/files/beta_vs_all_multinichenet.rds"))

# Make Circos plot of top 30 links ----------------------------------------
colors_sender = cluster_anno
colors_receiver = cluster_anno

pdf(file = here::here("data/nichenet/beta/figures/top30_circosplot_beta_vs_all.pdf"),
    height = 3,
    width = 3)
prioritized_tbl_oi_LFD_30 = multinichenet_output$prioritization_tables$group_prioritization_tbl %>% 
  dplyr::filter(group == "LFD" & fraction_expressing_ligand_receptor > 0) %>% 
  dplyr::distinct(group, sender, receiver, ligand, receptor, receiver, id, prioritization_score) %>% 
  dplyr::mutate(prioritization_rank = rank(-prioritization_score)) %>% 
  dplyr::filter(prioritization_rank <= 30)

make_circos_one_group(prioritized_tbl_oi_LFD_30, colors_sender, colors_receiver)

prioritized_tbl_oi_HFD1_30 = multinichenet_output$prioritization_tables$group_prioritization_tbl %>% 
  dplyr::filter(group == "HFD1" & fraction_expressing_ligand_receptor > 0) %>% 
  dplyr::distinct(group, sender, receiver, ligand, receptor, receiver, id, prioritization_score) %>% 
  dplyr::mutate(prioritization_rank = rank(-prioritization_score)) %>% 
  dplyr::filter(prioritization_rank <= 30)

make_circos_one_group(prioritized_tbl_oi_HFD1_30, colors_sender, colors_receiver)

prioritized_tbl_oi_HFD3_30 = multinichenet_output$prioritization_tables$group_prioritization_tbl %>% 
  dplyr::filter(group == "HFD3" & fraction_expressing_ligand_receptor > 0) %>% 
  dplyr::distinct(group, sender, receiver, ligand, receptor, receiver, id, prioritization_score) %>% 
  dplyr::mutate(prioritization_rank = rank(-prioritization_score)) %>% 
  dplyr::filter(prioritization_rank <= 30)

make_circos_one_group(prioritized_tbl_oi_HFD3_30, colors_sender, colors_receiver)

dev.off()


# Target genes TNF & IFNB1 ------------------------------------------------
# Before, we had calculated the correlation in expression between ligand-receptor pairs and DE genes.
# Now we will filter out correlated ligand-receptor –> target links that both show high expression
# correlation (spearman or activity > 0.50 in this example)
# and have some prior knowledge to support their link.' (lr_target_prior_cor)

## HFD 1 & HFD 3 TNF + IFNB1
group_oi = c("HFD1", "HFD3")
receiver_oi = "Beta"
ligand_oi = c("Tnf", "Ifnb1")
sender_oi = "Immune"

# get top 30 ligand and receptor links for tnf and ifnb1
hfd1_tnf_ifnb1 <- prioritized_tbl_oi_HFD1_30 %>% 
  dplyr::filter(ligand == "Tnf" | ligand == "Ifnb1") %>% 
  dplyr::pull("id")

hfd3_tnf_ifnb1 <- prioritized_tbl_oi_HFD3_30 %>% 
  dplyr::filter(ligand == "Tnf" | ligand == "Ifnb1") %>% 
  dplyr::pull("id")


# Extract information about ligand of interest and group of interest
lr_target_prior_cor_filtered = multinichenet_output$lr_target_prior_cor %>%
  dplyr::inner_join(multinichenet_output$ligand_activities_targets_DEgenes$ligand_activities %>%
                      dplyr:: distinct(ligand, target, direction_regulation, contrast)) %>%
  dplyr::inner_join(contrast_tbl) %>%
  dplyr::filter(group %in% group_oi,
                receiver %in% receiver_oi,
                sender %in% sender_oi,
                ligand %in% ligand_oi,
                id %in% c(hfd1_tnf_ifnb1, hfd3_tnf_ifnb1))

# Get upregulated target genes, which have a correlation above 0.5 between target gene expression and
# the product of ligand and recepter expression
lr_target_prior_cor_filtered_up = lr_target_prior_cor_filtered %>%
  dplyr::filter(direction_regulation == "up") %>%
  dplyr::filter( (rank_of_target < top_n_target) & (pearson > 0.50 | spearman > 0.50))

# Down regulated genes
lr_target_prior_cor_filtered_down = lr_target_prior_cor_filtered %>%
  dplyr::filter(direction_regulation == "down") %>%
  dplyr::filter( (rank_of_target < top_n_target) & (pearson < -0.50 | spearman < -0.50)) # downregulation -- negative correlation

# combine dataframes
lr_target_prior_cor_filtered = dplyr::bind_rows(lr_target_prior_cor_filtered_up, lr_target_prior_cor_filtered_down)

# also save information about receptor - ligand product as well as pseudobulk expression of target genes
sample_data = multinichenet_output$prioritization_tables$sample_prioritization_tbl %>%
  dplyr::filter(id %in% lr_target_prior_cor_filtered$id) %>%
  dplyr::select(- keep_receiver,
                - keep_sender,
                - keep_sender_receiver)

# pseudobulk expression of target genes
pd <- multinichenet_output$celltype_info$pb_df %>%
  dplyr::filter(celltype == "Beta") %>%
  dplyr::rename(target = gene) %>%
  dplyr::select(-celltype)

# combine sample data and pseudobulk
sample_data <- dplyr::left_join(x = pd,
                                y = sample_data,
                                by = c("sample"))

# extract all genes tested
universe <- unique(multinichenet_output$celltype_de$gene)
# save --------------------------------------------------------------------
saveRDS(lr_target_prior_cor_filtered, here::here("data/nichenet/beta/files/hfd1_hfd3_tnf_ifnb1_target_genes.rds"))
saveRDS(sample_data, here::here("data/nichenet/beta/files/hfd1_hfd3_tnf_ifnb1_sample_data.rds"))
saveRDS(universe, here::here("data/nichenet/beta/files/nichenet_universe.rds"))


# Plot ligand activities --------------------------------------------------
norm_counts <- readRDS(here::here("data/deseq2/rna_clusterbycond/files/norm_pseudobulk_counts_per_cluster.rds"))
beta <- norm_counts[["Beta"]]

keep <- lr_target_prior_cor_filtered %>%
  dplyr::filter(direction_regulation == "up")

keep_2 <- keep$target %>% unique()

myCol <- colorRampPalette(c('#004B7A', 'white', '#A83708'))(100)
myBreaks <- seq(-1.5, 1.5, length.out = 100)

pdf(here::here("data/nichenet/beta/figures/nichenet_expression_of_target_genes_pseudobulk.pdf"),
    width = 3,
    height = 4) %>% 
beta %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  dplyr::filter(gene %in% keep_2) %>%
  tibble::column_to_rownames("gene") %>%
  as.matrix() %>%
  pheatmap::pheatmap(scale = "row", cluster_rows = TRUE,
                     cluster_cols = FALSE, 
                     color = myCol,
                     border_color = NA,
                     breaks = myBreaks, fontsize_row = 2)
dev.off()


# save prioritization table -----------------------------------------------
openxlsx::write.xlsx(multinichenet_output$prioritization_tables$group_prioritization_tbl, here::here("data/nichenet/beta/files/ligang_receptor_interactions_beta_vs_other.xlsx"))
openxlsx::write.xlsx(lr_target_prior_cor_filtered, here::here("data/nichenet/beta/files/hfd1_hfd3_tnf_ifnb1_target_genes.xlsx"))
