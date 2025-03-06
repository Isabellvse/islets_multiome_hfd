# Description -------------------------------------------------------------
# Keep high quality eregulons - filtering using TF gene expression and motif activity (chomrvar)

# Setup -------------------------------------------------------------------
base::source(here::here("R/set_up.R"))
set.seed(1000)

# Load --------------------------------------------------------------------
# eregulons
eregulon <- read.csv(here::here("data/scenicplus/beta/results/scenicplus/export/eRegulon_metadata_filtered.csv"))

# motif annotation
motif_anno <- RcisTarget::importAnnotations(here::here("data/pwms/v10nr_clust_public/snapshots/motifs-v10-nr.mgi-m0.00001-o0.0.tbl"))

# beta cells
meta <- read.csv(here::here("data/scenicplus/beta/files/atac_meta.csv"), sep = "\t")

# AUC genes
auc_gene <- read.csv(here::here("data/scenicplus/beta/results/scenicplus/export/auc_gene.csv"),
                     check.names=FALSE) %>%
  dplyr::mutate(Cell = str_replace(Cell, "___cisTopic", "")) %>%
  dplyr::mutate(Cell = str_replace(Cell, "___cisTopic", "")) %>%
  tibble::column_to_rownames("Cell") %>%
  base::t() %>%
  BiocGenerics::as.data.frame() %>%
  # add eregulon and tf column
  tibble::rownames_to_column("eregulon") %>%
  dplyr::mutate(eregulon = base::sub("_\\(\\d+[rg]\\)$", "", eregulon),
                eregulon = stringr::str_replace(eregulon, "_extended_", "_"),
                tf = stringr::str_split(eregulon, "_", simplify = TRUE)[,1]) %>%
  dplyr::relocate(tf, .after = "eregulon")

# AUC region
auc_region <- read.csv(here::here("data/scenicplus/beta/results/scenicplus/export/auc_region.csv"),
                       check.names=FALSE) %>%
  dplyr::mutate(Cell = str_replace(Cell, "___cisTopic", "")) %>%
  tibble::column_to_rownames("Cell") %>%
  base::t() %>%
  BiocGenerics::as.data.frame() %>%
  # add eregulon and tf column
  tibble::rownames_to_column("eregulon") %>%
  dplyr::mutate(eregulon = sub("_\\(\\d+[rg]\\)$", "", eregulon),
                eregulon = stringr::str_replace(eregulon, "_extended_", "_"),
                tf = str_split(eregulon, "_", simplify = TRUE)[,1]) %>%
  dplyr::relocate(tf, .after = "eregulon")


# motifs associated with each eregulon
eregu_motifs <- read.csv(here::here("data/scenicplus/beta/results/scenicplus/export/motifs_associated_with_eregulon_TFs.csv"),
                         header = FALSE,
                         row.names = 1,
                         col.names = c("genes", "motifs"))

# Seurat object
seurat_5 <- base::readRDS(here::here("data/seurat_objects/seurat_5.rds"))

# Chromvar deviations ## look into this
dev_motifs_beta <- base::readRDS(here::here("data/scenicplus/beta/results/chromvar_deviations_motifs_beta.rds"))


# Prepare data ------------------------------------------------------------

## Motifs associated with each eregulon ----
# Remove the { and } signs
eregu_motifs$motifs <- base::gsub("\\{|\\}", "", eregu_motifs$motifs)

# Remove the single quotes
eregu_motifs$motifs <- base::gsub("'", "", eregu_motifs$motifs)

# Split the column into multiple columns - tfs with only one motifs will be repeated, will be removed later
eregu_motifs <- BiocGenerics::cbind(eregu_motifs, do.call('rbind', strsplit(eregu_motifs$motifs, ', ')))

# Remove the original column
eregu_motifs$motifs <- NULL

# Make dataframe into long format
eregu_motifs <- eregu_motifs %>%
  tibble::rownames_to_column("tf") %>%
  tidyr::pivot_longer(cols = -tf, names_to = "delete", values_to = "motifs") %>%
  dplyr::select(-delete) %>%
  distinct()   # because tfs with 1 motif were repeated

## Save
base::saveRDS(eregu_motifs, here::here("data/scenicplus/beta/results/eregulon_motifs.rds"))

## Transcription factor gene expression ----
# Transcription factors with eRegulons
tfs <- unique(eregu_motifs$tf)

# Extract normalized counts for each TF in the seurat object
rna_exp_tfs <- seurat_5 %>%
  Seurat::GetAssayData(assay = "RNA",
                       slot = "data") %>%
  as.data.frame() %>%
  dplyr::select(rownames(meta)) %>%
  base::t() %>%
  BiocGenerics::as.data.frame() %>%
  dplyr::select(all_of(tfs)) %>%
  base::t() %>%
  BiocGenerics::as.data.frame()


# Calculate pseudobulk - 100 pseudocells ----------------------------------

## Pseudobulk generation ----

# add eregulons to row names
BiocGenerics::rownames(auc_region) <- auc_region$eregulon
BiocGenerics::rownames(auc_gene) <- auc_gene$eregulon

# generate pseudobulk
pseudo_counts <- generate_pseudobulks(counts_gene = base::t(rna_exp_tfs),
                                      auc_region = base::t(dplyr::select(auc_region, -eregulon, -tf)),
                                      auc_gene = base::t(dplyr::select(auc_gene, -eregulon, -tf)),
                                      z_scores = base::t(dev_motifs_beta),
                                      meta = meta,
                                      variable = "condition",
                                      seed = 555,
                                      nr_cells = 10,
                                      nr_pseudobulks = 100)

# save --------------------------------------------------------------------
base::saveRDS(pseudo_counts, here::here("data/scenicplus/beta/results/pseudo_counts.rds"))

# Mean pseudobulk z-score deviations --------------------------------------

# Create a list of TFs contaning a list of each motif it is associated with that TF
tf_list <- eregu_motifs %>%
  base::split(f = as.factor(.$tf)) %>%
  purrr::modify_depth(1, ~ dplyr::select(., motifs) %>%
                        base::split(f = base::as.factor(.$motifs))) %>%
  purrr::modify_depth(2, ~ purrr::as_vector(.) %>%
                        base::unname())

# To the tf_list, all the individual motifs associated with this tf (unpacking meta-clusters, which are found in dev_motifs)
# Some TFs only have transfec_pro motifs, so they cannot be included in this analysis
tf_list <- BiocGenerics::Map(function(df, dev_motif){
  
  df_3 <- BiocGenerics::Map(function(df_2){
    
    # grep motifs from deviation matrix that contains motifs present in each eregulon
    output <- rownames(dev_motif)[grep(df_2, rownames(dev_motif))]
    return(output)
  },
  df_2 = df)
  
  output <- do.call(c, df_3) %>% unname()
  
  return(output)
  
},
df = tf_list,
dev_motif = rep(list(dev_motifs_beta), length(tf_list)))

# Calculate the mean motif activity for each tf
pseudo_mean_tf_activity <- BiocGenerics::Map(function(vec, dev_score, name){
  output <- dev_score %>%
    as.data.frame() %>%
    tibble::rownames_to_column("motif") %>%
    dplyr::filter(motif %in% vec) %>%
    tidyr::pivot_longer(-motif, names_to = "barcode", values_to = "values") %>%
    dplyr::group_by(barcode) %>%
    dplyr::summarise(mean = mean(values)) %>%
    tidyr::pivot_wider(names_from = barcode, values_from = mean) %>%
    dplyr::mutate(tf = name) %>%
    tibble::column_to_rownames("tf")
  
  return(output)
  
},
vec = tf_list,
dev_score = rep(list(pseudo_counts[["z_scores_agg"]]), length(tf_list)),
name = base::names(tf_list)) %>%
  dplyr::bind_rows()

# scale TF gene expression and motif activity -----------------------------
### gene expression ----
pseudo_tf_scale <- pseudo_counts[["rna_counts_agg"]] %>%
  base::t() %>%
  BiocGenerics::as.data.frame() %>%
  dplyr::mutate(dplyr::across(dplyr::everything(), ~ (. - base::min(.)) / (base::max(.) - base::min(.)))) %>%
  base::t() %>%
  BiocGenerics::as.data.frame() %>%
  tibble::rownames_to_column("tf") %>% # make rownames a column
  dplyr::relocate("tf")

### deviation expression ----
# scale
pseudo_dev_scale <- pseudo_mean_tf_activity %>%
  base::t() %>%
  BiocGenerics::as.data.frame() %>%
  dplyr::mutate(dplyr::across(dplyr::everything(), ~ (. - base::min(.)) / (base::max(.) - base::min(.)))) %>%
  base::t() %>%
  BiocGenerics::as.data.frame() %>%
  tibble::rownames_to_column("tf") %>% # make rownames a column
  dplyr::relocate("tf")

### Add as column to the rest of the dataframes ----
# AUC REGION
pseudo_auc_region <- pseudo_counts[["auc_region_agg"]] %>%
  BiocGenerics::as.data.frame() %>%
  tibble::rownames_to_column("eregulon") %>%
  dplyr::mutate(tf = stringr::str_split(eregulon, "_", simplify = TRUE)[,1]) %>%
  dplyr::relocate("tf")

# AUC GENE
pseudo_auc_gene <- pseudo_counts[["auc_gene_agg"]] %>%
  BiocGenerics::as.data.frame() %>%
  tibble::rownames_to_column("eregulon") %>%
  dplyr::mutate(tf = stringr::str_split(eregulon, "_", simplify = TRUE)[,1]) %>%
  dplyr::relocate("tf")


## convert to long format ----
pseudo_df_list <- BiocGenerics::Map(function(df, names, values){
  
  # columns to exclude
  col <- base::names(df)[BiocGenerics::sapply(df, base::is.character)]
  
  output <- df %>%
    tidyr::pivot_longer(-col, names_to = names, values_to = values)
  return(output)
},
df = base::list(pseudo_tf_scale, pseudo_dev_scale, pseudo_auc_region, pseudo_auc_gene),
names = base::rep(list("cell"), 4),
values = base::list("tf_exp", "tf_dev", "auc_region", "auc_gene"))

# merge dataframes, and add tf exp and motif dev sum
pseudo_merged_df <- pseudo_df_list %>% purrr::reduce(dplyr::full_join, by = c("tf", "cell")) %>%
  tidyr::drop_na() %>%
  dplyr::rename(eregulon = eregulon.x) %>%
  dplyr::select(-eregulon.y) %>% 
  dplyr::relocate(eregulon) %>%
  dplyr::mutate(product = tf_exp * tf_dev,
                condition = stringr::str_extract(cell, "LFD|HFD_1|HFD_3"),
                condition = base::factor(condition, levels = condition_levels))

# Correlation -------------------------------------------------------------
# prepare data data
df_long <- pseudo_merged_df %>% 
  dplyr::select(eregulon, tf, product, auc_region, auc_gene, condition) %>%
  tidyr::pivot_longer(c(-eregulon, -tf, -product, -condition), names_to = "auc", values_to = "value")

# Step 2: Compute correlation p-values for each group
cor_pvalues <- df_long %>%
  dplyr::group_by(eregulon, tf, auc) %>%
  dplyr::summarise(cor_test = base::list(stats::cor.test(value, product))) %>%
  dplyr::mutate(p_value = base::sapply(cor_test, function(x) x$p.value),
                r = base::sapply(cor_test, function(x) x$estimate)) %>%
  dplyr::ungroup()

# Step 3: Adjust the p-values for multiple comparisons
cor_pvalues <- cor_pvalues %>%
  dplyr::mutate(fdr = stats::p.adjust(p_value, "fdr"))  # or use "fdr" for false discovery rate

# Step 4: Merge adjusted p-values with your original data
df_long_2 <- df_long %>%
  dplyr::left_join(cor_pvalues %>% dplyr::select(eregulon, tf, auc, r, fdr), 
            by = c("tf", "auc", "eregulon"))

# Step 5: Plot with adjusted p-values
pseudo_cor_plot <- df_long_2 %>%
  ggplot2::ggplot(aes(x = value, y = product)) +
  ggplot2::geom_point(aes(color = condition), size = 0.1) +
  ggplot2::scale_color_manual(values = condition_color) +
  ggplot2::geom_text(data = cor_pvalues, 
            aes(x = Inf, y = Inf, label = base::paste("R =",format(r, digits = 2), "FDR =", format(fdr, digits = 2))),
            inherit.aes = FALSE, hjust = 1.1, vjust = 1.1, size = 1) +
  ggplot2::labs(y = "TF gene expression * TF motif activity",
                x = "Activity score (AUCell)") +
  ggplot2::facet_wrap(~ tf + auc, scales = "free", ncol = 6,  labeller = 
                        labeller(.multi_line = FALSE)) +
  ggplot2::scale_y_continuous(breaks=scales::pretty_breaks(n = 2)) +
  ggplot2::scale_x_continuous(breaks=scales::pretty_breaks(n = 2)) +
  my_theme() +
  ggplot2::theme(legend.position = "none")

pdf(file = here::here("data/scenicplus/beta/results/eregulon_correlation_semi_pseudobulk.pdf"),
    height = 8,
    width = 5)
pseudo_cor_plot
dev.off()

# Filter ------------------------------------------------------------------
enhan <- df_long_2 %>%
  dplyr::filter(auc == "auc_region") %>%
  dplyr::filter(fdr <= 0.05 & abs(r) >= 0.6)

gen <- df_long_2 %>%
  dplyr::filter(auc == "auc_gene") %>%
  dplyr::filter(fdr <= 0.05 & abs(r) >= 0.6)

high_q <- intersect(enhan$eregulon, gen$eregulon)

# plot high quality -------------------------------------------------------
pseudo_cor_plot <- df_long_2 %>%
  dplyr::filter(eregulon %in% high_q) %>% 
  ggplot2::ggplot(aes(x = value, y = product)) +
  ggplot2::geom_point(aes(color = condition), size = 0.1) +
  ggplot2::scale_color_manual(values = condition_color) +
  ggplot2::geom_text(data = cor_pvalues %>% dplyr::filter(eregulon %in% high_q), 
                     aes(x = Inf, y = Inf, label = base::paste("R =",format(r, digits = 2), "FDR =", format(fdr, digits = 2))),
                     inherit.aes = FALSE, hjust = 1.1, vjust = 1.1, size = 1) +
  ggplot2::labs(y = "TF gene expression * TF motif activity",
                x = "Activity score (AUCell)") +
  ggplot2::facet_wrap(~ tf + auc, scales = "free", ncol = 4,  labeller = 
                        labeller(.multi_line = FALSE)) +
  ggplot2::scale_y_continuous(breaks=scales::pretty_breaks(n = 2)) +
  ggplot2::scale_x_continuous(breaks=scales::pretty_breaks(n = 2)) +
  my_theme() +
  ggplot2::theme(legend.position = "bottom")

pdf(file = here::here("data/scenicplus/beta/results/eregulon_correlation_semi_pseudobulk_high_q.pdf"),
    height = 8,
    width = 6)
pseudo_cor_plot
dev.off()


# extract correlations ----------------------------------------------------
cor_pvalues <- df_long %>%
  dplyr::group_by(eregulon, tf, auc) %>%
  dplyr::summarise(cor_test = list(stats::cor.test(value, product))) %>%
  dplyr::mutate(
    Correlation_coefficient = sapply(cor_test, function(x) x$estimate),
    T_statistic = sapply(cor_test, function(x) x$statistic),
    P_value = sapply(cor_test, function(x) x$p.value),
    FDR = stats::p.adjust(P_value , "fdr"), 
    Confidence_interval_lower = sapply(cor_test, function(x) x$conf.int[1]),
    Confidence_interval_upper = sapply(cor_test, function(x) x$conf.int[2]),
    Method = sapply(cor_test, function(x) x$method),
    Alternative = sapply(cor_test, function(x) x$alternative)
  ) %>%
  dplyr::ungroup() %>% 
  dplyr::select(-cor_test)

# Save --------------------------------------------------------------------
saveRDS(pseudo_merged_df, here::here("data/scenicplus/beta/results/pseudo_counts_scaled.rds"))
saveRDS(df_long_2, here::here("data/scenicplus/beta/results/pseudo_counts_correlation.rds"))
saveRDS(cor_pvalues, here::here("data/scenicplus/beta/results/pseudo_counts_correlation_test.rds"))
saveRDS(high_q, here::here("data/scenicplus/beta/results/high_quality_eregulons.rds"))
