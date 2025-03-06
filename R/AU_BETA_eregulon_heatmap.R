# Description -------------------------------------------------------------
# High quality eregulon heatmaps

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


# Seurat object
seurat_5 <- readRDS(here::here("data/seurat_objects/seurat_5.rds"))

# Chromvar deviations
dev_motifs_beta <- base::readRDS(here::here("data/scenicplus/beta/results/chromvar_deviations_motifs_beta.rds"))

# eregu_motifs
eregu_motifs <- base::readRDS(here::here("data/scenicplus/beta/results/eregulon_motifs.rds"))

high_quality_eregulon <- base::readRDS(here::here("data/scenicplus/beta/results/high_quality_eregulons.rds"))


# Extract data ------------------------------------------------------------

## Transcription factor gene expression ----
# Transcription factors with eRegulons
tfs <- unique(eregu_motifs$tf)

# Extract normalized counts for each TF in each cell in the seurat object
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


## Deviations for each TF per cell

# Create a list of eregulongs, contaning a list of each motif it is associatated with
eregu_list <- eregu_motifs %>%
  base::split(f = as.factor(.$tf)) %>%
  purrr::modify_depth(1, ~ dplyr::select(., motifs) %>%
                        base::split(f = as.factor(.$motifs))) %>%
  purrr::modify_depth(2, ~ purrr::as_vector(.) %>%
                        unname())

# to the eregu_list, keep motifs associated with this tf
# Some TFs only have transfec_pro motifs, so they cannot be included in this analysis
eregu_motifs_list <- BiocGenerics::Map(function(df, dev_motif){
  
  df_3 <- BiocGenerics::Map(function(df_2){
    
    # grep motifs from deviation matrix that contains motifs present in each eregulon
    output <- rownames(dev_motif)[grep(df_2, rownames(dev_motif))]
    return(output)
  },
  df_2 = df)
  
  output <- do.call(c, df_3) %>% unname()
  
  return(output)
  
},
df = eregu_list,
dev_motif = rep(list(dev_motifs_beta), length(eregu_list)))

# Calculate mean motif deviations per e-regulon for each cell
mean_activity_eregu <- BiocGenerics::Map(function(vec, dev_score, name){
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
vec = eregu_motifs_list,
dev_score = rep(list(dev_motifs_beta), length(eregu_motifs_list)),
name = names(eregu_motifs_list)) %>%
  dplyr::bind_rows()

# not scaled yet
### tf expression ----
tf_scale <- rna_exp_tfs %>%
  tibble::rownames_to_column("tf")

### deviation expression ----
dev_scale <- mean_activity_eregu %>%
  tibble::rownames_to_column("tf")

## Combine data: tf expression, deviation score, region and gene auc score into 1 datframe ----

# Convert dataframes to long format
df_list <- BiocGenerics::Map(function(df, names, values){
  
  # columns to exclude
  col <- names(df)[sapply(df, is.character)]
  
  output <- df %>%
    tidyr::pivot_longer(-col, names_to = names, values_to = values)
  return(output)
},
df = list(tf_scale, dev_scale, auc_region, auc_gene),
names = rep(list("cell"), 4),
values = list("tf_exp", "tf_dev", "auc_region", "auc_gene"))

# merge dataframes, and add tf exp and motif dev sum
merged_df <- df_list %>% purrr::reduce(dplyr::full_join, by = c("tf", "cell")) %>%
  tidyr::drop_na() %>%
  dplyr::select(-eregulon.x) %>%
  dplyr::rename(eregulon = eregulon.y) %>%
  dplyr::relocate(eregulon) %>%
  dplyr::mutate(condition = str_extract(cell, "LFD|HFD_1|HFD_3"),
                condition = factor(condition, levels = condition_levels))

mean_merged_df <- merged_df %>%
  dplyr::group_by(eregulon, condition) %>%
  dplyr::summarise(mean_exp = mean(base::expm1(tf_exp)),
                   mean_dev = mean(tf_dev),
                   mean_r_auc = mean(auc_region),
                   mean_g_auc = mean(auc_gene)) %>%
  dplyr::filter(eregulon %in% high_quality_eregulon)

## Scale from 0-1

# TF expression
scale_exp <- mean_merged_df %>% dplyr::select(eregulon, condition, mean_exp) %>%
  tidyr::pivot_wider(names_from = condition, values_from = mean_exp) %>%
  tibble::column_to_rownames("eregulon") %>%
  t() %>%
  as.data.frame() %>%
  dplyr::mutate(dplyr::across(dplyr::everything(), ~ (. - min(.)) / (max(.) - min(.)))) %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("eregulon") %>%
  pivot_longer(cols = -eregulon, names_to = "condition", values_to = "scale_exp")

# deviation score
scale_dev <- mean_merged_df %>% dplyr::select(eregulon, condition, mean_dev) %>%
  tidyr::pivot_wider(names_from = condition, values_from = mean_dev) %>%
  tibble::column_to_rownames("eregulon") %>%
  t() %>%
  as.data.frame() %>%
  dplyr::mutate(dplyr::across(dplyr::everything(), ~ (. - min(.)) / (max(.) - min(.)))) %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("eregulon") %>%
  pivot_longer(cols = -eregulon, names_to = "condition", values_to = "scale_dev")

# region AUC
scale_region <- mean_merged_df %>% dplyr::select(eregulon, condition, mean_r_auc) %>%
  tidyr::pivot_wider(names_from = condition, values_from = mean_r_auc) %>%
  tibble::column_to_rownames("eregulon") %>%
  t() %>%
  as.data.frame() %>%
  dplyr::mutate(dplyr::across(dplyr::everything(), ~ (. - min(.)) / (max(.) - min(.)))) %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("eregulon") %>%
  pivot_longer(cols = -eregulon, names_to = "condition", values_to = "scale_r_auc")

# Gene AUC
scale_gene <- mean_merged_df %>% dplyr::select(eregulon, condition, mean_g_auc) %>%
  tidyr::pivot_wider(names_from = condition, values_from = mean_g_auc) %>%
  tibble::column_to_rownames("eregulon") %>%
  t() %>%
  as.data.frame() %>%
  dplyr::mutate(dplyr::across(dplyr::everything(), ~ (. - min(.)) / (max(.) - min(.)))) %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("eregulon") %>%
  pivot_longer(cols = -eregulon, names_to = "condition", values_to = "scale_g_auc")

# put into a dataframe
mean_scale_df <- mean_merged_df %>%
  dplyr::inner_join(y = scale_exp, by = c("eregulon", "condition")) %>%
  dplyr::inner_join(y = scale_dev, by = c("eregulon", "condition")) %>%
  dplyr::inner_join(y = scale_region, by = c("eregulon", "condition")) %>%
  dplyr::inner_join(y = scale_gene, by = c("eregulon", "condition")) %>%
  dplyr::mutate(product = scale_exp * scale_dev,
                condition = factor(condition, levels = condition_levels))


# Heatmap  -----------------------------------------------------------------

# Heatmap of eregulon with clusters ---------------------------------------
# create heatmap where rows are clustered based on the product
test <- mean_scale_df %>%
  dplyr::select(eregulon, condition, scale_r_auc) %>%
  tidyr::pivot_wider(id_cols = eregulon, names_from = condition, values_from = scale_r_auc) %>%
  tibble::column_to_rownames("eregulon") %>%
  pheatmap::pheatmap(cluster_cols = FALSE, scale = "none")
test
# create vector of the order observed in test
target <- test[["gtable"]][["grobs"]][[4]][["label"]]

myCol <- colorRampPalette(c('#004B7A', 'white', '#A83708'))(100)
myBreaks <- seq(-1.5, 1.5, length.out = 100)

# now create heatmaps based on the row order of the product
pdf(here::here("data/scenicplus/beta/results/high_quality_eregulon_exp_dev_heatmap_dotplot.pdf"),
    width = 2,
    height = 3)
# TF expression
mean_scale_df %>%
  dplyr::select(eregulon, condition, mean_exp) %>%
  tidyr::pivot_wider(id_cols = eregulon, names_from = condition, values_from = mean_exp) %>%
  dplyr::arrange(factor(eregulon, levels = unique(target))) %>%
  tibble::column_to_rownames("eregulon") %>%
  pheatmap::pheatmap(cluster_cols = FALSE, cluster_rows = FALSE, scale = "row", border_color = "grey", fontsize = 6,
                     color = myCol,
                     breaks = myBreaks)

# TF deviation
mean_scale_df %>%
  dplyr::select(eregulon, condition, mean_dev) %>%
  tidyr::pivot_wider(id_cols = eregulon, names_from = condition, values_from = mean_dev) %>%
  dplyr::arrange(factor(eregulon, levels = unique(target))) %>%
  tibble::column_to_rownames("eregulon") %>%
  pheatmap::pheatmap(cluster_cols = FALSE, cluster_rows = FALSE, scale = "row", border_color = "grey", fontsize = 6,
                     color = myCol,
                     breaks = myBreaks)

# dotplot gene AUC
mean_scale_df %>%
  dplyr::mutate(eregulon = factor(eregulon, levels = rev(target))) %>%
  ggplot(aes(x = condition, y = eregulon)) +
  geom_point(aes(size = scale_g_auc)) +
  ggplot2::scale_size(range = c(0, 5)) +
  labs(y = "High Quality Eregulons",
       x = "Conditions",
       size = "Gene AUC") +
  my_theme() +
  theme(axis.ticks = element_blank(),
        legend.title = element_text(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)) +
  ggplot2::guides(colour = "none", size = guide_legend(keyheight = 0.5, keywidth = 0.5))

# dotplot region AUC
mean_scale_df %>%
  dplyr::mutate(eregulon = factor(eregulon, levels = rev(target))) %>%
  ggplot(aes(x = condition, y = eregulon)) +
  geom_point(aes(size = scale_r_auc)) +
  ggplot2::scale_size(range = c(0, 5)) +
  labs(y = "High Quality Eregulons",
       x = "Conditions",
       size = "Region AUC") +
  my_theme() +
  theme(axis.ticks = element_blank(),
        legend.title = element_text(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)) +
  ggplot2::guides(colour = "none", size = guide_legend(keyheight = 0.5, keywidth = 0.5))

dev.off()
