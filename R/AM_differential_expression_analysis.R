# Description -------------------------------------------------------------
# Differential gene expression - likelihood ratio test and wald
# Setup -------------------------------------------------------------------
base::source(here::here("R/set_up.R"))
set.seed(1000)
create_directories(c(here::here("data/deseq2/rna_clusterbycond/"),
                     here::here("data/deseq2/rna_clusterbycond/figures/"),
                     here::here("data/deseq2/rna_clusterbycond/files/")))
# Load  -------------------------------------------------------------------
seurat_5 <- base::readRDS(here::here("data/seurat_objects/seurat_5.rds"))
augur <- base::readRDS(here::here("data/cell_prioritization/augur_across_all_conditions.rds"))

# Load variables ----------------------------------------------------------
padj_thres_lrt <- 0.05

# Pseudobulk counts -------------------------------------------------------
counts_list <- pseudobulk_list(
  object = seurat_5,
  sample = "orig.ident",
  cluster = "manual_anno",
  assay = "RNA",
  sample_levels = sample_levels,
  anno_levels = manual_anno_levels
)

# meta data file ----------------------------------------------------------
# Here we consider gene expression of mice fed LFD, HFD_1 or HFD 3. These conditions are
# classified using a group factor which contains three levels: LFD, HFD_1 and HFD_3.
meta <- base::data.frame(condition = base::gsub("_R[1-3]", "", sample_levels),
                   row.names = sample_levels) %>%
  dplyr::mutate(condition = base::factor(condition, levels = condition_levels))

# create dds object -------------------------------------------------------
# We will not include an intercept as we want to make multiple different comparisons
# and the results are basically the same when including the intercept,
# conceptually for me this is easier to understand
# https://github.com/tavareshugo/tutorial_DESeq2_contrasts/blob/main/DESeq2_contrasts.md
dds_list <- purrr::map2(
  counts_list, 
  base::rep(base::list(meta), base::length(manual_anno_levels)),
  ~ create_dds(counts = .x, meta = .y, design = ~ 0 + condition)
)

# PCA ---------------------------------------------------------------------
## log transform data
rld_list <- purrr::map(dds_list, get_rld)

## extract rld dataframe
rld_df <- purrr::map(rld_list, extract_rld)

### plot pca
gg_pca <- purrr::map2(rld_list, base::names(rld_list), ~ plot_pca(rld = .x, cluster = .y))

# Normalized counts -------------------------------------------------------
norm_counts <- purrr::map(dds_list, get_norm_counts)

# LRT TEST ----------------------------------------------------------------
dds_lrt <- purrr::map2(
  dds_list, 
  base::rep( base::list(~ 1),  base::length(dds_list)),
  ~ dds_test_lrt(dds = .x, test_use = "LRT", reduced = .y)
)

# Extract results from LRT ------------------------------------------------
res_lrt <- purrr::map(dds_lrt, get_results_lrt)

# Significantly expressed genes -------------------------------------------
res_lrt_sig <- pmap(
  base::list(res_lrt, base::rep(base::list(padj_thres_lrt), base::length(res_lrt)), base::rep("yes", base::length(res_lrt))),
  ~ get_sig(df = ..1, padj_thres = ..2, only_padj = ..3))

# Augur cell type prioritization plot -------------------------------------
## plot augur ----

## Number of significant genes ----
numb_genes <- purrr::map2(res_lrt_sig, base::names(res_lrt_sig), num_genes) %>% 
  bind_rows(.id = "cell") %>% 
  dplyr::arrange(-n)

# add steallte to make plot work
aug_df <- augur$AUC %>% 
  dplyr::add_row(cell_type = "Stellate", auc = 0)

# augur barplot ordered by AUC
gg_augur <- gg_augur_score(df = aug_df, cell_type_levels = numb_genes$cell)

## number of genes barplot
gg_num_bar <- gg_num_genes(df = numb_genes, cell_type_levels = numb_genes$cell)

# Cluster genes -----------------------------------------------------------
# Get significant genes
genes <- res_lrt_sig %>%
  purrr::modify_depth(1, ~ BiocGenerics::rownames(.))

# remove stellate cells, as they only have 1 diff gene. 
norm_counts <- norm_counts[-7]
genes <- genes[-7]

base::all.equal(base::names(norm_counts), base::names(genes))

# Determine parameters, remove stellate cells, as there is only one diff gene
fuzzy_pam <- purrr::pmap(base::list(genes, norm_counts, base::names(norm_counts)), 
                         ~ mfuzzy_params(genes = ..1, df = ..2, cluster = ..3))
# wrap plots
patchwork::wrap_plots(purrr::map(fuzzy_pam, dmin_plot))

# number of clusters
n_cl <- list("Beta" = 5,
             "Alpha" = 5,
             "Delta" = 5,
             "Gamma" = 5,
             "Acinar" = 5,
             "Endothalial" = 4,
             "Immune" = 4)
# cluster
cl <- purrr::map2(fuzzy_pam, n_cl, ~ mfuzz_cluster(exp = .x, n = .y))

# Plot heatmap ------------------------------------------------------------
myCol <- colorRampPalette(c('#004B7A', 'white', '#A83708'))(100)
myBreaks <- seq(-1.5, 1.5, length.out = 100)

# generate dataframe of membership for heat cluster
heat_cluster <- purrr::map(cl, heatmap_cluster)
base::all.equal(base::names(heat_cluster), base::names(norm_counts))

# create row annotation
row_anno <- purrr::map2(heat_cluster, norm_counts, heat_row_anno)

# annotation color
anno_color <- base::list(condition = condition_color,
                   membership = member_color)

# genes to cluster
genes_to_clust <- purrr::map(heat_cluster, ~ .x$gene)

# heatmap
gg_heat <- purrr::pmap(
  base::list(
    df = norm_counts,
    heat_color = rep(list(myCol), length(norm_counts)),
    breaks = rep(list(myBreaks), length(norm_counts)),
    title = names(norm_counts),
    gene = genes_to_clust,
    row_anno = row_anno,
    col_anno = rep(list(meta), length(norm_counts)),
    anno_color = rep(list(anno_color), length(norm_counts))
  ),
  gg_fuzzy_cluster_heatmap
)

# Pathway analysis using clusterprofiler ----------------------------------
## Create background dataset for hypergeometric testing using all
## genes tested for significance in the results

# GO-TERM ----

## Prepare data ----
# Generate dataframe of membership
gene_cluster_list <- purrr::map(cl, heatmap_cluster) %>%
  purrr::modify_depth(1, ~ column_to_rownames(.x, "gene"))

# Go term analysis
# remove stellate
res_lrt <- res_lrt[-7]
base::all.equal(base::names(gene_cluster_list), base::names(res_lrt))

go_term <- purrr::pmap(list(test_genes = gene_cluster_list,
                 test_cell = base::names(gene_cluster_list),
                 bg_genes = res_lrt,
                 bg_cell = base::names(res_lrt)), function(test_genes, test_cell, bg_genes, bg_cell) {
                   
                   base::print(base::paste0("testing:", test_cell))
                   
                   # Add entrez IDs to genes
                   test_genes$entrez <- AnnotationDbi::mapIds(org.Mm.eg.db::org.Mm.eg.db,
                                                              keys = BiocGenerics::rownames(test_genes),
                                                              column = "ENTREZID",
                                                              keytype = "SYMBOL",
                                                              multiVals = "first")
                   
                   # Add entrez IDs to background genes
                   bg_genes$entrez <- AnnotationDbi::mapIds(org.Mm.eg.db::org.Mm.eg.db,
                                                            keys = BiocGenerics::rownames(bg_genes),
                                                            column = "ENTREZID",
                                                            keytype = "SYMBOL",
                                                            multiVals = "first")
                   
                   base::print(base::paste0("background:", bg_cell))
                   
                   # Create universe
                   my_universe <- bg_genes %>%
                     BiocGenerics::as.data.frame() %>%
                     dplyr::filter(!base::is.na(entrez)) %>%
                     dplyr::pull(entrez) %>%
                     BiocGenerics::unlist() %>%
                     base::unique()
                   
                   # Split genes into clusters
                   gene_list <- test_genes %>%
                     base::split(f = base::as.factor(.$membership))
                   
                   # Perform pathway analysis
                   output <- purrr::map(gene_list, function(df) {
                     df %>%
                       dplyr::pull(entrez) %>%
                       base::unique() %>%
                       clusterProfiler::enrichGO(OrgDb = "org.Mm.eg.db",
                                                 pAdjustMethod = "fdr",
                                                 ont = "BP",
                                                 pvalueCutoff = 0.2,
                                                 qvalueCutoff = 0.2,
                                                 readable = TRUE,
                                                 universe = my_universe) %>%
                       dplyr::mutate(k = base::as.numeric(base::sub("/\\d+$", "", base::as.character(GeneRatio))),
                                     n = base::as.numeric(base::sub("^\\d+/", "", base::as.character(GeneRatio))),
                                     M = base::as.numeric(base::sub("/\\d+$", "", base::as.character(BgRatio))),
                                     N = base::as.numeric(base::sub("^\\d+/", "", base::as.character(BgRatio))),
                                     j_path = k / ((n + M) - k),
                                     j_bg = M / ((N + M) - M),
                                     is_sig = dplyr::case_when(p.adjust <= 0.05 ~ "*",
                                                               p.adjust > 0.05 ~ "ns"))
                   })
                   
                   return(output)
                 })


# Plot --------------------------------------------------------------------
# Get top 5 enriched pathways based on p-value
## GO TERM ----
go_term_top5 <- purrr::map(go_term, function(res) {
  
  purrr::map(res, function(res_2) {
    # Top 5 most significant pathways
    top_5_pathways <- res_2@result %>%
      dplyr::arrange(p.adjust) %>%
      utils::head(n = 5)  # Select top 5 most significant pathways
    
    return(top_5_pathways)
  })
})


# Save gene list for each cluster -----------------------------------------
genesets_deseq <- gene_cluster_list %>%
  purrr::modify_depth(1, ~ base::split(., f = as.factor(.$membership))) %>%
  purrr::modify_depth(2, ~ BiocGenerics::rownames(.))

base::saveRDS(genesets_deseq, file = here::here("data/deseq2/rna_clusterbycond/files/genesets_deseq_analysis.rds"))

# Save data ---------------------------------------------------------------
## meta ----
base::saveRDS(meta, file = here::here("data/deseq2/rna_clusterbycond/files/meta.rds"))
openxlsx::write.xlsx(meta %>% tibble::rownames_to_column("replicate"), here::here("data/deseq2/rna_clusterbycond/files/meta.xlsx"))

## pseudobulk counts ----
base::saveRDS(counts_list, file = here::here("data/deseq2/rna_clusterbycond/files/pseudobulk_counts.rds"))

## rlog counts ----
base::saveRDS(rld_df, file = here::here("data/deseq2/rna_clusterbycond/files/rlog_counts.rds"))

## normalized counts ----
base::saveRDS(norm_counts, file = here::here("data/deseq2/rna_clusterbycond/files/norm_pseudobulk_counts_per_cluster.rds"))

## Results LRT ----
base::saveRDS(res_lrt, file = here::here("data/deseq2/rna_clusterbycond/files/res_lrt.rds"))
openxlsx::write.xlsx(res_lrt %>% 
                       purrr::modify_depth(1, ~BiocGenerics::as.data.frame(.) %>% 
                                             tibble::rownames_to_column("gene")), 
                     here::here("data/deseq2/rna_clusterbycond/files/res_lrt.xlsx"))

## Significant results LRT ----
base::saveRDS(res_lrt_sig, file = here::here("data/deseq2/rna_clusterbycond/files/res_lrt_significant.rds"))
openxlsx::write.xlsx(res_lrt_sig %>% 
                       purrr::modify_depth(1, ~BiocGenerics::as.data.frame(.) %>% 
                                             tibble::rownames_to_column("gene")), here::here("data/deseq2/rna_clusterbycond/files/res_lrt_significant.xlsx"))

## Fuzzy cluster parameters ----
base::saveRDS(fuzzy_pam, file = here::here("data/deseq2/rna_clusterbycond/files/fuzzy_parameters.rds"))

## Identified clusters ----
base::saveRDS(cl, file = here::here("data/deseq2/rna_clusterbycond/files/clustering_new.rds"))

## GO-term results ----
base::saveRDS(go_term, file = here::here("data/deseq2/rna_clusterbycond/files/clusterprofiler_go_biolgoical_process.rds"))

# save plots
## PCA ----
pdf(here::here("data/deseq2/rna_clusterbycond/figures/pca.pdf"),
    width = 2,
    height = 2)
gg_pca
dev.off()

## Heatmap LRT ----
pdf(here::here("data/deseq2/rna_clusterbycond/figures/gene_clusters_heatmap.pdf"),
    width = 9.18,
    height = 9.18)
gg_heat
dev.off()

## AUGUR ----
pdf(here::here("data/deseq2/rna_clusterbycond/figures/augur_and_number_of_genes.pdf"),
    width = 6,
    height = 2)
gg_augur | gg_num_bar | gg_pca[["Beta"]]
dev.off()

## Fuzzy parameters ----
pdf(here::here("data/deseq2/rna_clusterbycond/figures/minimum_centroid_distance_fuzzy.pdf"),
    width = 2,
    height = 2)
Map(dmin_plot,
    fuzzy_pam)
dev.off()

# Plot gene examples as bar plots -----------------------------------------
# for beta-cells
norm_beta <- norm_counts[["Beta"]]

gene_examp <- list("gene_cluster_1" = c("Slc2a2", "Nkx6-1", "Mafa", "Rfx6", "Glp1r"),
                   "gene_cluster_2" = c("Stat1", "Stat2", "Irf7", "Ifit1", "Cxcl10"),
                   "gene_cluster_3" = c("Ndufs1", "Ndufs3", "Smim20", "Pet100", "Tmem223"),
                   "gene_cluster 4" = c("Hmgcr", "Sqle", "Dhcr7", "Cyp51", "Fdft1"),
                   "gene_cluster_5" = c("Clock", "Per1", "Per2", "Cry1", "Rorc"))


pdf(here::here("data/deseq2/rna_clusterbycond/figures/gene_cluster_examples.pdf"),
    height = 6,
    width = 3)
gene_examp %>%
  tibble::enframe() %>%
  tidyr::unnest(value) %>%
  dplyr::rename(gene = value,
                cluster = name) %>%
  dplyr::right_join(x = tibble::rownames_to_column(BiocGenerics::as.data.frame(norm_beta), "gene"), by = "gene") %>%
  tidyr::pivot_longer(c(-gene, - cluster), names_to = "sample", values_to = "exp") %>%
  dplyr::mutate(condition = stringr::str_extract(sample, "LFD|HFD_1|HFD_3"),
                condition = factor(condition, levels = condition_levels),
                cluster = factor(cluster, levels = names(gene_examp))) %>%
  ggplot2::ggplot(aes(x = condition, y = exp, fill = condition)) +
  ggplot2::geom_bar(stat="summary", fun=mean, aes(fill=factor(condition))) +
  ggplot2::geom_point(size = 0.5, color = "black")+
  ggplot2::facet_wrap(~ cluster + gene, scales = "free", ncol = 5) +
  ggplot2::scale_fill_manual(values = condition_color) +
  ggplot2::labs(x = "condition",
                y = "Pseudo Bulk Normalized Counts") +
  my_theme() +
  ggplot2::theme(legend.position = "none", strip.background = element_blank(),
                 axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

# Total cell count vs number of diff genes --------------------------------
n_cells <- seurat_5@meta.data %>%
  as.data.frame() %>%
  dplyr::group_by(manual_anno) %>%
  dplyr::summarize(n_cluster = n())

n_genes <- res_lrt_sig %>%
  purrr::modify_depth(1, ~ dim(.) %>%
                        .[1]) %>%
  unlist() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("manual_anno") %>%
  dplyr::rename("n_genes" = '.')

# augur vs number of cells ------------------------------------------------
augur_auc <- augur$AUC %>%
  dplyr::rename(manual_anno = cell_type)

df <- merge(n_cells, augur_auc)
stats::cor.test(df$n_cluster, df$auc)

df_2 <- merge(n_cells, n_genes)
stats::cor.test(df_2$n_cluster, df_2$n_genes)

df_3 <- merge(n_genes, augur_auc)
stats::cor.test(df_3$n_genes, df_3$auc)

pdf(here::here("data/deseq2/rna_clusterbycond/figures/n_cells_vs_n_genes.pdf"),
    height = 2,
    width = 2)
merge(n_cells, n_genes) %>%
  ggplot2::ggplot(aes(x = n_cluster, y = n_genes)) +
  ggplot2::geom_point(aes(color = manual_anno), size = 1) +
  ggpubr::stat_cor() +
  ggplot2::labs(x = "#Cells",
                y = "#DEGs") +
  ggplot2::scale_color_manual(values = cluster_anno) +
  my_theme() +
  ggplot2::theme(legend.position = "none")

merge(n_cells, augur_auc) %>%
  ggplot2::ggplot(aes(x = n_cluster, y = auc)) +
  ggplot2::geom_point(aes(color = manual_anno), size = 1) +
  ggpubr::stat_cor() +
  ggplot2::labs(x = "#Cells",
                y = "AUC") +
  ggplot2::scale_color_manual(values = cluster_anno) +
  ggplot2::ylim(0.5, 1) +
  my_theme() +
  ggplot2::theme(legend.position = "none")

merge(n_genes, augur_auc) %>%
  ggplot2::ggplot(aes(x = n_genes, y = auc)) +
  ggplot2::geom_point(aes(color = manual_anno), size = 1) +
  ggpubr::stat_cor() +
  ggplot2::labs(x = "#DEGs",
                y = "AUC") +
  ggplot2::scale_color_manual(values = cluster_anno) +
  ggplot2::ylim(0.5, 1) +
  my_theme() +
  ggplot2::theme(legend.position = "none")
dev.off()
