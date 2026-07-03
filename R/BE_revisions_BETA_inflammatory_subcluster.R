# Description -------------------------------------------------------------
# Adress reviewer question about:
# "Examine heterogeneity in inflammatory response in greater detail using existing data (Reviewer 2 point 7)." 

# Setup -------------------------------------------------------------------
base::source(here::here("R/set_up_revisions.R"))
set.seed(1000)
create_directories(c(here::here("data/revisions/")))

# Load --------------------------------------------------------------------
seurat_beta <-  readRDS(here::here("data/inflammatory_cluster/files/seurat_beta.rds"))

# Visualize on umap -------------------------------------------------------
pdf(here::here("data/revisions/high_inflammatory_cells_subcluster.pdf"), height = 2, width = 3)
p1 <- DimPlot(seurat_beta, group.by = "is_inflammatory_ucell", reduction = "umap.wnn", pt.size = 0.1) +
  ggplot2::scale_color_manual(values = c("low" = "grey", "high" = "#FA8231")) + my_theme_void()
p2 <- DimPlot(seurat_beta, group.by = "seurat_clusters", reduction = "umap.wnn", pt.size = 0.1) + my_theme_void()

p1 + p2
dev.off()

pdf(here::here("data/revisions/stat1_stat2_expression.pdf"), height = 5, width = 6)
FeaturePlot(seurat_beta, features = c("Stat1", "Stat2"), reduction = "umap.wnn", keep.scale = "all", split.by = "condition", order = TRUE, pt.size = 0.1) & my_theme_void() &
  ggplot2::theme(legend.position = "right") &
  patchwork::plot_layout(guides = "collect")
dev.off()

# Subcluster --------------------------------------------------------------
seurat_beta <- seurat_beta  |> 
  Seurat::FindMultiModalNeighbors(
    reduction.list = list("jointly", "liger_embeddings"),
    dims.list = list(1:15, 1:20))

seurat_beta <- Seurat::FindClusters(
  seurat_beta,
  graph.name = "wsnn",
  res = 0.2,
  random.seed = 1000,
  verbose = TRUE)

cluster_composition <- seurat_beta@meta.data |> 
  dplyr::group_by(seurat_clusters, is_inflammatory_ucell) |> 
  dplyr::summarise(n = n(), .groups = "drop") |> 
  dplyr::group_by(seurat_clusters) |> 
  dplyr::mutate(pct = n / sum(n) * 100) |> 
  dplyr::filter(is_inflammatory_ucell == "high")


# DEGs cluster 3 ----------------------------------------------------------
seurat_beta$cluster_status <- dplyr::if_else(
  seurat_beta$seurat_clusters == "3", 
  "cluster_3", 
  "other"
)

counts_cluster3 <- seurat_beta %>%
  edgeR::Seurat2PB(sample = "orig.ident",
                   cluster = "cluster_status")

counts_cluster3_df <- counts_cluster3 %>%
  as.data.frame() %>%
  dplyr::rename_with(~gsub("cluster_status", "", .x)) %>%
  tibble::column_to_rownames("gene")

meta_data_cluster3 <- counts_cluster3[["samples"]] %>%
  magrittr::set_rownames(base::gsub("cluster_status", "", BiocGenerics::rownames(.))) %>%
  dplyr::mutate(
    # Remove everything from "_cluster" onwards
    sample = stringr::str_remove(BiocGenerics::rownames(.), "_cluster.*"),
    cluster = cluster
  )

all.equal(BiocGenerics::colnames(counts_cluster3_df), BiocGenerics::rownames(meta_data_cluster3))

dds_cluster3 <- DESeq2::DESeqDataSetFromMatrix(
  countData = counts_cluster3_df,
  colData = meta_data_cluster3,
  design = ~ sample + cluster
)

dds_cluster3 <- DESeq2::DESeq(dds_cluster3)

res_cluster3 <- DESeq2::results(
  dds_cluster3, 
  contrast = c("cluster", "cluster_3", "other")
)

# Percentage of genes expressed
expr_mat <- seurat_beta@assays$RNA@data

# Identify cluster 3 and other cells
cluster3_cells <- WhichCells(seurat_beta, expression = seurat_clusters == "3")
other_cells <- WhichCells(seurat_beta, expression = seurat_clusters != "3")

# Calculate percentage expressing (> 0) for each gene
pct_cluster3 <- rowSums(expr_mat[, cluster3_cells] > 0) / length(cluster3_cells) * 100
pct_other <- rowSums(expr_mat[, other_cells] > 0) / length(other_cells) * 100

# Create a dataframe with percentages
pct_df <- tibble::tibble(
  gene = names(pct_cluster3),
  pct_cluster3 = pct_cluster3,
  pct_other = pct_other,
  pct_diff = pct_cluster3 - pct_other
)

res_cluster3_df <- res_cluster3_df %>%
  dplyr::left_join(pct_df, by = "gene")

markers <- res_cluster3_df %>%
  dplyr::filter(pct_cluster3 > 50) %>% 
  dplyr::filter(log2FoldChange > 0, padj < 0.05) %>%
  dplyr::arrange(padj, dplyr::desc(log2FoldChange), dplyr::desc(pct_cluster3)) %>% 
  dplyr::mutate("tnf_ifn_target" = dplyr::case_when(gene %in% inf_markers ~ TRUE,
                                                    .default = FALSE)) %>% 
  tibble::rowid_to_column("rank") %>% 
  dplyr::relocate(rank, .after = "tnf_ifn_target")

# View with percentages
res_top_10 <- res_cluster3_df %>%
  dplyr::filter(pct_cluster3 > 50) %>% 
  dplyr::filter(log2FoldChange > 0, padj < 0.05) %>%
  dplyr::arrange(padj, dplyr::desc(log2FoldChange), dplyr::desc(pct_cluster3)) %>% 
  head(n = 10) %>%
  dplyr::select(gene, log2FoldChange, padj, pct_cluster3, pct_other, pct_diff)

p1 <- Seurat::DotPlot(seurat_beta, group.by = "seurat_clusters",
                      features = res_top_10$gene,
                      assay = "RNA",
                      cluster.idents = FALSE) +
  ggplot2::scale_colour_gradient2(low = "#004B7AFF", mid = "#FDFDFCFF", high = "#A83708FF") +
  my_theme() +
  ggplot2::coord_flip() +
  ggplot2::theme(panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(), 
                 axis.title = element_blank(),
                 legend.position = "none")
p2 <- Seurat::DotPlot(seurat_beta, group.by = "seurat_clusters",
                      features = res_top_10$gene,
                      assay = "activity",
                      cluster.idents = FALSE) +
  ggplot2::scale_colour_gradient2(low = "#004B7AFF", mid = "#FDFDFCFF", high = "#A83708FF") +
  my_theme() +
  ggplot2::coord_flip() +
  ggplot2::theme(panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(), 
                 axis.title = element_blank(), legend.position = "none",
                 axis.text.y = element_blank(),
                 axis.ticks.y = element_blank())

pdf(here::here("data/revisions/high_inflammatory_cells_subcluster_top10_dotplot.pdf"), height = 2, width = 2.5)
p1 + p2
dev.off()


# Save --------------------------------------------------------------------
vroom::vroom_write(res_top_10, here::here("data/revisions/diff_genes_cluster3_top10.csv"))
vroom::vroom_write(res_cluster3_df, here::here("data/revisions/diff_genes_cluster3.csv"))
vroom::vroom_write(markers, here::here("data/revisions/cluster3_markers_all.csv"))
saveRDS(seurat_beta, here::here("data/revisions/seurat_beta_subcluster.rds"))
vroom::vroom_write(cluster_composition, here::here("data/revisions/cluster_composisiton.csv"))



