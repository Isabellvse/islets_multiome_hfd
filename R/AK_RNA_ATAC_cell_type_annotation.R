# Description -------------------------------------------------------------
# Annotate clusters manually

# Setup -------------------------------------------------------------------
base::source(here::here("R/set_up.R"))
set.seed(1000)

# Setting a Genome and GeneAnnotation
ArchR::addArchRGenome("mm10")

create_directories(c(here::here("data/dimensional_reduction/wnn/celltype_annotation/manual/"),
                     here::here("data/deseq2/celltype/files/")))

# Load --------------------------------------------------------------------
seurat_5 <- base::readRDS(here::here("data/seurat_objects/seurat_5.rds"))
multi_islets_4 <- ArchR::loadArchRProject(path = here::here("data/archr_projects/save_multi_islets_4/"))

# Manual annotation -------------------------------------------------------
Seurat::DefaultAssay(seurat_5) <- "RNA"
Seurat::Idents(seurat_5) <- "wsnn_res.0.8"

### Dotplot of initial clustering ----
# gene expression
p1 <- Seurat::DotPlot(seurat_5, group.by = "wsnn_res.0.8",
        features = markers_short,
        assay = "RNA",
        cluster.idents = FALSE) +
  ggplot2::scale_colour_gradient2(low = "#004B7AFF", mid = "#FDFDFCFF", high = "#A83708FF") +
  my_theme() +
  ggplot2::theme(panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(), 
                 axis.title.x = element_blank())
ggsave(filename = here::here("data/dimensional_reduction/wnn/celltype_annotation/manual/azimuth_marker_gene_dotplot_wsnn_res.0.8_gene.pdf"),
       plot = p1,
       height = 6/1,
       width = 6,
       dpi = 1000)

# gene activity
p2 <- Seurat::DotPlot(seurat_5, group.by = "wsnn_res.0.8",
                features = markers_short,
                assay = "activity",
                cluster.idents = FALSE) +
  ggplot2::scale_colour_gradient2(low = "#004B7AFF", mid = "#FDFDFCFF", high = "#A83708FF") +
  my_theme() +
  ggplot2::theme(panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(), 
                 axis.title.x = element_blank())

ggsave(filename = here::here("data/dimensional_reduction/wnn/celltype_annotation/manual/azimuth_marker_gene_dotplot_wsnn_res.0.8_activity.pdf"),
       plot = p2,
       height = 6/1,
       width = 6,
       dpi = 1000)

### Add manual annotation -----
seurat_5@meta.data <- seurat_5@meta.data %>%
  dplyr::mutate(manual_anno = dplyr::case_when(wsnn_res.0.8 == "0" ~ "Beta",
                                               wsnn_res.0.8 == "1" ~ "Beta",
                                               wsnn_res.0.8 == "2" ~ "Beta",
                                               wsnn_res.0.8 == "3" ~ "Beta",
                                               wsnn_res.0.8 == "4" ~ "Beta",
                                               wsnn_res.0.8 == "5" ~ "Beta",
                                               wsnn_res.0.8 == "6" ~ "Alpha",
                                               wsnn_res.0.8 == "7" ~ "Beta",
                                               wsnn_res.0.8 == "8" ~ "Beta",
                                               wsnn_res.0.8 == "9" ~ "Beta",
                                               wsnn_res.0.8 == "10" ~ "Beta",
                                               wsnn_res.0.8 == "11" ~ "Delta",
                                               wsnn_res.0.8 == "12" ~ "Alpha",
                                               wsnn_res.0.8 == "13" ~ "Gamma",
                                               wsnn_res.0.8 == "14" ~ "Endothelial",
                                               wsnn_res.0.8 == "15" ~ "Beta",
                                               wsnn_res.0.8 == "16" ~ "Endothelial",
                                               wsnn_res.0.8 == "17" ~ "Endothelial",
                                               wsnn_res.0.8 == "18" ~ "Acinar",
                                               wsnn_res.0.8 == "19" ~ "Immune",
                                               wsnn_res.0.8 == "20" ~ "Immune",
                                               wsnn_res.0.8 == "21" ~ "Stellate",
                                               wsnn_res.0.8 == "22" ~ "Beta"),
                manual_anno = factor(manual_anno, levels = c("Beta", "Alpha", "Delta", "Gamma", "Acinar",
                                                             "Endothelial", "Stellate", "Immune")))

### Dotplot of manual annotation ----
# gene expression
p3 <- DotPlot(seurat_5, group.by = "manual_anno",
        features = markers_short,
        assay = "RNA",
        cluster.idents = FALSE, scale.min = 0, col.max = 100,
        dot.scale = 1.5) +
  ggplot2::scale_colour_gradient2(low = "#004B7AFF", mid = "#FDFDFCFF", high = "#A83708FF") +
  my_theme() +
  ggplot2::theme(panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(), 
                 axis.title.x = element_blank(), 
                 strip.text = ggplot2::element_blank(),
                 axis.text.x = element_blank(),
                 axis.title.y = element_blank()) +
  ggplot2::guides(colour = "none", size = guide_legend(keyheight = 0.5, keywidth = 0.5))
# gene activity
p4 <- DotPlot(seurat_5, group.by = "manual_anno",
        features = markers_short,
        assay = "activity",
        cluster.idents = FALSE,
        scale.min = 0, col.max = 100,
        dot.scale = 1.5) +
  ggplot2::scale_colour_gradient2(low = "#004B7AFF", mid = "#FDFDFCFF", high = "#A83708FF") +
  my_theme() +
  ggplot2::theme(panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(), 
                 axis.title.x = element_blank(),
                 axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                 strip.text = ggplot2::element_blank(),
                 axis.title.y = element_blank()) +
  ggplot2::guides(colour = ggplot2::guide_colourbar(barwidth = 0.5, barheight = 6/2),
                  size = "none")

ggsave(filename = here::here("data/dimensional_reduction/wnn/celltype_annotation/manual/manual_annotation_dotplot_marker_genes.pdf"),
       plot = p3/p4,
       height = 6/3,
       width = 3.5,
       dpi = 1000)

### Dimplot - manual clusters ----
p5 <- seurat_5 %>%
  Seurat::DimPlot(
    reduction = "umap.wnn",
    group.by = "manual_anno",
    cols = cluster_anno, 
    pt.size = 3, 
    raster = TRUE, 
    raster.dpi = c(800, 800)) +
  my_theme_void(remove_strip_text = FALSE) +
  ggplot2::theme(title = element_blank(),
                 legend.position = "none")

p6 <- seurat_5 %>%
  Seurat::DimPlot(
    reduction = "umap.atac",
    group.by = "manual_anno",
    cols = cluster_anno, 
    pt.size = 3, 
    raster = TRUE, 
    raster.dpi = c(800, 800)) +
  my_theme_void(remove_strip_text = FALSE) +
  ggplot2::theme(title = element_blank(),
                 legend.position = "none")


p7 <-seurat_5 %>%
  Seurat::DimPlot(
    reduction = "umap",
    group.by = "manual_anno",
    cols = cluster_anno, 
    pt.size = 3, 
    raster = TRUE, 
    raster.dpi = c(800, 800)) +
  my_theme_void(remove_strip_text = FALSE) +
  ggplot2::theme(title = element_blank(),
                 legend.position = "none")

ggsave(filename = here::here("data/dimensional_reduction/wnn/celltype_annotation/manual/manual_annotation_dimplot_wnn_atac_rna.pdf"),
       plot = p5+p6+p7,
       height = 6/3,
       width =6,
       dpi = 1000)

### RNA and ATAC weights ----
pdf(file = here::here("data/dimensional_reduction/wnn/celltype_annotation/manual/manual_annotation_vlnplot_weights.pdf"),
    height = 3.5,
    width = 1.5)

seurat_5 %>%
  Seurat::VlnPlot(group.by = "manual_anno",
                  features = c("RNA.weight", "peaks.weight"),
                  cols = cluster_anno,
                  pt.size = 0, ncol = 1) &
  ggplot2::geom_boxplot(width = 0.3, alpha = 0, lwd = 0.5) &
  my_theme() +
  ggplot2::theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
dev.off()

### Number of cells in each cluster ----
df <- seurat_5@meta.data %>%
  as.data.frame() %>%
  dplyr::group_by(manual_anno) %>%
  dplyr::summarise(n = n())

utils::write.csv(df, file = here::here("data/dimensional_reduction/wnn/celltype_annotation/manual/manual_annotation_numberofcellspercluster.csv"))

### Perc of endocrine cells ----
df <- seurat_5@meta.data %>%
  as.data.frame() %>%
  dplyr::filter(manual_anno == "Beta"|
                  manual_anno == "Alpha" |
                  manual_anno == "Delta" |
                  manual_anno == "Gamma") %>%
  dplyr::group_by(manual_anno) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(perc = round((n/sum(n))*100, 2))

utils::write.csv(df, file = here::here("data/dimensional_reduction/wnn/celltype_annotation/manual/manual_annotation_percofendocrinecells.csv"))

### Number of cells in each cluster - per condition ----
df <- seurat_5@meta.data %>%
  as.data.frame() %>%
  dplyr::select(manual_anno, condition) %>%
  dplyr::group_by(condition, manual_anno) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(sum = sum(n),
                perc = round((n/sum)*100, 2))

utils::write.csv(df, file = here::here("data/dimensional_reduction/wnn/celltype_annotation/manual/manual_annotation_numberofcellspercluster_condition.csv"))

### Number of cells in each cluster - per condition bar plot ----
p8 <- df %>%
  ggplot2::ggplot(aes(x = condition, y = perc, fill = manual_anno)) +
  ggplot2::geom_bar(stat = "identity") +
  ggplot2::scale_fill_manual(values = cluster_anno) +
  ggplot2::labs(y = "% of cells per condition") +
  my_theme() +
  ggplot2::theme(legend.position = "none",
                 axis.title.x = element_blank())
ggplot2::ggsave(plot = p8,
                file = here::here("data/dimensional_reduction/wnn/celltype_annotation/manual/manual_annotation_perc_all_cell.pdf"),
                height = 2,
                width = 1,
                dpi = 1000)

### Perc of endocrine cells - per condition ----
df <- seurat_5@meta.data %>%
  as.data.frame() %>%
  dplyr::filter(manual_anno == "Beta"|
                  manual_anno == "Alpha" |
                  manual_anno == "Delta" |
                  manual_anno == "Gamma") %>%
  dplyr::group_by(condition, manual_anno) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(perc = round((n/sum(n))*100, 2))

utils::write.csv(df, file = here::here("data/dimensional_reduction/wnn/celltype_annotation/manual/manual_annotation_percofendocrinecells_condition.csv"))

### Number of cells per replicate
### Number of cells in each cluster - per condition ----
df <- seurat_5@meta.data %>%
  as.data.frame() %>%
  dplyr::select(manual_anno, orig.ident) %>%
  dplyr::group_by(orig.ident, manual_anno) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(sum = sum(n),
                perc = round((n/sum)*100, 3))

utils::write.csv(df, file = here::here("data/dimensional_reduction/wnn/celltype_annotation/manual/manual_annotation_percofeallcells_sample.csv"))

p9 <- df %>%
  ggplot2::ggplot(aes(x = orig.ident, y = perc, fill = manual_anno)) +
  ggplot2::geom_bar(stat = "identity") +
  ggplot2::scale_fill_manual(values = cluster_anno) +
  ggplot2::labs(y = "% of total number of cells") +
  my_theme() +
  ggplot2::theme(legend.position = "none",
                 axis.title.x = element_blank())
ggplot2::ggsave(plot = p9,
                file = here::here("data/dimensional_reduction/wnn/celltype_annotation/manual/manual_annotation_perc_all_cell_sample.pdf"),
                height = 2,
                width = 2,
                dpi = 1000)

### Number of endocrine cells per replicate
df <- seurat_5@meta.data %>%
  as.data.frame() %>%
  dplyr::filter(manual_anno == "Beta"|
                  manual_anno == "Alpha" |
                  manual_anno == "Delta" |
                  manual_anno == "Gamma") %>%
  dplyr::select(manual_anno, orig.ident) %>%
  dplyr::group_by(orig.ident, manual_anno) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(sum = sum(n),
                perc = round((n/sum)*100, 3))

utils::write.csv(df, file = here::here("data/dimensional_reduction/wnn/celltype_annotation/manual/manual_annotation_percofeendocrine_sample.csv"))


### Coverage plot marker genes ----
# set identity
Idents(seurat_5) <- "manual_anno"
pdf(file = here::here("data/dimensional_reduction/wnn/celltype_annotation/manual/manual_annotation_coverage_markers.pdf"),
    height = 9.18,
    width = 5)
Map(coverage_plot,
    obj = rep(list(seurat_5), length(markers_short_tn)),
    region = markers_short_tn,
    fill = rep(list(cluster_anno), length(markers_short_tn)))
dev.off()

cov_plot <- Map(coverage_plot,
                obj = rep(list(seurat_5), length(markers_short_tn)),
                region = markers_short_tn,
                fill = rep(list(cluster_anno), length(markers_short_tn)))

cov_plot_2 <- wrap_plots(cov_plot, nrow = 1)
ggsave(filename = here::here("data/dimensional_reduction/wnn/celltype_annotation/manual/manual_annotation_coverage_markers.pdf"),
       plot = cov_plot_2,
       height = 6/2,
       width =3,
       dpi = 1000)
# Save seurat object ------------------------------------------------------
base::saveRDS(seurat_5, file = here::here("data/seurat_objects/seurat_5.rds"))

# Differential gene expression per cell type ------------------------------
### Prepare data ----

# Set default assay
Seurat::DefaultAssay(seurat_5) <- "RNA"

# Create columns to use for comparisons
seurat_5@meta.data <- seurat_5@meta.data %>%
  BiocGenerics::as.data.frame() %>%
  dplyr::mutate(Beta_vs_other = case_when(manual_anno == "Beta" ~ "beta",
                                          manual_anno != "Beta" ~ "other"),
                Alpha_vs_other = case_when(manual_anno == "Alpha" ~ "alpha",
                                           manual_anno != "Alpha" ~ "other"),
                Delta_vs_other = case_when(manual_anno == "Delta" ~ "delta",
                                           manual_anno != "Delta" ~ "other"),
                Gamma_vs_other = case_when(manual_anno == "Gamma" ~ "gamma",
                                           manual_anno != "Gamma" ~ "other"),
                Acinar_vs_other = case_when(manual_anno == "Acinar" ~ "acinar",
                                            manual_anno != "Acinar" ~ "other"),
                Stellate_vs_other = case_when(manual_anno == "Stellate" ~ "stellate",
                                              manual_anno != "Stellate" ~ "other"),
                Endothelial_vs_other = case_when(manual_anno == "Endothelial" ~ "endothelial",
                                                 manual_anno != "Endothelial" ~ "other"),
                Immune_vs_other = case_when(manual_anno == "Immune" ~ "immune",
                                            manual_anno != "Immune" ~ "other"))



# Find differentially expressed genes -------------------------------------
counts_list <- list()
meta_list <- list()

for (i in base::names(cluster_anno)) {
  
  # Pseudobulked counts
  counts <- seurat_5 %>%
    edgeR::Seurat2PB(sample = "orig.ident",
                     cluster = base::paste0(i, "_vs_other"))
  
  # Get counts
  counts_2 <- counts[["counts"]] %>%
    BiocGenerics::as.data.frame() %>%
    dplyr::rename_with(~base::gsub("cluster", "", .x))
  
  # Get meta data
  meta_data <- counts[["samples"]] %>%
    magrittr::set_rownames(base::gsub("cluster", "", BiocGenerics::rownames(.))) %>%
    dplyr::mutate(condition = stringr::str_extract(rownames(.), "LFD|HFD_1|HFD_3"),
                  condition_inflammatory = base::factor(base::paste0(condition, "_", cluster)))
  
  # check that colnames in counts_2 and rownames in meta_data are the same order
  base::print(base::all.equal(BiocGenerics::colnames(counts_2), BiocGenerics::rownames(meta_data)))
  
  counts_list[[i]] <- counts_2
  meta_list[[i]] <- meta_data
}

# Define the comparisons to be performed
comparisons <- base::list(
  "Beta" = list("design" = "~ sample + cluster", "contrast" = c("cluster", "beta", "other")),
  "Alpha" = list("design" = "~ sample + cluster", "contrast" = c("cluster", "alpha", "other")),
  "Delta" = list("design" = "~ sample + cluster", "contrast" = c("cluster", "delta", "other")),
  "Gamma" = list("design" = "~ sample + cluster", "contrast" = c("cluster", "gamma", "other")),
  "Acinar" = list("design" = "~ sample + cluster", "contrast" = c("cluster", "acinar", "other")),
  "Stellate" = list("design" = "~ sample + cluster", "contrast" = c("cluster", "stellate", "other")),
  "Endothelial" = list("design" = "~ sample + cluster", "contrast" = c("cluster", "endothelial", "other")),
  "Immune" = list("design" = "~ sample + cluster", "contrast" = c("cluster", "immune", "other"))
)

# Loop through the comparisons
results <- list()
dds_list <- list()
for (i in base::names(counts_list)) {
  gene_counts <- counts_list[[i]]
  meta_data <- meta_list[[i]]
  comparison <- comparisons[[i]]  # Use 'i' to index the comparisons list
  
  # Create DESeq2 dataset with the specified design
  dds <- DESeqDataSetFromMatrix(
    countData = gene_counts,
    colData = meta_data,
    design = stats::as.formula(comparison$design)
  )
  
  # Normalize and run DESeq
  dds <- DESeq(dds)
  
  # Run Wald test for the specific contrast
  res <- results(dds, contrast = comparison$contrast)
  
  # Store the result in the list
  results[[i]] <- res  # Use 'i' instead of 'comp_name'
  dds_list[[i]] <- dds  # Use 'i' instead of 'comp_name'
}

# Save --------------------------------------------------------------------
base::saveRDS(results, here::here("data/deseq2/celltype/files/marker_genes_results.rds"))

# Save as table -----------------------------------------------------------
# Filter significant results
sig_res <- purrr::map2(results, names(results), function(df, name) {
  # Convert df to tibble, add rownames as "gene" column, and filter significant results
  df_2 <- df %>%
    BiocGenerics::as.data.frame() %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::arrange(padj, desc(log2FoldChange)) %>% 
    dplyr::filter(padj <= 0.05)
  
  # Write the significant results to CSV
  utils::write.csv(df_2, file = base::paste0(here::here("data/deseq2/celltype/files/"), name, "_sig_results.csv"), row.names = FALSE)
  return(df_2)
  })

# save as excel -----------------------------------------------------------
openxlsx::write.xlsx(sig_res, here::here("data/deseq2/celltype/files/celltype_vs_other_sig_results.xlsx"))


# save top 40 genes -------------------------------------------------------

sig_res_top <- purrr::map(sig_res, function(df){
  df %>% dplyr::filter(log2FoldChange > 0) %>% 
    dplyr::arrange(padj, desc(log2FoldChange)) %>% 
    head(n = 40)
})

openxlsx::write.xlsx(sig_res_top, here::here("data/deseq2/celltype/files/celltype_vs_other_sig_results_top_40.xlsx"))

# Quality control by diet -------------------------------------------------
# also run UMAP based on PCA and iLISI
Seurat::DefaultAssay(seurat_5) <- "RNA"
seurat_test <- seurat_5
seurat_test <- Seurat::FindVariableFeatures(seurat_test, selection.method = "vst", nfeatures = 1000)
seurat_test <- Seurat::ScaleData(seurat_test)
seurat_test <- Seurat::RunPCA(seurat_test, features = Seurat::VariableFeatures(object = seurat_test))

# Get lisi embeddings

svd <- multi_islets_4@reducedDims$IterativeLSI_peaks@listData[["matSVD"]]
BiocGenerics::rownames(svd) <- atac_to_rna_syntax(BiocGenerics::rownames(svd))

## split by diet ----
metadata <- seurat_5@meta.data %>%
  base::split(f = factor(.$condition))

embeddings <- list("pca" = seurat_test@reductions$pca@cell.embeddings[,1:15],
                   "jointly" = seurat_5@reductions$jointly@cell.embeddings,
                   "svd" = svd,
                   "liger" = seurat_5@reductions$liger_embeddings@cell.embeddings[,1:20])
q <- list()
for (i in 1:length(metadata)){
  for (j in 1:length(embeddings)){
    # Define variables
    batch_var <- "orig.ident"
    label_var <- "manual_anno"
    metadata_keep <- BiocGenerics::rownames(metadata[[i]])
    
    embed <- embeddings[[j]][metadata_keep,]
    embed_name <- base::names(embeddings)[j]
    
    dims <- base::length(BiocGenerics::colnames(embed))
    
    # run test
    test <- evaluateEmbedding(x = embed,
                              metadata = metadata[[i]],
                              dims = dims,
                              batch_var = batch_var,
                              label_var = label_var)
    
    # create dataframe of results
    df <- base::data.frame("condition" = base::names(metadata)[i],
                           "method" = embed_name,
                           test[["summary"]]["globaliLISI"],
                           test[["summary"]]["globalbASW"])
    
    q[[base::paste0(embed_name, base::names(metadata)[i])]] <- df
  }
}

df <- dplyr::bind_rows(q)
openxlsx::write.xlsx(df, here::here("data/dimensional_reduction/wnn/lisi_basw.xlsx"))

p10 <- df %>%
  tidyr::pivot_longer(c(-condition, -method), names_to = "test", values_to = "score") %>%
  dplyr::mutate(condition = factor(condition, levels = condition_levels),
                method = factor(method, levels = c("pca", "jointly", "svd", "liger")),
                test = factor(test, c("globaliLISI", "globalbASW")),
                modality = case_when(method %in% c("pca", "jointly") ~ "rna",
                                     method %in% c("svd", "liger") ~ "atac")) %>%
  ggplot2::ggplot(aes(x = condition, y = score, fill = method)) +
  ggplot2::geom_bar(stat="identity", position = position_dodge(), show.legend = TRUE) +
  ggplot2::facet_wrap(~ modality + test, ncol = 6) +
  my_theme() +
  ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ggplot2::coord_cartesian(clip = "off")

ggplot2::ggsave(plot = p10,
                file = here::here("data/dimensional_reduction/wnn/lisi_basw.pdf"),
                height = 2,
                width = 3,
                dpi = 1000)
