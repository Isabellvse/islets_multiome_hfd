# Description -------------------------------------------------------------
# In this section we will keep barcodes which has passed quality control in both
# RNA and ATAC. The seurat object will be subset to contain these cells.
# hereafter I will integrate cells using JOINTLY
# here after we will use this cells to identify polyhormone cells, which could be doublets

# Setup -------------------------------------------------------------------
base::source(here::here("R/set_up.R"))
# we don't set seed yet

# Create diretories
create_directories(c(here::here("data/dimensional_reduction"), 
                     here::here("data/dimensional_reduction/rna"),
                     here::here("data/export")))

# Load --------------------------------------------------------------------
# Seurat object list
seurat_list <- base::readRDS(file = here::here("data/seurat_objects/seurat_object_1.rds"))

# Barcodes to keep
bar_keep <-  base::readRDS(here::here("data/quality_control/rna_atac_barcode_keep.rds"))

# Subset seurat object ----------------------------------------------------
# Create a list of barcodes to keep
bar_keep_list <- bar_keep %>%
  BiocGenerics::as.data.frame() %>%
  dplyr::rename("barcode" = ".") %>%
  dplyr::mutate(sample = gsub( "_[^_]*$", "", barcode),
                sample = factor(sample, levels = sample_levels)) %>%
  dplyr::mutate(split_by = sample) %>%
  collapse::rsplit(~ split_by) %>% 
  purrr::modify_depth(1, ~ dplyr::pull(., barcode))

# subset seurat object
seurat_list_sub <- purrr::map2(seurat_list, base::names(seurat_list), function(s_obj, names){
  output <- base::subset(s_obj,
                         cells = bar_keep_list[[names]])
  return(output)
})

# Prepare data for batch correction suing JOINTLY -------------------------
# Need a seurat object containing raw counts
# Remove everying except RNA assay
seurat_list_sub <- purrr::map(seurat_list_sub, function(s_obj){
  
  Seurat::DefaultAssay(s_obj) <- "RNA"
  
  output <- Seurat::DietSeurat(object = s_obj,
                               assays = "RNA")
  return(output)
  
})

# Merge seurat objects
seurat_2 <- purrr::reduce(seurat_list_sub, base::merge)

# clean up metadata and add condition
seurat_2@meta.data <- seurat_2@meta.data %>%
  BiocGenerics::as.data.frame() %>%
  dplyr::select(-seurat_clusters,
                - RNA_snn_res.0.7,
                - dplyr::starts_with("pANN"),
                - dplyr::starts_with("DF.classifications")) %>%
  dplyr::mutate(condition = gsub( "_[^_]*$", "", orig.ident),
                condition = factor(condition, levels = condition_levels),
                orig.ident = factor(orig.ident, levels = sample_levels))

# save seurat object
saveRDS(seurat_2, file = here::here("data/seurat_objects/seurat_object_2.rds"))

# Run JOUNTLY -------------------------------------------------------------
# Run Jointly
proc <- JOINTLY::preprocess(data = seurat_2,
                            batch.var = "orig.ident")

cpca <- JOINTLY::cpca(dataset.list = proc,
                      bpparam = BiocParallel::MulticoreParam())

inputs <- JOINTLY::prepareData(dataset.list = cpca$cpca)

solved <- list()

# Do jointly 5 time, each time with a new seed
repeat {
  set.seed(Sys.time())
  new_element <- JOINTLY::JOINTLYsolve(
    kernel.list = inputs$kernels,
    snn.list = inputs$snn,
    rare.list = inputs$rareity,
    cpca.result = cpca,
    k = 15,
    bpparam = BiocParallel::MulticoreParam()
  )
  solved[[length(solved) +1 ]] <- new_element
  
  if (length(solved) == 5)
    break
}

## save objects ----
base::saveRDS(proc, here::here("data/dimensional_reduction/rna/jointly_proc_seurat_2.rds"))
base::saveRDS(cpca, here::here("data/dimensional_reduction/rna/jointly_cpca_seurat_2.rds"))
base::saveRDS(inputs, here::here("data/dimensional_reduction/rna/jointly_inputs_seurat_2.rds"))
base::saveRDS(solved, here::here("data/dimensional_reduction/rna/jointly_5_itr_solved_seurat_2.rds"))

# set seed again ----------------------------------------------------------
set.seed(1000)

# Normalize RNA data ------------------------------------------------------
seurat_2 <- seurat_2 %>%
  Seurat::NormalizeData()

# Add jointly embedding to umap -------------------------------------------

# Add embedding to seurat object and run umap
s_list <- base::list()
for (i in 1:length(solved)){
  H <- solved[[i]]$Hmat.scaled
  colnames(H) <- paste("JOINTLY", 1:ncol(H), sep="_")
  seurat_2[["jointly"]] <- Seurat::CreateDimReducObject(as.matrix(H), assay = "RNA")
  new_element <- Seurat::RunUMAP(seurat_2,
                                 reduction = "jointly",
                                 dims = 1:15)
  s_list[[length(s_list) + 1]] <- new_element
}

# plot UMAPS
p1 <- purrr::map(s_list, function(s_obj){
  output <- s_obj %>% Seurat::DimPlot(reduction = "umap",
                                      group.by = "condition",
                                      cols = condition_color, 
                                      pt.size = 3, 
                                      raster = TRUE, 
                                      raster.dpi = c(800, 800)) +
    my_theme_void() +
    ggplot2::theme(legend.position = "none", 
                   title = element_blank())
  return(output)
})

# wrap umaps
p1 <- patchwork::wrap_plots(p1, ncol = 5)

# save plot
ggplot2::ggsave(p1, filename = here::here("data/dimensional_reduction/rna/jointly_umap_5_itr_seurat_2.pdf"),
                height = 1.2,
                width = 6)


# Evaluate embedding ------------------------------------------------------
## split by diet ----
metadata <-  purrr::map(s_list, ~ .x@meta.data %>% 
                          dplyr::mutate(., split_by = condition) %>% 
                          collapse::rsplit(~ split_by))
embeddings <- purrr::map(s_list, ~ .x@reductions$jointly@cell.embeddings) 

embeddings_sub <- purrr::map2(embeddings, metadata, function(embedding, metadata_list) {
  purrr::map(metadata_list, function(md) {
    rownames_md <- rownames(md)
    embedding[rownames_md, , drop = FALSE]
  })
  })

# Define variables
batch_var <- "orig.ident"

# Apply the ilisi_comp function to each pair of metadata and embeddings
results <- purrr::map2(embeddings_sub, metadata,~ purrr::map2(.x,.y,
                ~ ilisi_comp(.x,.y,dims = 15, batch_var = batch_var)))

# Extract results
global_ilisi <- purrr::map_dfr(seq_along(results), function(i) {
  purrr::map_dfr(base::names(results[[i]]), function(name) {
    tibble::tibble(iteration = i, condition = name, global_iLISI = results[[i]][[name]]$global_iLISI)
  })
}) 

## Save glocal results ----
openxlsx::write.xlsx(global_ilisi, here::here("data/export/ilisi_jointly.xlsx"))

## save seurat list object ----
base::saveRDS(s_list, file = here::here("data/dimensional_reduction/rna/s_list.rds"))

# closer look at iteration 5 ---------------------------------------------
# extract seurat object of interest
seurat_3 <- s_list[[5]]

## change these plots 
p2 <- seurat_3 %>%
  Seurat::DimPlot(reduction = "umap",
                  split.by = "condition",
                  group.by = "orig.ident",
                  cols = sample_color,
                  pt.size = 5, 
                  raster = TRUE, 
                  raster.dpi = c(800, 800)) +
  my_theme_void(remove_strip_text = FALSE) +
  ggplot2::theme(title = element_blank(),
                 legend.position = "none")



p3 <- seurat_3 %>%
  Seurat::DimPlot(reduction = "umap",
                  split.by = "orig.ident",
                  group.by = "orig.ident",
                  cols = sample_color,
                  pt.size = 5, 
                  raster = TRUE, 
                  raster.dpi = c(800, 800)) +
  my_theme_void(remove_strip_text = FALSE) +
  ggplot2::theme(title = element_blank(),
                 legend.position = "none")

# save plots
ggplot2::ggsave(p2, filename = here::here("data/dimensional_reduction/rna/jointly_umap_itr_nr_5_split_by_condition.pdf"),
                height = 6/3,
                width = 6,
                dpi = 500)
ggplot2::ggsave(p3, filename = here::here("data/dimensional_reduction/rna/jointly_umap_itr_nr_5_split_by_sample.pdf"),
                height = 6/6,
                width = 6, dpi = 500)

# Interpret W matrix ------------------------------------------------------
# Interpret W matrix - for iteration nr 5

W <- solved[[5]]$Wmat

for (i in 1:length(W)) {
  W.tmp <- W[[i]]
  W.tmp <- scale(W.tmp)
  W.tmp <- t(W.tmp)
  W.tmp <- scale(W.tmp)
  W.tmp <- t(W.tmp)
  W[[i]] <- W.tmp
}

# Sum the W matrix
W.sum <- W[[1]]
for (i in 2:length(W)) {
  W.sum <- W.sum + W[[i]]
}

W.sum <- W.sum / length(solved[[5]]$Wmat)

rownames(W.sum) <- rownames(cpca$normalized[[1]])

# Scale the sum matrix
W.tmp <- W.sum
W.tmp <- scale(W.tmp)
W.tmp <- t(W.tmp)
W.tmp <- scale(W.tmp)
W.tmp <- t(W.tmp)
W.sum <- W.tmp

# Get module genes
modules <- list()
for (i in 1:15) {
  modules[[length(modules) + 1]] <-
    base::names(sort(W.sum[, i], decreasing = TRUE))[1:inflection::uik(y = sort(W.sum[, i], decreasing = TRUE),
                                                                 x = seq(1, nrow(W.sum), 1))]
  base::names(modules)[base::length(modules)] <-
    paste("factor_", i, sep = "")
}

## Save modules ----
base::saveRDS(modules, file = here::here("data/dimensional_reduction/rna/jointly_modules_itr_nr_5_seurat_2.rds"))

# Deep clustering for doublet identificaiton ------------------------------
# Find clusters
seurat_3 <- seurat_3 %>%
  Seurat::FindNeighbors(reduction = "jointly",
                        dims = 1:15) %>%
  Seurat::FindClusters(resolution = 20,
                       algorithm = 1)

# plot clusteirng
p4 <- seurat_3 %>% Seurat::DimPlot(reduction = "umap",
                                   group.by = "seurat_clusters",
                                   label = FALSE,
                                   pt.size = 3, 
                                   raster = TRUE, 
                                   raster.dpi = c(800, 800)) +
  my_theme_void(remove_strip_text = FALSE) +
  ggplot2::theme(title = element_blank(),
                 legend.position = "none")

# save cluster
ggplot2::ggsave(p4,
                filename = here::here("data/quality_control/rna/ALL_dimplot_clusters_used_for_doublet_identification.pdf"),
                height = 6,
                width = 6,
                dpi = 500)


# Polyhormone detection - marker genes ------------------------------------

# Define marker genes to use for polyhormone detection
markers <- c("Ins2", "Sst", "Ppy", "Gcg")

# find doublet clusters
# get average expression per cluster
avg <- Seurat::AverageExpression(seurat_3, assay = "RNA",
                                 group.by = "RNA_snn_res.20")$RNA

# Scale gene expression column-wise (genes in columns)
avg.scaled <- t(scale(t(avg)))

# get polyhormone clusters

# Get scaled expression of canonical marker genes, and keep only clusters (columns) which have an
# average scaled expression above 0.6 (the expression is 0.6 standard deviation to the right of the mean on a bell curve (normal distribution).
# if the sum of these expression values within a cluster (column) is above 1 it means they expression
# more than 1 conical marker to a high degree, and thus could be doublets.
db_cluster <-
  names(which(colSums(avg.scaled[rownames(avg.scaled) %in% markers, ] > 0.6) > 1))

# save expression of marker genes
avg.scaled[rownames(avg.scaled) %in% markers, ] %>%
  as.data.frame() %>%
  dplyr::select(all_of(db_cluster))

# Plot doublets -----------------------------------------------------------
seurat_3@meta.data <- seurat_3@meta.data %>% 
  dplyr::mutate(multiplet = dplyr::case_when(RNA_snn_res.20 %in% db_cluster ~ "multiplet",
                                             !RNA_snn_res.20 %in% db_cluster ~ "singlet"))
p5 <- seurat_3 %>%
  Seurat::DimPlot(reduction = "umap",
                  group.by = "multiplet",
                  pt.size = 3,
                  cols = c("multiplet" = "#A83708",
                           "singlet" = "grey"),
                  order = TRUE,
                  raster = TRUE, 
                  raster.dpi = c(800, 800)) +
  my_theme_void(remove_strip_text = FALSE) +
  ggplot2::theme(title = element_blank(),
                 legend.position = "none")

ggplot2::ggsave(p5,
                filename = here::here("data/quality_control/rna/ALL_dimplot_doublet_identification.pdf"),
                height = 6,
                width = 6,
                dpi = 500)

# Plot feature plot of expression of marker genes
p6 <- seurat_3 %>%
  Seurat::FeaturePlot(features = markers,
                      pt.size = 3, 
                      raster = TRUE, 
                      raster.dpi = c(800, 800), 
                      cols = base::rev(RColorBrewer::brewer.pal(n = 11, name = "Spectral")),
                      ncol = 4) +
  patchwork::plot_layout(guides = "collect") &
  my_theme_void(remove_strip_text = FALSE) &
  ggplot2::guides(colour = ggplot2::guide_colourbar(barwidth = 0.5, barheight = 6))

# save plot
ggplot2::ggsave(p6,
                filename = here::here("data/quality_control/rna/ALL_feature_plot_marker_genes_used_for_doublet_identification.pdf"),
                height = 6/4,
                width = 6,
                dpi = 500)


# Calculate module scores -------------------------------------------------
seurat_3 <- UCell::AddModuleScore_UCell(seurat_3, features = modules)

# Plot module scores ------------------------------------------------------

# Feature plot
# Get names of modules
seurat_3 %>%
  Seurat::FeaturePlot(
    reduction = "umap",
    features =  base::names(modules) %>% stringi::stri_join("_UCell"),
    pt.size = 3, 
    raster = TRUE, 
    raster.dpi = c(800, 800), 
    cols = base::rev(RColorBrewer::brewer.pal(n = 11, name = "Spectral")),
    ncol = 5) +
  patchwork::plot_layout(guides = "collect") &
  my_theme_void(remove_strip_text = FALSE) &
  ggplot2::guides(colour = ggplot2::guide_colourbar(barwidth = 0.5, barheight = 6))

# Find clusters with lower resolution -------------------------------------
# Find clusters
seurat_3 <- seurat_3 %>%
  Seurat::FindClusters(resolution = 0.2,
                       algorithm = 1)

seurat_3 %>% Seurat::DimPlot(reduction = "umap", group.by = "seurat_clusters", cols = cluster_color) +
  my_theme_void(remove_strip_text = TRUE)


# Plot modules scores in heatmap ------------------------------------------
# extract module scores - and find median
module_df <- seurat_3@meta.data %>%
  BiocGenerics::as.data.frame() %>%
  dplyr::select(starts_with("factor_"), seurat_clusters) %>%
  dplyr::group_by(seurat_clusters) %>%
  dplyr::summarise(across(starts_with("factor_"), median)) %>%
  dplyr::rename_with(snakecase::to_title_case, starts_with("factor_")) %>%
  tibble::column_to_rownames("seurat_clusters") %>%
  dplyr::rename_all(~stringr::str_remove(., "u Cell"))

anno_row <- module_df %>%
  dplyr::mutate(seurat_clusters = BiocGenerics::rownames(module_df)) %>%
  dplyr::select(seurat_clusters)

module_df %>%
  pheatmap::pheatmap(
    scale = "row",
    color = rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")),
    annotation_row = anno_row,
    annotation_colors = list("seurat_clusters" = cluster_color),
    annotation_legend = FALSE,
    border_color = "black")

saveRDS(seurat_3, here::here("data/seurat_objects/seurat_3.rds"))

# Remove doublets ---------------------------------------------------------
db_remove <- seurat_3@meta.data %>%
  dplyr::filter(RNA_snn_res.20 %in% db_cluster) %>%
  BiocGenerics::rownames()

# remove doublet cells
seurat_4 <- base::subset(seurat_3,
                         cells = db_remove,
                         invert = TRUE)

# redo normalization
seurat_4 <- seurat_4 %>% Seurat::NormalizeData()

# redo umap on jointly embeddings
seurat_4 <- seurat_4 %>%
  Seurat::RunUMAP(reduction = "jointly",
                  dims = 1:15)

# Dimplots ----------------------------------------------------------------

p7 <- seurat_4 %>%
  Seurat::DimPlot(reduction = "umap",
                  group.by = "orig.ident",
                  cols = sample_color,
                  pt.size = 3, 
                  raster = TRUE, 
                  raster.dpi = c(800, 800)) +
  my_theme_void(remove_strip_text = TRUE)


p8 <- seurat_4 %>%
  Seurat::DimPlot(reduction = "umap",
                  split.by = "condition",
                  group.by = "orig.ident",
                  cols = sample_color,
                  pt.size = 3, 
                  raster = TRUE, 
                  raster.dpi = c(800, 800)) +
  my_theme_void(remove_strip_text = FALSE) +
  ggplot2::theme(legend.position = "none",
                 title = element_blank())


p9 <- seurat_4 %>%
  Seurat::DimPlot(reduction = "umap",
                  split.by = "orig.ident",
                  group.by = "orig.ident",
                  cols = sample_color,
                  pt.size = 3, 
                  raster = TRUE, 
                  raster.dpi = c(800, 800)) +
  my_theme_void(remove_strip_text = FALSE) +
  ggplot2::theme(legend.position = "none",
                 title = element_blank())

# save plots
ggplot2::ggsave(p7, filename = here::here("data/dimensional_reduction/rna/seurat4_dimplot.pdf"),
                height = 6,
                width = 6)
ggplot2::ggsave(p8, filename = here::here("data/dimensional_reduction/rna/seurat4_dimplot_condition_split_sample_color.pdf"),
                height = 6/3,
                width = 6)
ggplot2::ggsave(p9, filename = here::here("data/dimensional_reduction/rna/seurat4_dimplot_sample_split_sample_color.pdf"),
                height = 6/7,
                width = 6)

# save
base::saveRDS(seurat_4, file = here::here("data/seurat_objects/seurat_4.rds"))
base::saveRDS(db_remove, file = here::here("data/quality_control/rna/ALL_polyhormone_cells_remove.rds"))
