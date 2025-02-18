# Description -------------------------------------------------------------
# Cluster identification on liger embeddings performed on tiles.
# These clusters will be used for peak detection
# Setup -------------------------------------------------------------------
base::source(here::here("R/set_up.R"))
set.seed(1000)

# Setting a Genome and GeneAnnotation
ArchR::addArchRGenome("mm10")

create_directories(c(here::here("data/dimensional_reduction"), 
                     here::here("data/dimensional_reduction/atac")))

# Load --------------------------------------------------------------------
multi_islets_3 <- ArchR::loadArchRProject(path = here::here("archr_projects/save_multi_islets_3/"))
liger_tiles <- base::readRDS(here::here("data/dimensional_reduction/atac/liger_integration_5000_bp_tiles.rds"))
seurat_4 <- base::readRDS(here::here("data/seurat_objects/seurat_4.rds"))

# Find clusters and UMAP --------------------------------------------------
# Run lovastin clusterstering, using all 20 dimensions.
liger_tiles <- liger_tiles %>%
  rliger::louvainCluster(
    resolution = 2,
    k = 20,
    prune = 1/15,
    eps = 0.1,
    nRandomStarts = 10,
    nIterations = 100,
    random.seed = 1000,
    verbose = TRUE)


# Run umap
liger_tiles <- rliger::runUMAP(
  liger_tiles,
  use.raw = FALSE,
  dims.use = 1:ncol(liger_tiles@H.norm),
  k = 2,
  distance = "cosine",
  n_neighbors = 30,
  min_dist = 0.3,
  rand.seed = 42
)

# Add UMAP embeddings and liger clusters to archr project -----------------
# Approach from this thread: https://github.com/GreenleafLab/ArchR/issues/455
# "you create a new dataframe in the exact format shown below,
# importing dimension 1 and 2 of your embedding.
# You then add this to the embeddings present in your project."

# get UMAP embeddings
embeddings <- S4Vectors::DataFrame(liger_tiles@tsne.coords) %>%
  BiocGenerics::as.data.frame() %>%
  tibble::rownames_to_column("barcode")

# match order to barcodes
embeddings <- embeddings[base::match(BiocGenerics::rownames(multi_islets_3@cellColData), embeddings$barcode), ]

# add barcodes back as rownmaes

BiocGenerics::rownames(embeddings) <- NULL
embeddings <- embeddings %>%
  tibble::column_to_rownames("barcode")

# add new colnames
BiocGenerics::colnames(embeddings) <- c("ligertiles#UMAP_1",
                          "ligertiles#UMAP_2")

# Check that barcode order is equal between dataframe and archr project
base::all.equal(BiocGenerics::rownames(embeddings), BiocGenerics::rownames(multi_islets_3@cellColData))

# Add UMAP embeddings to archr
multi_islets_3@embeddings$liger_tiles_umap <- SimpleList(df = embeddings, params = list())

# Add clusters - join the two dataframes together to ensure that we are adding clusters to the right barcode.
# Extract liger clusters
liger_clusters <- liger_tiles@clusters %>%
  BiocGenerics::as.data.frame() %>%
  tibble::rownames_to_column("barcode") %>%
  dplyr::rename(liger_clusters_res_2 = ".")

# Combine liger clusters with cellcoldata from archr proejct
cellcoldata <- multi_islets_3@cellColData %>%
  BiocGenerics::as.data.frame() %>%
  tibble::rownames_to_column("barcode") %>%
  dplyr::left_join(y = liger_clusters,
                   by = "barcode") %>%
  dplyr::mutate(liger_clusters_res_2 = as.character(liger_clusters_res_2))

# add liger clusters to archr project using the cellcoldata dataframe
multi_islets_3<- ArchR::addCellColData(ArchRProj = multi_islets_3,
                                       data = cellcoldata$liger_clusters_res_2,
                                       cells = cellcoldata$barcode,
                                       name = "liger_clusters_res_2",
                                       force = TRUE)

# UMAP plots --------------------------------------------------------------
# Plot UMAP with ggplot - color by sample
tile_sample_plot <- multi_islets_3@cellColData %>%
  BiocGenerics::as.data.frame () %>%
  tibble::rownames_to_column("barcode") %>%
  dplyr::select(barcode, Sample, liger_clusters_res_2) %>%
  dplyr::left_join(multi_islets_3@embeddings$liger_tiles_umap$df %>%
                     BiocGenerics::as.data.frame() %>%
                     tibble::rownames_to_column("barcode"), by = "barcode") %>%
  dplyr::rename("UMAP1" = "ligertiles#UMAP_1", "UMAP2" = "ligertiles#UMAP_2") %>%
  dplyr::mutate(Sample = factor(Sample, levels = sample_levels)) %>%
  ggplot2::ggplot(aes(x = UMAP1,
                      y = UMAP2,
                      color = Sample))+
  ggplot2::scale_color_manual(values = sample_color) +
  ggplot2::geom_point(size = 0.1)+
  ggplot2::labs(x = "UMAP 1 (Tiles)",
                y = "UMAP 2 (Tiles)") +
  ggplot2::guides(color=guide_legend(ncol = 3,
                                     override.aes = list(size = 2))) +
  my_theme_void() +
  ggplot2::coord_cartesian(clip = "off") +
  ggplot2::theme(legend.title = element_blank(),
                 legend.position = c(.5, .8),
                 legend.justification = c("left", "bottom"),
                 legend.direction = "horizontal")

# Plot UMAP with ggplot - color by sample - split by sample
tile_sample_split_plot <- multi_islets_3@cellColData %>%
  BiocGenerics::as.data.frame () %>%
  tibble::rownames_to_column("barcode") %>%
  dplyr::select(barcode, Sample, liger_clusters_res_2) %>%
  dplyr::left_join(multi_islets_3@embeddings$liger_tiles_umap$df %>%
                     BiocGenerics::as.data.frame() %>%
                     tibble::rownames_to_column("barcode"), by = "barcode") %>%
  dplyr::rename("UMAP1" = "ligertiles#UMAP_1", "UMAP2" = "ligertiles#UMAP_2") %>%
  dplyr::mutate(Sample = factor(Sample, levels = sample_levels)) %>%
  ggplot2::ggplot(aes(x = UMAP1,
                      y = UMAP2,
                      color = Sample))+
  ggplot2::scale_color_manual(values = sample_color) +
  ggplot2::geom_point(size = 0.1)+
  ggplot2::labs(x = "UMAP 1 (Tiles)",
                y = "UMAP 2 (Tiles)") +
  ggplot2::facet_wrap(~ Sample, ncol = 7) +
  my_theme_void() +
  ggplot2::coord_cartesian(clip = "off") +
  ggplot2::theme(legend.position = "none")

# Plot UMAP with ggplot - color by sample - split by diet
tile_diet_split_plot <- multi_islets_3@cellColData %>%
  BiocGenerics::as.data.frame () %>%
  tibble::rownames_to_column("barcode") %>%
  dplyr::select(barcode, Sample, liger_clusters_res_2) %>%
  dplyr::left_join(multi_islets_3@embeddings$liger_tiles_umap$df %>%
                     BiocGenerics::as.data.frame() %>%
                     tibble::rownames_to_column("barcode"), by = "barcode") %>%
  dplyr::rename("UMAP1" = "ligertiles#UMAP_1", "UMAP2" = "ligertiles#UMAP_2") %>%
  dplyr::mutate(diet = gsub('.{3}$', '', Sample),
                diet = factor(diet, levels = c("LFD", "HFD_1", "HFD_3")),
                Sample = factor(Sample, levels = sample_levels)) %>%
  ggplot2::ggplot(aes(x = UMAP1,
                      y = UMAP2,
                      color = Sample))+
  ggplot2::scale_color_manual(values = sample_color) +
  ggplot2::geom_point(size = 0.1)+
  ggplot2::labs(x = "UMAP 1 (Tiles)",
                y = "UMAP 2 (Tiles)") +
  ggplot2::facet_wrap(~ diet, ncol = 3) +
  ggplot2::guides(color=guide_legend(ncol = 3,
                                     override.aes = list(size = 2))) +
  my_theme_void() +
  ggplot2::coord_cartesian(clip = "off") +
  ggplot2::theme(legend.title = element_blank(),
                 legend.position = c(.7, .7),
                 legend.justification = c("left", "bottom"),
                 legend.direction = "horizontal")

# Plot UMAP with ggplot - color by liger cluster
tile_cluster_plot <- multi_islets_3@cellColData %>%
  BiocGenerics::as.data.frame () %>%
  tibble::rownames_to_column("barcode") %>%
  dplyr::select(barcode, Sample, liger_clusters_res_2) %>%
  dplyr::left_join(multi_islets_3@embeddings$liger_tiles_umap$df %>%
                     BiocGenerics::as.data.frame() %>%
                     tibble::rownames_to_column("barcode"), by = "barcode") %>%
  dplyr::mutate(liger_clusters_res_2 = factor(liger_clusters_res_2,
                                              levels = as.character(c(0:100)))) %>%
  dplyr::rename("UMAP1" = "ligertiles#UMAP_1", "UMAP2" = "ligertiles#UMAP_2") %>%
  ggplot2::ggplot(aes(x = UMAP1,
                      y = UMAP2,
                      color = liger_clusters_res_2))+
  ggplot2::scale_color_manual(values = cluster_color) +
  ggplot2::geom_point(size = 0.1)+
  ggplot2::labs(x = "UMAP 1 (Tiles)",
                y = "UMAP 2 (Tiles)") +
  ggplot2::guides(color=guide_legend(ncol = 4,
                                     override.aes = list(size = 2))) +
  my_theme_void() +
  ggplot2::coord_cartesian(clip = "off") +
  ggplot2::theme(legend.title = element_blank(),
                 legend.position = c(.5, .7),
                 legend.justification = c("left", "bottom"),
                 legend.direction = "horizontal")

# create list of plots
plot_list <- base::list(tile_sample_plot, tile_cluster_plot)


# Save plots --------------------------------------------------------------
pdf(file = here::here("data/dimensional_reduction/atac/liger_tiles_integration_umap.pdf"),
    height = 6,
    width = 6)
plot_list
dev.off()

pdf(file = here::here("data/dimensional_reduction/atac/liger_tiles_integration_umap_diet.pdf"),
    height = 2,
    width = 6)
tile_diet_split_plot
dev.off()

pdf(file = here::here("data/dimensional_reduction/atac/liger_tiles_integration_umap_sample.pdf"),
    height = 0.85,
    width = 6)
tile_sample_split_plot
dev.off()


# Add Gene score ----------------------------------------------------------
# Add gene score matrix, default parameters expect for tile size is 5000
multi_islets_3 <-
  ArchR::addGeneScoreMatrix(
    input = multi_islets_3,
    force = TRUE,
    logFile = createLogFile("addGeneScoreMatrix")
  )

# Impute (smooth) gene scores with MAGIC
multi_islets_3 <- ArchR::addImputeWeights(multi_islets_3)

# Gene score plots --------------------------------------------------------
# Plot gene scores
marker_plot <- ArchR::plotEmbedding(
  ArchRProj = multi_islets_3,
  colorBy = "GeneScoreMatrix",
  name = purrr::as_vector(markers_short),
  embedding = "liger_tiles_umap",
  imputeWeights = getImputeWeights(multi_islets_3)
)

# Add new theme and title
marker_plot <- BiocGenerics::Map(archr_dim_gene_score,
                                 plot = marker_plot,
                                 title = names(marker_plot))
# Save plots --------------------------------------------------------------
pdf(file = here::here("data/dimensional_reduction/atac/liger_tiles_integration_gene_score_markers.pdf"),
    height = 9.18,
    width = 9.18)
marker_plot
dev.off()

# Gene expression + Gene score in seurat ----------------------------------
## Add liger clusters to seurat object
cellcoldata <- cellcoldata %>%
  mutate(barcode = cellcoldata$barcode %>% atac_to_rna_syntax())

# Add clusters  to meta data
seurat_4@meta.data <- seurat_4@meta.data %>%
  tibble::rownames_to_column("barcode") %>%
  full_join(y = cellcoldata,
            by = "barcode") %>%
  tibble::column_to_rownames("barcode")

# change liger clusters to numeric
seurat_4@meta.data$liger_clusters_res_2 <- as.numeric(seurat_4@meta.data$liger_clusters_res_2)

## Add gene scores to seurat object
# Extract gene score
gene_score <- ArchR::getMatrixFromProject(multi_islets_3, useMatrix='GeneScoreMatrix')

# Gene score matrix
gene_score_matrix <- assays(gene_score)$GeneScoreMatrix

# colnames
colnames(gene_score_matrix) <- colnames(gene_score_matrix) %>%
  atac_to_rna_syntax()

# row data
rowdata <- gene_score %>%
  SummarizedExperiment::rowData()
rownames(gene_score_matrix) <- rowdata$name

# Add gene score to new assay
seurat_4[["tiles"]] <- CreateAssayObject(counts = gene_score_matrix)

# Normalize gene scores
DefaultAssay(seurat_4) <- "tiles"
seurat_4 <- seurat_4 %>%
  Seurat::NormalizeData()

# Violin plot -------------------------------------------------------------
# Gene expression
gene_exp <-  Seurat::VlnPlot(seurat_4,
                             features = purrr::as_vector(markers_short),
                             stack = TRUE,
                             flip = TRUE,
                             fill.by = "ident",
                             assay = "RNA",
                             group.by = "liger_clusters_res_2",
                             cols = cluster_color) +
  ggplot2::labs(x = "Clusters",
                y = "ln(NormCounts +1)",
                title = "Gene Expression") +
  ggprism::theme_prism(border = T,
                       base_fontface = "plain",
                       base_size = 12) +
  ggplot2::theme(legend.position = "none") &
  ggplot2::geom_boxplot(width = 0.3,
                        outlier.shape = NA)

# Gene score
gene_sc <- Seurat::VlnPlot(seurat_4,
                           features = purrr::as_vector(markers_short),
                           stack = TRUE,
                           flip = TRUE,
                           fill.by = "ident",
                           assay = "tiles",
                           group.by = "liger_clusters_res_2",
                           cols = cluster_color) +
  ggplot2::labs(x = "Clusters",
                y = "ln(NormCounts +1)",
                title = "Gene Score") +
  ggprism::theme_prism(border = T,
                       base_fontface = "plain",
                       base_size = 12) +
  ggplot2::theme(legend.position = "none") &
  ggplot2::geom_boxplot(width = 0.3,
                        outlier.shape = NA)

# save plot
pdf(height = 20,
    width = 15,
    file = here::here("data/dimensional_reduction/atac/liger_tiles_integration_tiles_gene_score_exp_vln.pdf"))
gene_sc + gene_exp
dev.off()


# Add cluster identity ----------------------------------------------------
cellcoldata_clu <- multi_islets_3@cellColData %>%
  as.data.frame() %>%
  tibble::rownames_to_column("barcode") %>%
  dplyr::mutate(cluster_anno = dplyr::case_when(liger_clusters_res_2 == 0 ~ "Beta",
                                                liger_clusters_res_2 == 1 ~ "Beta",
                                                liger_clusters_res_2 == 2 ~ "Beta",
                                                liger_clusters_res_2 == 3 ~ "Beta",
                                                liger_clusters_res_2 == 4 ~ "Beta",
                                                liger_clusters_res_2 == 5 ~ "Beta",
                                                liger_clusters_res_2 == 6 ~ "Alpha",
                                                liger_clusters_res_2 == 7 ~ "Beta",
                                                liger_clusters_res_2 == 8 ~ "Beta",
                                                liger_clusters_res_2 == 9 ~ "Beta",
                                                liger_clusters_res_2 == 10 ~ "Beta",
                                                liger_clusters_res_2 == 11 ~ "Beta",
                                                liger_clusters_res_2 == 12 ~ "Beta",
                                                liger_clusters_res_2 == 13 ~ "Delta",
                                                liger_clusters_res_2 == 14 ~ "Beta",
                                                liger_clusters_res_2 == 15 ~ "Beta",
                                                liger_clusters_res_2 == 16 ~ "Alpha",
                                                liger_clusters_res_2 == 17 ~ "Beta",
                                                liger_clusters_res_2 == 18 ~ "Beta",
                                                liger_clusters_res_2 == 19 ~ "Beta",
                                                liger_clusters_res_2 == 20 ~ "Alpha",
                                                liger_clusters_res_2 == 21 ~ "Beta",
                                                liger_clusters_res_2 == 22 ~ "Endothelial",
                                                liger_clusters_res_2 == 23 ~ "Beta",
                                                liger_clusters_res_2 == 24 ~ "Gamma",
                                                liger_clusters_res_2 == 25 ~ "Delta",
                                                liger_clusters_res_2 == 26 ~ "Stellate",
                                                liger_clusters_res_2 == 27 ~ "Immune",
                                                liger_clusters_res_2 == 28 ~ "Acinar",
                                                liger_clusters_res_2 == 29 ~ "Stellate",
                                                liger_clusters_res_2 == 30 ~ "Endothelial",
                                                liger_clusters_res_2 == 31 ~ "Alpha",
                                                liger_clusters_res_2 == 32 ~ "Stellate",
                                                liger_clusters_res_2 == 33 ~ "Immune"
  ))

# add new annotation to archr project using the cellcoldata dataframe
multi_islets_3 <- ArchR::addCellColData(ArchRProj = multi_islets_3,
                                        data = cellcoldata_clu$cluster_anno,
                                        cells = cellcoldata_clu$barcode,
                                        name = "cluster_anno",
                                        force = TRUE)


# UMAP plot ---------------------------------------------------------------
tile_cluster_anno_plot <- multi_islets_3@cellColData %>%
  BiocGenerics::as.data.frame () %>%
  tibble::rownames_to_column("barcode") %>%
  dplyr::select(barcode, Sample, cluster_anno) %>%
  dplyr::left_join(multi_islets_3@embeddings$liger_tiles_umap$df %>%
                     BiocGenerics::as.data.frame() %>%
                     tibble::rownames_to_column("barcode"), by = "barcode") %>%
  dplyr::rename("UMAP1" = "ligertiles#UMAP_1", "UMAP2" = "ligertiles#UMAP_2") %>%
  ggplot2::ggplot(aes(x = UMAP1,
                      y = UMAP2,
                      color = cluster_anno))+
  ggplot2::scale_color_manual(values = cluster_anno) +
  ggplot2::geom_point(size = 0.1)+
  ggplot2::labs(x = "UMAP 1 (Tiles)",
                y = "UMAP 2 (Tiles)") +
  ggplot2::guides(color=guide_legend(ncol = 4,
                                     override.aes = list(size = 2))) +
  ggprism::theme_prism(border = T,
                       base_fontface = "plain",
                       base_size = 12) +
  ggplot2::coord_cartesian(clip = "off") +
  ggplot2::theme(legend.title = element_blank(),
                 legend.position = c(.5, .7),
                 legend.justification = c("left", "bottom"),
                 legend.direction = "horizontal")

# save plot
pdf(file = here::here("data/dimensional_reduction/atac/liger_tiles_integration_umap_cluster_anno.pdf"),
    height = 9.18,
    width = 9.18)
tile_cluster_anno_plot
dev.off()


# Add annotation to seurat ------------------------------------------------
## Add liger clusters to seurat object
cellcoldata_clu <- cellcoldata_clu %>%
  mutate(barcode = cellcoldata_clu$barcode %>% atac_to_rna_syntax()) %>%
  dplyr::select(barcode, cluster_anno)

# Add clusters  to meta data
seurat_4@meta.data <- seurat_4@meta.data %>%
  tibble::rownames_to_column("barcode") %>%
  left_join(y = cellcoldata_clu,
            by = "barcode") %>%
  tibble::column_to_rownames("barcode") %>%
  dplyr::mutate(cluster_anno = factor(cluster_anno, levels = c("Beta",
                                                               "Alpha",
                                                               "Delta",
                                                               "Gamma",
                                                               "Acinar",
                                                               "Endothelial",
                                                               "Stellate",
                                                               "Immune")))
# gene score + gene expresssion cluster anno ------------------------------
# Violin plot -------------------------------------------------------------
# Gene expression
gene_exp <-  Seurat::VlnPlot(seurat_4,
                             features = purrr::as_vector(markers_short),
                             stack = TRUE,
                             flip = TRUE,
                             fill.by = "ident",
                             assay = "RNA",
                             group.by = "cluster_anno",
                             cols = cluster_anno) +
  ggplot2::labs(x = "Clusters",
                y = "ln(NormCounts +1)",
                title = "Gene Expression") +
  ggprism::theme_prism(border = T,
                       base_fontface = "plain",
                       base_size = 12) +
  ggplot2::theme(legend.position = "none")  &
  ggplot2::geom_boxplot(width = 0.3,
                        outlier.shape = NA)


# Gene score
gene_sc <- Seurat::VlnPlot(seurat_4,
                           features = purrr::as_vector(markers_short),
                           stack = TRUE,
                           flip = TRUE,
                           fill.by = "ident",
                           assay = "tiles",
                           group.by = "cluster_anno",
                           cols = cluster_anno) +
  ggplot2::labs(x = "Clusters",
                y = "ln(NormCounts +1)",
                title = "Gene Score") +
  ggprism::theme_prism(border = T,
                       base_fontface = "plain",
                       base_size = 12) +
  ggplot2::theme(legend.position = "none")  &
  ggplot2::geom_boxplot(width = 0.3,
                        outlier.shape = NA)

# save plot
pdf(height = 20,
    width = 9.18,
    file = here::here("data/dimensional_reduction/atac/liger_tiles_integration_cluster_anno_gene_score_exp_vln.pdf"))
gene_sc + gene_exp
dev.off()


# save --------------------------------------------------------------------
# save liger tiles
saveRDS(liger_tiles, file = here::here("data/dimensional_reduction/atac/liger_integration_5000_bp_tiles.rds"))

# Save Archr project
saveArchRProject(ArchRProj = multi_islets_3,
                 outputDirectory = here::here("archr_projects/save_multi_islets_3"),
                 load = FALSE)

# save seurat obj
saveRDS(seurat_4, file = here::here("data/dimensional_reduction/atac/liger_tiles_integration_seurat_4_rna.rds"))

