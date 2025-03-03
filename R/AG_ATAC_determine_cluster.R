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
multi_islets_3 <- ArchR::loadArchRProject(path = here::here("data/archr_projects/save_multi_islets_3/"))
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
multi_islets_3@embeddings$liger_tiles_umap <- S4Vectors::SimpleList(df = embeddings, params = list())

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
  ggrastr::geom_point_rast(size = 0.1, raster.dpi = 500)+
  ggplot2::labs(x = "UMAP 1 (Tiles)",
                y = "UMAP 2 (Tiles)") +
  ggplot2::guides(color=guide_legend(ncol = 3,
                                     override.aes = list(size = 2))) +
  my_theme_void() +
  ggplot2::theme(title = element_blank(),
                 legend.position = "none")
ggplot2::ggsave(tile_sample_plot, filename = here::here("data/dimensional_reduction/atac/liger_tiles_integration_umap.pdf"),
                height = 6,
                width = 6,
                dpi = 500)

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
  ggrastr::geom_point_rast(size = 1, raster.dpi = 500, scale = 0.2)+
  ggplot2::labs(x = "UMAP 1 (Tiles)",
                y = "UMAP 2 (Tiles)") +
  ggplot2::facet_wrap(~ Sample, ncol = 7) +
  my_theme_void(remove_strip_text = FALSE) +
  ggplot2::theme(title = element_blank(),
                 legend.position = "none")

ggplot2::ggsave(tile_sample_split_plot, filename = here::here("data/dimensional_reduction/atac/tile_sample_split_umap.pdf"),
                height = 6/7,
                width = 6,
                dpi = 500)

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
  ggrastr::geom_point_rast(size = 2, raster.dpi = 500, scale = 0.2) +
  ggplot2::labs(x = "UMAP 1 (Tiles)",
                y = "UMAP 2 (Tiles)") +
  ggplot2::facet_wrap(~ diet, ncol = 3) +
  ggplot2::guides(color=guide_legend(ncol = 3,
                                     override.aes = list(size = 2))) +
  my_theme_void(remove_strip_text = FALSE) +
  ggplot2::theme(title = element_blank(),
                 legend.position = "none")

ggplot2::ggsave(tile_diet_split_plot, filename = here::here("data/dimensional_reduction/atac/tile_diet_split_umap.pdf"),
                height = 6/3,
                width = 6,
                dpi = 500)


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
  ggrastr::geom_point_rast(size = 0.1, raster.dpi = 500)+
  ggplot2::labs(x = "UMAP 1 (Tiles)",
                y = "UMAP 2 (Tiles)") +
  ggplot2::guides(color=guide_legend(ncol = 4,
                                     override.aes = list(size = 2))) +
  my_theme_void(remove_strip_text = FALSE) +
  ggplot2::theme(title = element_blank(),
                 legend.position = "none")

ggplot2::ggsave(tile_cluster_plot, filename = here::here("data/dimensional_reduction/atac/tile_cluster_plot _umap.pdf"),
                height = 6,
                width = 6,
                dpi = 500)

# Add Gene score ----------------------------------------------------------
# Add gene score matrix, default parameters expect for tile size is 5000
multi_islets_3 <-
  ArchR::addGeneScoreMatrix(
    input = multi_islets_3,
    force = TRUE, 
    threads = parallel::detectCores() - 1,
    logFile = createLogFile("addGeneScoreMatrix")
  )

# Impute (smooth) gene scores with MAGIC
multi_islets_3 <- ArchR::addImputeWeights(multi_islets_3, threads = parallel::detectCores() - 1)

# Gene score plots --------------------------------------------------------
# Plot gene scores
marker_plot <- ArchR::plotEmbedding(
  ArchRProj = multi_islets_3,
  threads = parallel::detectCores() - 1,
  colorBy = "GeneScoreMatrix",
  name = purrr::as_vector(markers_short),
  embedding = "liger_tiles_umap",
  imputeWeights = ArchR::getImputeWeights(multi_islets_3)
)

# Add new theme and title
marker_plot <- BiocGenerics::Map(archr_dim_gene_score,
                                 plot = marker_plot,
                                 title = base::names(marker_plot))
# Save plots --------------------------------------------------------------
pdf(file = here::here("data/dimensional_reduction/atac/liger_tiles_integration_gene_score_markers.pdf"),
    height = 6,
    width = 6)
marker_plot
dev.off()

# Gene expression + Gene score in seurat ----------------------------------
## Add liger clusters to seurat object
cellcoldata <- cellcoldata %>%
  dplyr::mutate(barcode = cellcoldata$barcode %>% atac_to_rna_syntax())

# Add clusters  to meta data
seurat_4@meta.data <- seurat_4@meta.data %>%
  tibble::rownames_to_column("barcode") %>%
  dplyr::full_join(y = cellcoldata,
            by = "barcode") %>%
  tibble::column_to_rownames("barcode")

# change liger clusters to numeric
seurat_4@meta.data$liger_clusters_res_2 <- base::as.numeric(seurat_4@meta.data$liger_clusters_res_2)

## Add gene scores to seurat object
# Extract gene score
gene_score <- ArchR::getMatrixFromProject(multi_islets_3, useMatrix='GeneScoreMatrix')

# Gene score matrix
gene_score_matrix <- SummarizedExperiment::assays(gene_score)$GeneScoreMatrix

# colnames
colnames(gene_score_matrix) <- BiocGenerics::colnames(gene_score_matrix) %>%
  atac_to_rna_syntax()

# row data
rowdata <- gene_score %>%
  SummarizedExperiment::rowData()

BiocGenerics::rownames(gene_score_matrix) <- rowdata$name

# Add gene score to new assay
seurat_4[["tiles"]] <- Seurat::CreateAssayObject(counts = gene_score_matrix)

# Normalize gene scores
Seurat::DefaultAssay(seurat_4) <- "tiles"
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
  my_theme() +
  ggplot2::theme(title = element_blank(), 
                 legend.position = "none")

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
  my_theme() +
  ggplot2::theme(title = element_blank(),
                 legend.position = "none")


gene_exp_sc <- gene_sc + gene_exp
ggplot2::ggsave(gene_exp_sc, filename=here::here("data/dimensional_reduction/atac/liger_tiles_integration_tiles_gene_score_exp_vln.pdf"),
       width = 6,
       height = 8, 
       dpi = 500)

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
                                                liger_clusters_res_2 == 13 ~ "Beta",
                                                liger_clusters_res_2 == 14 ~ "Beta",
                                                liger_clusters_res_2 == 15 ~ "Delta",
                                                liger_clusters_res_2 == 16 ~ "Beta",
                                                liger_clusters_res_2 == 17 ~ "Alpha",
                                                liger_clusters_res_2 == 18 ~ "Beta",
                                                liger_clusters_res_2 == 19 ~ "Beta",
                                                liger_clusters_res_2 == 20 ~ "Beta",
                                                liger_clusters_res_2 == 21 ~ "Alpha",
                                                liger_clusters_res_2 == 22 ~ "Gamma",
                                                liger_clusters_res_2 == 23 ~ "Endothelial",
                                                liger_clusters_res_2 == 24 ~ "Beta",
                                                liger_clusters_res_2 == 25 ~ "Delta",
                                                liger_clusters_res_2 == 26 ~ "Endothelial",
                                                liger_clusters_res_2 == 27 ~ "Endothelial",
                                                liger_clusters_res_2 == 28 ~ "Immune",
                                                liger_clusters_res_2 == 29 ~ "Acinar",
                                                liger_clusters_res_2 == 30 ~ "Beta",
                                                liger_clusters_res_2 == 31 ~ "Stellate",
                                                liger_clusters_res_2 == 32 ~ "Beta",
                                                liger_clusters_res_2 == 33 ~ "Beta",
                                                liger_clusters_res_2 == 34 ~ "Immune"
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
  ggrastr::geom_point_rast(size = 0.1, raster.dpi = 500)+
  ggplot2::labs(x = "UMAP 1 (Tiles)",
                y = "UMAP 2 (Tiles)") +
  my_theme_void(remove_strip_text = FALSE) +
  ggplot2::theme(title = element_blank(),
                 legend.position = "right")

ggplot2::ggsave(tile_cluster_anno_plot, filename = here::here("data/dimensional_reduction/atac/liger_tiles_integration_umap_cluster_anno.pdf"),
                height = 6,
                width = 6,
                dpi = 500)

# Add annotation to seurat ------------------------------------------------
## Add liger clusters to seurat object
cellcoldata_clu <- cellcoldata_clu %>%
  dplyr::mutate(barcode = cellcoldata_clu$barcode %>% atac_to_rna_syntax()) %>%
  dplyr::select(barcode, cluster_anno)

# Add clusters  to meta data
seurat_4@meta.data <- seurat_4@meta.data %>%
  tibble::rownames_to_column("barcode") %>%
  dplyr::left_join(y = cellcoldata_clu,
            by = "barcode") %>%
  tibble::column_to_rownames("barcode") %>%
  dplyr::mutate(cluster_anno = base::factor(cluster_anno, levels = c("Beta",
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
  my_theme() +
  ggplot2::theme(title = element_blank(),
                 legend.position = "none")

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
  my_theme() +
  ggplot2::theme(title = element_blank(),
                 legend.position = "none")
# save plot
pdf(height = 20,
    width = 9.18,
    file = here::here("data/dimensional_reduction/atac/liger_tiles_integration_cluster_anno_gene_score_exp_vln.pdf"))
gene_sc + gene_exp
dev.off()


# save --------------------------------------------------------------------
# save liger tiles
base::saveRDS(liger_tiles, file = here::here("data/dimensional_reduction/atac/liger_integration_5000_bp_tiles.rds"))

# Save Archr project
ArchR::saveArchRProject(ArchRProj = multi_islets_3,
                 outputDirectory = here::here("data/archr_projects/save_multi_islets_3"),
                 threads = parallel::detectCores() - 1,
                 load = FALSE)

# save seurat obj
base::saveRDS(seurat_4, file = here::here("data/dimensional_reduction/atac/liger_tiles_integration_seurat_4_rna.rds"))

