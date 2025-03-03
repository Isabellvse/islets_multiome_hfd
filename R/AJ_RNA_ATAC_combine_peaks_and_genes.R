# Description -------------------------------------------------------------
# Combine RNA and ATAC
# Setup -------------------------------------------------------------------
base::source(here::here("R/set_up.R"))
set.seed(1000)

# Setting a Genome and GeneAnnotation
ArchR::addArchRGenome("mm10")

create_directories(here::here("data/dimensional_reduction/wnn/"))
# Load --------------------------------------------------------------------
multi_islets_4 <- ArchR::loadArchRProject(path = here::here("data/archr_projects/save_multi_islets_4/"))

# liger object
liger_peaks <- base::readRDS(here::here("data/dimensional_reduction/atac/liger_integration_peaks.rds"))

# seurat object with rna data
seurat_5 <- base::readRDS(here::here("data/seurat_objects/seurat_4.rds"))


# Fragment directory ------------------------------------------------------
#  A path to the cellranger-ATAC output.
# The directory contains all samples’ folders before '/outs/fragments.tsv.gz'.
fragments_dir <- here::here("data/cellranger_fragments/")

# Annotations -------------------------------------------------------------
annotation <- get_annotation(ensdb = EnsDb.Mmusculus.v79,
                             genome = "mm10")
# Convert archr to signac seurat object -----------------------------------

# Extract peak matrix
peak_matrix <- ArchRtoSignac::getPeakMatrix(ArchRProject= multi_islets_4)

# Create seurat object with atac assay
seurat_atac <- ArchR2Signac_multiome(
  ArchRProject = multi_islets_4,
  samples = sample_levels,
  fragments_dir = fragments_dir,
  pm = peak_matrix,
  output_dir = '/outs/',
  refversion = 'mm10',
  annotation = annotation
)

# save
base::saveRDS(seurat_atac, file = here::here("data/seurat_objects/seurat_atac.rds"))

# rename cells to fit rna syntax
seurat_atac <- Seurat::RenameCells(seurat_atac, new.names = BiocGenerics::colnames(seurat_atac) %>% stringr::str_remove("-1"))

# Combine meta-data -------------------------------------------------------
rna_meta <- seurat_5@meta.data %>%
  tibble::rownames_to_column("barcode")

atac_meta <- seurat_atac@meta.data %>%
  tibble::rownames_to_column("barcode") %>%
  dplyr::select(-"orig.ident", -condition)

# combined meta_data
meta_combined <- dplyr::left_join(rna_meta, atac_meta, by = "barcode") %>%
  tibble::column_to_rownames("barcode")

# Add meta data
seurat_5@meta.data <- meta_combined

# Add peaks assay ---------------------------------------------------------
seurat_5[["peaks"]] <- seurat_atac[["peaks"]]

# Add liger embedding -----------------------------------------------------
# Change name of embeddings to fit rna format
liger_h_norm <- liger_peaks@H.norm %>% BiocGenerics::as.data.frame() %>%
  tibble::rownames_to_column("barcode") %>%
  dplyr::mutate(barcode = atac_to_rna_syntax(barcode))

# Same order as seurat object
liger_h_norm_order <- liger_h_norm[BiocGenerics::match(BiocGenerics::colnames(seurat_5), liger_h_norm$barcode), ]

# Check that cell names are in the right order
base::all.equal(BiocGenerics::colnames(seurat_5), liger_h_norm_order$barcode)

# Add barcode back to rowname and create a matrix
BiocGenerics::rownames(liger_h_norm_order) <- NULL
liger_h_norm_order <- liger_h_norm_order %>%
  tibble::column_to_rownames("barcode") %>%
  base::as.matrix()

# Check order of barcodes again again
base::all.equal(BiocGenerics::colnames(seurat_5), BiocGenerics::rownames(liger_h_norm_order))

# Also change colnames
BiocGenerics::colnames(liger_h_norm_order) <- base::paste0("ligerembedding_", base::rep(1:20))

# Extract loadings and also change colnames
liger_loadings <- t(liger_peaks@W)
BiocGenerics::colnames(liger_loadings) <- base::paste0("ligerembedding_", base::rep(1:20))

# Extract liger embeddings
liger_embeddings <- Seurat::CreateDimReducObject(embeddings = liger_h_norm_order,
                                         key = "ligerembedding_",
                                         assay = "peaks",
                                         loadings = liger_loadings)

# add embedding
seurat_5@reductions[["liger_embeddings"]] <- liger_embeddings

# Save seurat  ------------------------------------------------------------
base::saveRDS(seurat_5, file = here::here("data/seurat_objects/seurat_5.rds"))

# UMAP liger embedding ----------------------------------------------------
seurat_5 <- seurat_5 %>%
  Seurat::RunUMAP(reduction = "liger_embeddings",
                  assay = "peaks",
                  dims = 1:20,
                  reduction.name = 'umap.atac',
                  reduction.key = 'atacUMAP_',
                  seed.use = 1000,
                  verbose = TRUE)


# Plot atac dim plot - to evaluate integration ----------------------------
p1 <- seurat_5 %>%
  Seurat::DimPlot(reduction = "umap.atac",
                  group.by = "orig.ident",
                  cols = sample_color,
                  pt.size = 5, 
                  raster = TRUE, 
                  raster.dpi = c(800, 800)) +
  my_theme_void() +
  ggplot2::theme(legend.position = "none")

p2 <- seurat_5 %>%
  Seurat::DimPlot(reduction = "umap.atac",
                  split.by = "condition",
                  group.by = "orig.ident",
                  cols = sample_color,
                  pt.size = 5, 
                  raster = TRUE, 
                  raster.dpi = c(800, 800)) +
  my_theme_void(remove_strip_text = FALSE) +
  ggplot2::theme(title = element_blank(),
                 legend.position = "none")



p3 <- seurat_5 %>%
  Seurat::DimPlot(reduction = "umap.atac",
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
ggplot2::ggsave(p1, filename = here::here("data/dimensional_reduction/atac/liger_peaks_integration_umap.pdf"),
                height = 6,
                width = 6,
                dpi = 500)
ggplot2::ggsave(p2, filename = here::here("data/dimensional_reduction/atac/liger_peaks_integration_umap_diet.pdf"),
                height = 6/3,
                width = 6,
                dpi = 500)
ggplot2::ggsave(p3, filename = here::here("data/dimensional_reduction/atac/liger_peaks_integration_umap_sample.pdf"),
                height = 6/6,
                width = 6,
                dpi = 500)

# WNN ---------------------------------------------------------------------
## we already ran umap for rna in script 4
seurat_5 <- seurat_5 %>%
  Seurat::FindMultiModalNeighbors(
    reduction.list = list("jointly", "liger_embeddings"),
    dims.list = list(1:15, 1:20))


# Run UMAP ----------------------------------------------------------------
seurat_5 <- seurat_5 %>%
  Seurat::RunUMAP(
    nn.name = "weighted.nn",
    reduction.name = "umap.wnn",
    reduction.key = "wnnUMAP_",
    seed.use = 42)

# Plot wnn dim plot - to evaluate integration ----------------------------
p1 <- seurat_5 %>%
  Seurat::DimPlot(reduction = "umap.wnn",
                  group.by = "orig.ident",
                  cols = sample_color,
                  pt.size = 3, 
                  raster = TRUE, 
                  raster.dpi = c(800, 800)) +
  my_theme_void() +
  ggplot2::theme(title = element_blank(),
                 legend.position = "none")

ggplot2::ggsave(p1, filename = here::here("data/dimensional_reduction/wnn/wnn_integration_umap.pdf"),
                height = 6,
                width = 6,
                dpi = 1000)


p2 <- seurat_5 %>%
  Seurat::DimPlot(reduction = "umap.wnn",
                  split.by = "condition",
                  group.by = "orig.ident",
                  cols = sample_color,
                  pt.size = 3, 
                  raster = TRUE, 
                  raster.dpi = c(800, 800)) +
  my_theme_void(remove_strip_text = F) +
  ggplot2::theme(title = element_blank(),
                 legend.position = "none")

ggplot2::ggsave(p2, filename = here::here("data/dimensional_reduction/wnn/wnn_integration_umap_diet.pdf"),
                height = 6/3,
                width = 6,
                dpi = 1000)

p3 <- seurat_5 %>%
  Seurat::DimPlot(reduction = "umap.wnn",
                  split.by = "orig.ident",
                  group.by = "orig.ident",
                  cols = sample_color,
                  pt.size = 5, 
                  raster = TRUE, 
                  raster.dpi = c(800, 800)) +
  my_theme_void(remove_strip_text = F) +
  ggplot2::theme(title = element_blank(),
                 legend.position = "none")

ggplot2::ggsave(p3, filename = here::here("data/dimensional_reduction/wnn/wnn_integration_umap_sample.pdf"),
                height = 6/6,
                width = 6,
                dpi = 1000)


# Find Clusters -----------------------------------------------------------
# using default resolution
seurat_5 <- seurat_5 %>%
  Seurat::FindClusters(
    graph.name = "wsnn",
    res = seq(0,0.8, by = 0.05),
    random.seed = 1000,
    verbose = TRUE)

pdf(file = here::here("data/dimensional_reduction/wnn/wnn_clusttree.pdf"),
    height = 8,
    width = 6)
clustree::clustree(seurat_5, prefix = "wsnn_res.")
dev.off()

# Add gene score from ArchR -----------------------------------------------

# Extract gene score
gene_score <- ArchR::getMatrixFromProject(multi_islets_4, useMatrix='GeneScoreMatrix')

# Gene score matrix
gene_score_matrix <- SummarizedExperiment::assays(gene_score)$GeneScoreMatrix

# colnames
BiocGenerics::colnames(gene_score_matrix) <- BiocGenerics::colnames(gene_score_matrix) %>%
  atac_to_rna_syntax()

# row data
rowdata <- gene_score %>%
  SummarizedExperiment::rowData()
BiocGenerics::rownames(gene_score_matrix) <- rowdata$name

# Add gene score to new assay
seurat_5[["activity"]] <- Seurat::CreateAssayObject(counts = gene_score_matrix)

# Normalize gene scores
Seurat::DefaultAssay(seurat_5) <- "activity"
seurat_5 <- seurat_5 %>%
  Seurat::NormalizeData()


# Save seurat  ------------------------------------------------------------
base::saveRDS(seurat_5, file = here::here("data/seurat_objects/seurat_5.rds"))
