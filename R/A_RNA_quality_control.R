# Description -------------------------------------------------------------
# In this section we will perform empty droplet removal, quality control using the valiDrops  package
# we will use raw data from Starsolo where counts are from exon + intron (genefull)
# and use a contrast where counts are only from exon (gene)
# Additionally we will find doublets (but not remove) using doubletfinder

# Setup -------------------------------------------------------------------
base::source(here::here("R/set_up.R"))
set.seed(1000)

# Create diretories
create_directories(c(here::here("data/quality_control"), 
                     here::here("data/quality_control/rna"),
                     here::here("data/seurat_objects/")))
  
# Load data ---------------------------------------------------------------
# Meta data file
meta <- BiocGenerics::as.data.frame(readxl::read_excel(here::here("data/meta.xlsx"),
                                         sheet = "Sheet1")) %>%
  dplyr::rename_with(snakecase::to_snake_case) %>%
  dplyr::filter(library == "RNA") # keep rows from RNA

# GeneFull mtx ------------------------------------------------------------
# Extract file paths for matrix.mtx, barcodes.tsv and features.tsv
matrix_files <- base::list.files(path = base::Sys.glob(here::here("data/star_output/*/Solo.out/GeneFull_Ex50pAS/raw")), pattern = "matrix.mtx", recursive = TRUE, full.names = TRUE)
barcode_files <- base::list.files(path = base::Sys.glob(here::here("data/star_output/*/Solo.out/GeneFull_Ex50pAS/raw")), pattern = "barcodes.tsv", recursive = TRUE, full.names = TRUE)
feature_files <- base::list.files(path = base::Sys.glob(here::here("data/star_output/*/Solo.out/GeneFull_Ex50pAS/raw")), pattern = "features.tsv", recursive = TRUE, full.names = TRUE)

# Add sample names to each list
base::names(matrix_files) <- purrr::map_chr(matrix_files, ~stringr::str_split(.x, "/|_", simplify = TRUE)[,10])
base::names(barcode_files) <- purrr::map_chr(barcode_files, ~stringr::str_split(.x, "/|_", simplify = TRUE)[,10])
base::names(feature_files) <- purrr::map_chr(feature_files, ~stringr::str_split(.x, "/|_", simplify = TRUE)[,10])

# Check that they are in the same order
base::all.equal(base::names(matrix_files), base::names(barcode_files))
base::all.equal(base::names(matrix_files), base::names(feature_files))

# Read mtx file
mtx_genefull <- purrr::pmap(list(matrix_files, barcode_files, feature_files), function(mtx, cells, features){
  Seurat::ReadMtx(mtx = mtx,
                  cells = cells,
                  features = features)
})

# Replace names
base::names(mtx_genefull) <- stringi::stri_replace_all_regex(
  names(mtx_genefull),
  pattern = meta$seq_id,
  replacement = meta$sample,
  vectorize = FALSE
)

# Gene mtx ----------------------------------------------------------------
# Extract file paths for matrix.mtx, barcodes.tsv and features.tsv
matrix_files <- base::list.files(path = base::Sys.glob(here::here("data/star_output/*/Solo.out/Gene/raw")), pattern = "matrix.mtx", recursive = TRUE, full.names = TRUE)
barcode_files <- base::list.files(path = base::Sys.glob(here::here("data/star_output/*/Solo.out/Gene/raw")), pattern = "barcodes.tsv", recursive = TRUE, full.names = TRUE)
feature_files <- base::list.files(path = base::Sys.glob(here::here("data/star_output/*/Solo.out/Gene/raw")), pattern = "features.tsv", recursive = TRUE, full.names = TRUE)

# Add sample names to each list
base::names(matrix_files) <- purrr::map_chr(matrix_files, ~stringr::str_split(.x, "/|_", simplify = TRUE)[,10])
base::names(barcode_files) <- purrr::map_chr(barcode_files, ~stringr::str_split(.x, "/|_", simplify = TRUE)[,10])
base::names(feature_files) <- purrr::map_chr(feature_files, ~stringr::str_split(.x, "/|_", simplify = TRUE)[,10])

# Check that they are in the same order
base::all.equal(base::names(matrix_files), base::names(barcode_files))
base::all.equal(base::names(matrix_files), base::names(feature_files))

# Read mtx file
mtx_gene <- purrr::pmap(list(matrix_files, barcode_files, feature_files), function(mtx, cells, features){
  Seurat::ReadMtx(mtx = mtx,
                  cells = cells,
                  features = features)
})

# Replace names
base::names(mtx_gene) <- stringi::stri_replace_all_regex(
  names(mtx_gene),
  pattern = meta$seq_id,
  replacement = meta$sample,
  vectorize = FALSE
)


# Change order of files ---------------------------------------------------
mtx_genefull_order <- mtx_genefull[BiocGenerics::match(sample_levels, names(mtx_genefull))]
mtx_gene_order <- mtx_gene[BiocGenerics::match(sample_levels, names(mtx_gene))]

base::all.equal(base::names(mtx_genefull_order), base::names(mtx_gene_order))

# Remove empty droplets ---------------------------------------------------
rank <- purrr::map2(mtx_genefull_order, names(mtx_genefull_order), function(mat, sample) {
  grDevices::pdf(file = base::paste0(here::here("data/quality_control/rna/"), sample, "_rankplot.pdf"))
  threshold <- valiDrops::rank_barcodes(mat)
  dev.off()
  rank.pass <- BiocGenerics::rownames(threshold$ranks[ threshold$ranks$counts >= threshold$lower.threshold, ])
  counts.subset <- mat[, colnames(mat) %in% rank.pass]
  
  # Plot the data (replace with your actual plotting code)
  base::plot(threshold$ranks$counts, main = "Rank Plot", xlab = "Rank", ylab = "Counts")
  
  return(counts.subset) 
})
dev.off()

## Save nonempty barcodes --------------------------------------------------
base::saveRDS(rank, file = here::here("data/quality_control/rna/non_empty_droplet_list.rds"))

# Subset contrast ---------------------------------------------------------
# Make sure sample are in the same order in both count matrixs'
base::all.equal(base::names(rank), base::names(mtx_gene_order))

# Subset list_mtx_gene which will be used as contrast
contrast <- purrr::map2(rank, mtx_gene_order, function(mat, contrast) {
  output <- contrast[, colnames(contrast) %in% colnames(mat)]
  return(output)})

# Create quality metrics --------------------------------------------------
base::all.equal(base::names(rank), base::names(contrast))

quality <- purrr::map2(rank, contrast, function(mat, contrast){
  output <- valiDrops::quality_metrics(counts = mat,
                                       contrast = contrast,
                                       contrast_type = "denominator",
                                       verbose = TRUE)
  return(output)
})


## Save quality stats ------------------------------------------------------
base::saveRDS(quality, file = here::here("data/quality_control/rna/quality_metrics_list.rds"))

# Quality filter ----------------------------------------------------------
# Perform quality filter
quality_filt <- purrr::map2(quality, base::names(quality), function(mat, sample){
  grDevices::pdf(file = base::paste0(here::here("data/quality_control/rna/"), sample, "_quality_filter.pdf"))
  output <- valiDrops::quality_filter(metrics = mat$metrics,
                                      contrast = TRUE)
  dev.off()
  return(output)
})

## Save quality metric -----------------------------------------------------
base::saveRDS(quality_filt, file = here::here("data/quality_control/rna/quality_filter_metrics_list.rds"))

# Subset counts -----------------------------------------------------------
base::all.equal(base::names(rank), base::names(quality_filt))

# keep only high quality nuclei
quality_count <- purrr::map2(rank, quality_filt, function(mat, q_pass){
  output <- mat[, colnames(mat) %in% q_pass$final]
  return(output)
})

## Save counts which passed quality filtering ------------------------------
base::saveRDS(quality_count, file = here::here("data/quality_control/rna/high_quality_barcodes_list.rds"))

# Create seurat object ----------------------------------------------------
seurat_list <- purrr::map2(quality_count, names(quality_count), ~ Seurat::CreateSeuratObject(counts = .x,
                                                                                             project = .y))

# Add meta data -----------------------------------------------------------
base::all.equal(base::names(seurat_list), base::names(quality))

seurat_list <- purrr::map2(seurat_list, quality, function(sobj, q_mat){
  sobj@meta.data <- sobj@meta.data %>%
    BiocGenerics::as.data.frame() %>%
    tibble::rownames_to_column("barcode") %>%
    dplyr::left_join(y = q_mat$metrics, by = "barcode") %>%
    tibble::column_to_rownames("barcode")
  return(sobj)
})

# Rename barcodes ---------------------------------------------------------
# Rename barcodes using object names as prefix
seurat_list <- purrr::imap(seurat_list, ~ Seurat::RenameCells(.x, add.cell.id = .y))

# global filter -----------------------------------------------------------
# following automatic quality control using validrops we set a global mitochondrial fraction threshold of 2.5
# and an exon ratio of 1, as these reads only comes from exons (hence, ambiant rna)
seurat_list <- purrr::map(seurat_list, function(s_obj) {
  # round contrast_fraction to only have 1 decimal, otherwise we wont be able to remove barcodes..
  s_obj@meta.data <- s_obj@meta.data %>%
    BiocGenerics::as.data.frame() %>%
    dplyr::mutate(contrast_fraction_round = round(contrast_fraction, 1))
  
  output <- base::subset(s_obj,
                         subset = mitochondrial_fraction < 0.025 & contrast_fraction_round > 1)
  return(output)
})

purrr::map(seurat_list, function(s_obj){length(colnames(s_obj))}) %>% purrr::as_vector()
             
# Find doublets -----------------------------------------------------------
# here we use DoubletFinder:
# https://github.com/chris-mcginnis-ucsf/DoubletFinder/blob/master/README.md

# Prepare data for doublet identification ---------------------------------
# Normalize counts and scale counts and find top 1000 variable genes
seurat_list <- purrr::map(seurat_list, ~ .x %>%
                     Seurat::NormalizeData() %>%
                     Seurat::FindVariableFeatures(selection.method = "vst", nfeatures = 1000) %>%
                     Seurat::ScaleData() %>%
                     Seurat::RunPCA(assay = "RNA", seed.use = 1000, verbose = FALSE) %>%
                     Seurat::RunUMAP(reduction = "pca", dims = 1:20, verbose = FALSE) %>%
                     Seurat::FindNeighbors() %>%
                     Seurat::FindClusters(resolution = 0.7)
)


# Find doublets -----------------------------------------------------------
## Define expected doublet frequency ----
df_freq <- data.frame(
  "sample" = sample_levels,
  "freq" = c(0.077, 0.077, 0.023, 0.031, 0.077, 0.037, 0.054))

## Find doublets ----
seurat_list_doublets <- purrr::pmap(list(seurat_list, base::names(seurat_list), base::list(df_freq)), process_doublets)

## Plot ----
plot_list <- purrr::map2(seurat_list_doublets, base::names(seurat_list_doublets), function(s_obj, names){
  s_obj@meta.data <- s_obj@meta.data %>%
    BiocGenerics::as.data.frame() %>%
    dplyr::rename_with(~paste0("doubletfinder"), starts_with("DF.classifications"))
  
  plot <- Seurat::DimPlot(s_obj,
                          reduction = 'umap',
                          group.by = "doubletfinder",
                          cols = c("Doublet" = "red", "Singlet" = "grey"),
                          pt.size = 5,
                          raster = TRUE,
                          raster.dpi = c(2000, 2000)) +
    ggplot2::ggtitle(names) +
    ggplot2::labs(x = "UMAP 1",
                  y = "UMAP 2") +
    my_theme_void()
  
  ggplot2::ggsave(plot,
                  filename = paste0(here::here("data/quality_control/rna/"), names, "_DoubletFinder_dimplot.pdf"),
                  height = 6,
                  width = 6)
  
  return(plot)
  })


# Save seurat object ------------------------------------------------------
base::saveRDS(seurat_list_doublets, file = here::here("data/seurat_objects/seurat_object_1.rds"))
