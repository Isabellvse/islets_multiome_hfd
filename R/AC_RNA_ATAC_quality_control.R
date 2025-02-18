# Description -------------------------------------------------------------
# In this section we will keep barcodes which has passed quality control in both
# RNA and ATAC. So both the Archr and seurat object will be subset to contain these cells.
# and we will keep barcodes which are high quality in both assays

# Setup -------------------------------------------------------------------
base::source(here::here("R/set_up.R"))
set.seed(1000)

# Create diretories
create_directories(c(here::here("data/quality_control"), 
                     here::here("data/quality_control/atac"),
                     here::here("data/archr_projects")))
set.seed(1000)

# Setting a Genome and GeneAnnotation
ArchR::addArchRGenome("mm10")

# Load --------------------------------------------------------------------
# Archr project
multi_islets_1 <- ArchR::loadArchRProject(path = here::here("data/archr_projects/save_multi_islets_1/"))

# Seurat object list
seurat_list <- base::readRDS(file = here::here("data/seurat_objects/seurat_object_1.rds"))

# doublet_atac
db_atac <- base::readRDS(here::here("data/quality_control/atac/doublet_barcodes_vector.rds"))

# Get doublets identified in both assays ----------------------------------
# rna
db_rna <- purrr::map(seurat_list, function(s_obj){
  output <- s_obj@meta.data %>%
    as.data.frame() %>%
    dplyr::rename_with(~paste0("doubletfinder"), starts_with("DF.classifications")) %>%
    dplyr::filter(doubletfinder == "Doublet") %>%
    BiocGenerics::rownames()
  return(output)
}) %>% BiocGenerics::unlist() %>% base::unname()

# Common doublets between rna and atac
db_common <- BiocGenerics::intersect(db_rna, atac_to_rna_syntax(db_atac))

# Get barcodes that are good quality in both assays -----------------------
bar_rna <- purrr::map(seurat_list, function(s_obj){
  output <- SeuratObject::Cells(s_obj)
  return(output)
}) %>% unlist()

bar_atac <- BiocGenerics::rownames(multi_islets_1@cellColData)

# Common good quality nuclei between both assays
bar_common <- BiocGenerics::intersect(bar_rna, atac_to_rna_syntax(bar_atac))

# Remove doublets
bar_keep <- bar_common[!bar_common %in% db_common]

# save barcodes to keep ---------------------------------------------------
base::saveRDS(db_rna, file = here::here("data/quality_control/rna/rna_doublet_barcodes.rds"))
base::saveRDS(bar_keep, file = here::here("data/quality_control/rna_atac_barcode_keep.rds"))
base::saveRDS(db_common, file = here::here("data/quality_control/rna_atac_doublets.rds"))

