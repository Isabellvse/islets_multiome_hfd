# Description -------------------------------------------------------------
# for endocrine cells
# Create Anndata object
# and extract fragment counts
# Setup -------------------------------------------------------------------
base::source(here::here("R/set_up.R"))
set.seed(1000)
create_directories(c(here::here("data/scenicplus/beta/files"),
                     here::here("data/scenicplus/beta/figures")))

# Load --------------------------------------------------------------------
seurat_5 <- base::readRDS(here::here("data/seurat_objects/seurat_5.rds"))

# BETA CELLS -------------------------------------------------------------
# Create Ann data -----

# subset to only keep Beta-cells
Seurat::Idents(seurat_5) <- "manual_anno"
s_Beta <- base::subset(x = seurat_5, idents = "Beta")

Seurat::DefaultAssay(s_Beta) <- "RNA"

seurat_rna <- Seurat::DietSeurat(s_Beta,
                                 assays = "RNA",
                                 dimreducs = c("jointly", "jointly_pca", "umap.wnn"))

# change data slot with count data
seurat_rna@assays$RNA@data <- seurat_rna@assays$RNA@counts
base::all.equal(seurat_rna@assays$RNA@counts, seurat_rna@assays$RNA@data)

# meta data
seurat_rna@meta.data <- seurat_rna@meta.data %>%
  BiocGenerics::as.data.frame() %>%
  dplyr::rename(sample_id = orig.ident) %>%
  dplyr::select(sample_id, manual_anno, condition) %>%
  dplyr::mutate(dplyr::across(dplyr::everything(), base::as.character))

# save
SeuratDisk::SaveH5Seurat(seurat_rna, filename =here::here("data/scenicplus/beta/files/rna.h5Seurat"))

# convert to ann data object
SeuratDisk::Convert(source = here::here("data/scenicplus/beta/files/rna.h5Seurat"),
                    dest = "H5ad",
                    assay = "RNA")

# create cistopic object
meta_atac <- seurat_rna@meta.data

# save
utils::write.table(meta_atac, file=here::here("data/scenicplus/beta/files/atac_meta.csv"), quote=FALSE, sep='\t')

# get atac data
matrix <- Seurat::GetAssayData(s_Beta, assay = "peaks", slot = "counts")

# change rownames
BiocGenerics::rownames(matrix) <- stringi::stri_replace_first_regex(rownames(matrix), "-", ":")

# save as csv
scrattch.io::write_dgCMatrix_csv(matrix, here::here("data/scenicplus/beta/files/counts_atac.csv"), col1_name = "row_names",
                                 chunk_size = 100000)
