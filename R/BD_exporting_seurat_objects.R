# Description -------------------------------------------------------------
# Exporting seurat objects

# Setup -------------------------------------------------------------------
base::source(here::here("R/set_up.R"))
set.seed(1000)

# load --------------------------------------------------------------------
seurat_5 <- readRDS(here::here("data/seurat_objects/seurat_5.rds"))
seurat_beta <- readRDS(here::here("data/inflammatory_cluster_2/files/seurat_beta.rds"))

# save --------------------------------------------------------------------
SeuratDisk::SaveH5Seurat(seurat_5, filename ="data/export/islets_multiome.h5Seurat")
SeuratDisk::SaveH5Seurat(seurat_beta, filename ="data/export/islets_multiome_beta_cells.h5Seurat")
