# Description -------------------------------------------------------------
# chromvar analysis using motifs from scenic+
# Setup -------------------------------------------------------------------
base::source(here::here("R/set_up.R"))
set.seed(1000)

# Setting a Genome and GeneAnnotation
ArchR::addArchRGenome("mm10")

# Load data ---------------------------------------------------------------
# motifs
pwmlist <- base::readRDS(here::here("data/scenicplus/scenic_PWMLIST_mouse.rds"))

# archr object
multi_islets_4 <- ArchR::loadArchRProject(path = here::here("data/archr_projects/save_multi_islets_4/"))
multi_islets_5 <- ArchR::loadArchRProject(path = here::here("data/archr_projects/save_multi_islets_5/"))

# seurat 
seurat_5 <- base::readRDS(here::here("data/seurat_objects/seurat_5.rds"))

# Subset Archr project ----------------------------------------------------
# only contain beta-cells
cells_keep <- seurat_5@meta.data %>%
  dplyr::filter(manual_anno == "Beta") %>%
  rownames() %>%
  rna_to_atac_syntax()

multi_islets_5 <- ArchR::subsetArchRProject(
  ArchRProj = multi_islets_4,
  cells = cells_keep,
  outputDirectory = here::here("data/archr_projects/save_multi_islets_5"),
  dropCells = TRUE,
  logFile = NULL,
  threads = parallel::detectCores() - 1,
  force =  TRUE)

# Chromvar ----------------------------------------------------------------
### Add motif annotation to ArchR project ----
multi_islets_5 <- ArchR::addMotifAnnotations(
  ArchRProj = multi_islets_5,
  motifSet = "NULL",
  name = "Motif",
  motifPWMs = pwmlist,
  force = TRUE)

# Add set of background peaks
multi_islets_5 <- ArchR::addBgdPeaks(multi_islets_5)

# Calculate deviations
multi_islets_5<- ArchR::addDeviationsMatrix(
  ArchRProj = multi_islets_5,
  peakAnnotation = "Motif",
  force = TRUE,
  threads = parallel::detectCores() - 1)

# multi_islets_beta_scenic
ArchR::saveArchRProject(ArchRProj = multi_islets_5,
                        outputDirectory = here::here("data/archr_projects/save_multi_islets_5"),
                        load = TRUE,
                        threads = parallel::detectCores() - 1)

# Extract results ---------------------------------------------------------
# Extract chromVAR results from ArchR project
chromvar_res <- ArchR::getMatrixFromProject(
  ArchRProj = multi_islets_5,
  useMatrix = "MotifMatrix",
  threads = parallel::detectCores() - 1)

## Deviations ----
# Extract deviations
dev_motifs <- as.data.frame(assay(chromvar_res, "deviations"))

# Change barcodes to RNA syntax
colnames(dev_motifs) <- atac_to_rna_syntax(colnames(dev_motifs))


# save --------------------------------------------------------------------
base::saveRDS(dev_motifs, here::here("data/scenicplus/beta/results/chromvar_deviations_motifs_beta.rds"))


# add motif annotation manually -------------------------------------------
motif_summary <- readRDS(here::here("data/archr_projects/save_multi_islets_5/Annotations/Motif-In-Peaks-Summary.rds"))
motif_matches_in_peaks <- readRDS(here::here("data/archr_projects/save_multi_islets_5/Annotations/Motif-Matches-In-Peaks.rds"))
motif_positions_in_peaks <- readRDS(here::here("data/archr_projects/save_multi_islets_5/Annotations/Motif-Positions-In-Peaks.rds"))


multi_islets_5@peakAnnotation[["Motif"]] <- list(
  Name = "Motif",
  motifs = motif_summary$motifList,
  motifSummary = motif_summary$motifSummary,
  Positions = here::here("data/archr_projects/save_multi_islets_5/Annotations/Motif-Positions-In-Peaks.rds"),
  Matches = here::here("data/archr_projects/save_multi_islets_5/Annotations/Motif-Matches-In-Peaks.rds")
)
