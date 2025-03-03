# Description -------------------------------------------------------------
# Liger integration peaks
# Setup -------------------------------------------------------------------
base::source(here::here("R/set_up.R"))
set.seed(1000)

# Setting a Genome and GeneAnnotation
ArchR::addArchRGenome("mm10")

# Load --------------------------------------------------------------------
multi_islets_4 <- ArchR::loadArchRProject(path = here::here("data/archr_projects/save_multi_islets_3/"))

ArchR::saveArchRProject(ArchRProj = multi_islets_4,
                        outputDirectory = here::here("data/archr_projects/save_multi_islets_4"),
                        load = FALSE,
                        threads = parallel::detectCores() - 1)

multi_islets_4 <- ArchR::loadArchRProject(path = here::here("data/archr_projects/save_multi_islets_4/"))

# # Iterative LSI -----------------------------------------------------------

# Integrate Archr peaks using harmony(5000 bp) - default parameters expet for number of iterations
multi_islets_4 <-
  ArchR::addIterativeLSI(
    ArchRProj = multi_islets_4,
    useMatrix = "PeakMatrix",
    name = "IterativeLSI_peaks",
    iterations = 4,
    threads = parallel::detectCores() - 1,
    clusterParams = list( #See Seurat::FindClusters
      resolution = c(0.2),
      sampleCells = 10000,
      n.start = 10
    ),
    varFeatures = 150000,
    dimsToUse = 1:30
  )

#save
ArchR::saveArchRProject(ArchRProj = multi_islets_4,
                        outputDirectory = here::here("data/archr_projects/save_multi_islets_4"),
                        load = TRUE,
                        threads = parallel::detectCores() - 1)

# Load variable features used in LSI, make each fature the syntax "chr-start-end" and convert to vector
variable_features <- multi_islets_4@reducedDims$IterativeLSI_peaks@listData[["LSIFeatures"]] %>%
  BiocGenerics::as.data.frame() %>%
  dplyr::mutate(feature = base::paste0(seqnames, "-", start, "-", end)) %>% # create row data chr-start-end
  dplyr::pull(feature) # create vector

# Extract peak matrix
peaks <- ArchR::getMatrixFromProject(
  ArchRProj = multi_islets_4,
  useMatrix = "PeakMatrix",
  useSeqnames = NULL,
  verbose = TRUE,
  binarize = TRUE,
  threads = parallel::detectCores() - 1,
  logFile = createLogFile("getMatrixFromProject")
)

# add rownames to peaks
BiocGenerics::rownames(peaks) <- base::paste0(base::as.character(GenomeInfoDb::seqnames(multi_islets_4@peakSet)), '-', base::as.character(BiocGenerics::start(multi_islets_4@peakSet)), '-', base::as.character(BiocGenerics::end(multi_islets_4@peakSet)))

# Split peaks matrix according to sample
peaks_list <- list()
for (i in sample_levels) {
  # Select columns of tiles corresponding to the current sample level
  new_element <- peaks[, peaks$Sample == i]
  # Add the new element to the list with the sample level as the name
  peaks_list[[i]] <- new_element
}


# Extract peak matrix
peaks_matrix <- purrr::map(peaks_list, function(.x){
  # extract data
  matrix <- SummarizedExperiment::assays(.x)$PeakMatrix
  return(matrix)})

# Create liger object
# Using code form: https://github.com/theislab/scib-pipeline/blob/main/scripts/integration/integration.R
liger_peaks <- rliger::createLiger(peaks_matrix,
                                   remove.missing = FALSE)

liger_peaks@norm.data <- liger_peaks@raw.data

# save liger object
base::saveRDS(liger_peaks, file = here::here("data/dimensional_reduction/atac/liger_integration_peaks_pre.rds"))

# remove unused files to get more space -----------------------------------
rm(peaks_list, peaks, peaks_matrix, new_element)

# free unused memory
gc()

# liger integration -------------------------------------------------------

# Assign highly variable features
liger_peaks@var.genes <- variable_features

# scale data
liger_peaks <- rliger::scaleNotCenter(liger_peaks,
                                      remove.missing = F,
                                      verbose = TRUE)
liger_peaks <- liger_peaks %>%
  optimizeALS(k = 20,
              lamda = 5,
              thresh = 5e-5,
              nrep = 3)

# quantile align snf
liger_peaks <- liger_peaks %>%
  rliger::quantile_norm(resolution = 0.4,
                        small.clust.thresh = 20)
# save liger object
base::saveRDS(liger_peaks, file = here::here("data/dimensional_reduction/atac/liger_integration_peaks.rds"))
