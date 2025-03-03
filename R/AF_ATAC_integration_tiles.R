# Description -------------------------------------------------------------
# Dimensional reduction and liger integration of tiles
# Setup -------------------------------------------------------------------
base::source(here::here("R/set_up.R"))
set.seed(1000)

# Setting a Genome and GeneAnnotation
ArchR::addArchRGenome("mm10")

create_directories(c(here::here("data/dimensional_reduction"), 
                     here::here("data/dimensional_reduction/atac")))
# # Load --------------------------------------------------------------------
multi_islets_3 <- ArchR::loadArchRProject(path = here::here("data/archr_projects/save_multi_islets_2/"))

ArchR::saveArchRProject(ArchRProj = multi_islets_3,
                 outputDirectory = here::here("data/archr_projects/save_multi_islets_3"),
                 threads = parallel::detectCores() - 1,
                 load = FALSE)

multi_islets_3 <- ArchR::loadArchRProject(path = here::here("data/archr_projects/save_multi_islets_3/"))

# Create 5000 bp tile matrix ----------------------------------------------
# Add 5000 bp tile matrix
multi_islets_3 <- ArchR::addTileMatrix(
  input = multi_islets_3,
  tileSize = 5000,
  threads = parallel::detectCores() - 1,
  force = TRUE)

# Dimensional reduction ---------------------------------------------------

# Iterative LSI
multi_islets_3 <-
  ArchR::addIterativeLSI(
    ArchRProj = multi_islets_3,
    useMatrix = "TileMatrix",
    name = "IterativeLSI",
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

# save archer project
ArchR::saveArchRProject(ArchRProj = multi_islets_3,
                 outputDirectory = here::here("data/archr_projects/save_multi_islets_3"),
                 threads = parallel::detectCores() - 1,
                 load = FALSE)

# Get variable features ---------------------------------------------------
multi_islets_3 <- ArchR::loadArchRProject(path = here::here("data/archr_projects/save_multi_islets_3/"))

# Load variable features used in LSI, make each fature the syntax "chr-start-end" and convert to vector
variable_features <- multi_islets_3@reducedDims$IterativeLSI@listData[["LSIFeatures"]] %>%
  BiocGenerics::as.data.frame() %>%
  dplyr::mutate(end = (start + 5000) - 1, # Create column of end position of bin.
                feature = paste0(seqnames, "-", start, "-", end)) %>% # create row data chr-start-end
  dplyr::pull(feature) # create vector

# Extract tile matrix
tiles <- ArchR::getMatrixFromProject(
  ArchRProj = multi_islets_3,
  useMatrix = "TileMatrix",
  useSeqnames = NULL,
  verbose = TRUE,
  binarize = TRUE,
  threads = parallel::detectCores() - 1,
  logFile = createLogFile("getMatrixFromProject")
)

# Split tiles matrix according to sample

tiles_list <- list()
for (i in sample_levels) {
  # Select columns of tiles corresponding to the current sample level
  new_element <- tiles[, tiles$Sample == i]
  # Add the new element to the list with the sample level as the name
  tiles_list[[i]] <- new_element
}


# Transform tile matrix into dataframe and also create a feature column with the syntax "chr-start-end"
# Extract matrix
tiles_matrix <- purrr::map(tiles_list, function(.x){
  # extract data
  matrix <- SummarizedExperiment::assays(.x)$TileMatrix
  # row names
  tiles_row <- SummarizedExperiment::rowData(.x)
  tiles_row <- tiles_row %>%
    BiocGenerics::as.data.frame() %>%
    dplyr::mutate(end = (start + 5000) - 1, # Create column of end position of bin.
                  feature = paste0(seqnames, "-", start, "-", end)) # create row data chr-start-end
  BiocGenerics::rownames(matrix) <- tiles_row$feature
  return(matrix)})

# Create liger object
# Using code form: https://github.com/theislab/scib-pipeline/blob/main/scripts/integration/integration.R
liger_tiles <- rliger::createLiger(tiles_matrix,
                                   remove.missing = FALSE)

liger_tiles@norm.data <- liger_tiles@raw.data

# Assign highly variable features
liger_tiles@var.genes <- variable_features

# remove unused data to increase space
rm("multi_islets_3", "tiles", "tiles_matrix", "samples", "tiles_list")

# free unsued memory
gc()

# scale data
liger_tiles <- rliger::scaleNotCenter(liger_tiles,
                                      remove.missing = F,
                                      verbose = TRUE)
liger_tiles <- liger_tiles %>%
  rliger::optimizeALS(k = 20,
                      lamda = 5,
                      thresh = 5e-5,
                      nrep = 3)

# quantile align snf
liger_tiles <- liger_tiles %>%
  rliger::quantile_norm(resolution = 0.4,
                        small.clust.thresh = 20)
# save liger object
base::saveRDS(liger_tiles, file = here::here("data/dimensional_reduction/atac/liger_integration_5000_bp_tiles.rds"))
