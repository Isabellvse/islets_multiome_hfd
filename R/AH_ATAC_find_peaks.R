# Description -------------------------------------------------------------
# peak detection
# Setup -------------------------------------------------------------------
base::source(here::here("R/set_up.R"))
set.seed(1000)

# Setting a Genome and GeneAnnotation
ArchR::addArchRGenome("mm10")

# Load --------------------------------------------------------------------
multi_islets_3 <- ArchR::loadArchRProject(path = here::here("data/archr_projects/save_multi_islets_3/"))

# # Meta data
# meta <- BiocGenerics::as.data.frame(readxl::read_excel(here::here("data/meta.xlsx"),
#                                          sheet = "Sheet1")) %>%
#   dplyr::rename_with(snakecase::to_snake_case) %>%
#   dplyr::filter(library == "ATAC") %>% # keep rows from ATA
#   dplyr::mutate(condition = stringr::str_c(diet,
#                                            time,
#                                            sep = "_",
#                                            collapse = NULL)) %>%
#   dplyr::mutate(condition = stringr::str_replace(condition, "[L][F][D]_[0]", "LFD")) # add condition, will be used later
# 
# # Prepare data ------------------------------------------------------------
# 
# # # How many cells are in each cluster from each sample ---------------------
# #
# # # Confusion matrix
# cM <- ArchR::confusionMatrix(base::paste0(multi_islets_3$cluster_anno), base::paste0(multi_islets_3$Sample))
# cM_df <- cM %>% BiocGenerics::as.data.frame() %>% dplyr::relocate(dplyr::all_of(sample_levels))
# 
# # Prepare data
# # Add info on experimental condition (diet + time on diet)
# multi_islets_3 <- ArchR::addCellColData(
#   ArchRProj = multi_islets_3,
#   data = stringi::stri_replace_all_regex(
#     multi_islets_3$Sample,
#     pattern = meta$sample,
#     replacement = meta$condition,
#     vectorize = FALSE),
#   name = "condition",
#   ArchR::getCellNames(multi_islets_3),
#   force = TRUE
# )
# 
# # Add info on experimental condition matched with Cluster identity to project
# multi_islets_3 <- ArchR::addCellColData(
#   ArchRProj = multi_islets_3,
#   data = stringr::str_c(
#     multi_islets_3$condition,
#     multi_islets_3$cluster_anno,
#     sep = "_",
#     collapse = NULL
#   ),
#   name = "ClusterByCond",
#   cells = ArchR::getCellNames(multi_islets_3),
#   force = TRUE
# )
# 
# # Create pseudo replicates ------------------------------------------------
# multi_islets_3 <- ArchR::addGroupCoverages(
#   ArchRProj = multi_islets_3,
#   groupBy = "ClusterByCond",
#   useLabels = TRUE,
#   minCells = 20,
#   maxCells = 5213,
#   maxFragments = 25 * 10^6,
#   minReplicates = 2,
#   maxReplicates = 3,
#   sampleRatio = 0.8,
#   kmerLength = 6,
#   returnGroups = FALSE,
#   parallelParam = NULL,
#   force = TRUE,
#   verbose = TRUE,
#   logFile = createLogFile("addGroupCoverages"))
# 
# ArchR::saveArchRProject(ArchRProj = multi_islets_3,
#                         outputDirectory = here::here("data/archr_projects/save_multi_islets_3"),
#                         load = TRUE,
#                         threads = parallel::detectCores() - 1)

# Find peaks --------------------------------------------------------------
# Find peaks
multi_islets_3 <- ArchR::addReproduciblePeakSet(
  ArchRProj = multi_islets_3,
  groupBy = "ClusterByCond",
  threads = parallel::detectCores() - 1,
  pathToMacs2 = here::here("macs2_env/bin/macs2"))

# Compute peak counts -----------------------------------------------------

# compute counts for each peak per cell in the provided ArchRProject
multi_islets_3 <- ArchR::addPeakMatrix(
  ArchRProj = multi_islets_3,
  ceiling = 4,
  binarize = FALSE,
  verbose = TRUE,
  threads = parallel::detectCores() - 1,
  force = TRUE,
  logFile = createLogFile("addPeakMatrix")
)


# Save project ------------------------------------------------------------
ArchR::saveArchRProject(ArchRProj = multi_islets_3,
                        outputDirectory = here::here("data/archr_projects/save_multi_islets_3"),
                        threads = parallel::detectCores() - 1,
                        load = TRUE)


# FRIP score --------------------------------------------------------------
pdf(file = here::here("data/quality_control/atac/FRIP.pdf"),
    height = 6,
    width = 6, )

multi_islets_3@cellColData %>%
  as.data.frame() %>%
  rstatix::reorder_levels(Sample, order = sample_levels) %>%
  ggplot2::ggplot(aes(x = FRIP, fill = Sample)) +
  ggplot2::geom_histogram(bins = 50) +
  ggplot2::geom_vline(xintercept = 0.2, color = "red") +
  ggplot2::scale_fill_manual(values = condition_color) +
  ggplot2::labs(y = "Frequency",
                x = "Fraction of reads in peaks (FRIP) score") +
  ggplot2::facet_wrap(~Sample,
                      scales = "free") +
  my_theme() +
  ggplot2::theme(legend.position = "none")

multi_islets_3@cellColData %>%
  as.data.frame() %>%
  rstatix::reorder_levels(Sample, order = sample_levels) %>%
  ggplot2::ggplot(aes(y = FRIP, x = Sample, fill = Sample)) +
  ggplot2::geom_violin() +
  ggplot2::geom_boxplot(outlier.shape = NA, width = 0.5) +
  ggplot2::geom_hline(yintercept = 0.2, color = "red") +
  ggplot2::scale_fill_manual(values = sample_color) +
  ggplot2::labs(x = "Frequency",
                y = "Fraction of reads in peaks (FRIP) score") +
  my_theme() +
  ggplot2::theme(legend.position = "none")

dev.off()

