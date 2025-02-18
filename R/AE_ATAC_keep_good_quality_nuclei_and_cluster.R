# Description -------------------------------------------------------------
# In this section we will keep barcodes which has passed quality control in both
# RNA and ATAC(identified in script 3), as well as removing doublets (identified in script 4)

# Setup -------------------------------------------------------------------
base::source(here::here("R/set_up.R"))
set.seed(1000)

# Setting a Genome and GeneAnnotation
ArchR::addArchRGenome("mm10")

# Load --------------------------------------------------------------------
multi_islets_2 <- ArchR::loadArchRProject(path = here::here("data/archr_projects/save_multi_islets_1/"))

ArchR::saveArchRProject(ArchRProj = multi_islets_2,
                 outputDirectory = here::here("data/archr_projects/save_multi_islets_2"),
                 threads = parallel::detectCores() - 1,
                 load = FALSE)

multi_islets_2 <- ArchR::loadArchRProject(path = here::here("data/archr_projects/save_multi_islets_2/"))

# Barcodes to keep
bar_keep <- base::readRDS(here::here("data/quality_control/rna_atac_barcode_keep.rds"))

# Doublets - no integration
db <- base::readRDS(here::here("data/quality_control/rna/ALL_polyhormone_cells_remove.rds"))

# Meta data
meta <- BiocGenerics::as.data.frame(readxl::read_excel(here::here("data/meta.xlsx"),
                                         sheet = "Sheet1")) %>%
  dplyr::rename_with(snakecase::to_snake_case) %>%
  dplyr::filter(library == "ATAC") %>% # keep rows from ATA
  dplyr::mutate(condition = stringr::str_c(diet,
                                           time,
                                           sep = "_",
                                           collapse = NULL)) %>%
  dplyr::mutate(condition = stringr::str_replace(condition, "[L][F][D]_[0]", "LFD")) # add condition, will be used later



# subset archr project ----------------------------------------------------
# remove db cells from bar_keep
bar_keep <- bar_keep[!bar_keep %in% db]

multi_islets_2 <- subsetArchRProject(
  ArchRProj = multi_islets_2,
  cells = rna_to_atac_syntax(bar_keep),
  outputDirectory = "data/archr_projects/save_multi_islets_2",
  dropCells = TRUE,
  logFile = NULL,
  threads = parallel::detectCores() - 1,
  force =  TRUE
)

# save archr project ------------------------------------------------------
saveArchRProject(ArchRProj = multi_islets_2,
                 outputDirectory = here::here("data/archr_projects/save_multi_islets_2"),
                 threads = parallel::detectCores() - 1,
                 load = FALSE)


