# Description -------------------------------------------------------------
# In this section we will perform quality control on the ATAC data using ArchR.
# And also find doublets (but not remove them)
# Setup -------------------------------------------------------------------
base::source(here::here("R/set_up.R"))
set.seed(1000)

# Create diretories
create_directories(c(here::here("data/quality_control"), 
                     here::here("data/quality_control/atac"),
                     here::here("data/archr_projects")))

# Setting a Genome and GeneAnnotation
ArchR::addArchRGenome("mm10")

# Load data ---------------------------------------------------------------
# Import meta data
meta <- BiocGenerics::as.data.frame(readxl::read_excel(here::here("data/meta.xlsx"),
                                         sheet = "Sheet1")) %>%
  dplyr::rename_with(snakecase::to_snake_case) %>%
  dplyr::filter(library == "ATAC") # keep rows from ATAC

# Preparation for arrow file creation -------------------------------------
# Path to fragment files
inputfiles <-
  base::list.files(base::Sys.glob(here::here("data/cellranger_fragments/*/outs")),
                   pattern = "atac_fragments.tsv.gz",
                   recursive = TRUE,
                   full.names = TRUE)

# Remove .tbi files from inputfiles (these are not used for arrowfile creation)
inputfiles <- inputfiles[!base::grepl(".tbi", inputfiles)]

# Add names to inputfiles, by extracting the ATAC sequencing id.
base::names(inputfiles) <- purrr::map_chr(inputfiles, ~stringr::str_split(.x, "/", simplify = TRUE)[,7])

# Rename the sequence ids of each element to the sample name.
base::names(inputfiles) <- stringi::stri_replace_all_regex(names(inputfiles),
                                                     pattern = meta$seq_id,
                                                     replacement = meta$sample,
                                                     vectorize = FALSE)
# Reorder inputfiles
inputfiles_order <- inputfiles[BiocGenerics::match(sample_levels, names(inputfiles))]

# Create arrow files ------------------------------------------------------
# Using Archr default parameters.
# !These arrow files will be saved in the working directory, but I move them to a folder within the data folder!
ArrowFiles <- inputfiles_order %>%
  ArchR::createArrowFiles(
    sampleNames = names(inputfiles_order),
    minTSS = 4, # don't set this to high, we will increase later
    minFrags = 1000,
    addTileMat = TRUE,
    addGeneScoreMat = TRUE
  )

# Doublet identification --------------------------------------------------
# Archr doublet identification using default parameters

doubScores <- ArchR::addDoubletScores(
  input = ArrowFiles,
  k = 10, # Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
  LSIMethod = 1, 
  threads = parallel::detectCores() - 1
)

# Create the archr project ------------------------------------------------
multi_islets <- ArchR::ArchRProject(
  ArrowFiles = ArrowFiles, 
  threads = parallel::detectCores() - 1,
  outputDirectory =  here::here("data/archr_projects/save_multi_islets"),
  copyArrows = TRUE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)

# save project ------------------------------------------------------------
ArchR::saveArchRProject(ArchRProj = multi_islets,
                 outputDirectory = here::here("data/archr_projects/save_multi_islets"), 
                 threads = parallel::detectCores() - 1,
                 load = TRUE)

# Extract cellcoldata -----------------------------------------------------
CellColData_prefilter <- ArchR::getCellColData(
  multi_islets,
  select = c(
    "Sample",
    "log10(nFrags)",
    "TSSEnrichment",
    "ReadsInTSS",
    "ReadsInPromoter",
    "ReadsInBlacklist",
    "PromoterRatio",
    "PassQC",
    "NucleosomeRatio",
    "nMultiFrags",
    "nMonoFrags",
    "nFrags",
    "nDiFrags",
    "DoubletScore",
    "DoubletEnrichment",
    "BlacklistRatio"
  )
)
# Fragment size distribution plot -----------------------------------------
fragment_size <- ArchR::plotFragmentSizes(ArchRProj = multi_islets,
                                          threads = parallel::detectCores() - 1,
                                          groupBy = "Sample")

# TSS enrichment ----------------------------------------------------------
# Do not run this with 1 threads it will not generate the correct TSS plots
# most of them will be "flat", I run with more than 1 threads and it works.
TSS_enrichment <-
  ArchR::plotTSSEnrichment(ArchRProj = multi_islets,
                           groupBy = "Sample",
                           threads = parallel::detectCores() - 1)

# Quality plots -----------------------------------------------------------
# TSS plot ----------------------------------------------------------------
pdf(file = here::here("data/quality_control/atac/tss_enrichment.pdf"),
    height = 1,
    width = 6)
TSS_enrichment +
  ggplot2::scale_colour_manual(values = sample_color) +
  ggplot2::facet_wrap(~group,
                      scales = "free", ncol = 7) +
  my_theme() +
  ggplot2::theme(legend.position = "none")

dev.off()


# Per nuclei TSS enrichment score -----------------------------------------
pdf(file = here::here("data/quality_control/atac/tss_enrichment_per_nuclei.pdf"),
    height = 1,
    width = 6)

CellColData_prefilter %>%
  BiocGenerics::as.data.frame() %>%
  rstatix::reorder_levels(Sample, order = sample_levels) %>%
  ggplot2::ggplot(aes(x = TSSEnrichment, fill = Sample)) +
  ggplot2::geom_histogram(bins = 50) +
  ggplot2::geom_vline(xintercept = 15, color = "red", linetype = "dashed") +
  ggplot2::scale_fill_manual(values = sample_color) +
  ggplot2::labs(y = "Frequency",
                x = "Transcription start site enrichment score") +
  ggplot2::facet_wrap(~Sample,
                      scales = "free", ncol = 7) +
  my_theme() +
  ggplot2::theme(legend.position = "none")

dev.off()

# Fragment size distrubution plot -----------------------------------------
pdf(file = here::here("data/quality_control/atac/fragment_size_distribution.pdf"),
    height = 1,
    width = 6)

fragment_size +
  ggplot2::scale_colour_manual(values = sample_color) +
  ggplot2::facet_wrap(~group,
                      scales = "free", ncol = 7) +
  my_theme() +
  ggplot2::theme(legend.position = "none")

dev.off()


# Number of unique fragments ----------------------------------------------
pdf(file = here::here("data/quality_control/atac/number_of_unique_fragments.pdf"),
    height = 1,
    width = 6)

CellColData_prefilter %>%
  as.data.frame() %>%
  ggplot2::ggplot(aes(x = nFrags, fill = Sample)) +
  ggplot2::geom_histogram(bins = 50) +
  ggplot2::geom_vline(xintercept = 2500, color = "red", linetype = "dashed") +
  ggplot2::scale_x_log10(labels = scales::comma) +
  ggplot2::scale_fill_manual(values = sample_color) +
  ggplot2::labs(y = "Frequency",
                x = "Number of unique fragments (log10 scale)") +
  ggplot2::facet_wrap(~Sample,
                      scales = "free", ncol = 7) +
  my_theme() +
  ggplot2::theme(legend.position = "none")

dev.off()


# Black list ratio --------------------------------------------------------
pdf(file = here::here("data/quality_control/atac/black_list_ratio.pdf"),
    height = 1,
    width = 6)

CellColData_prefilter %>%
  as.data.frame() %>%
  ggplot2::ggplot(aes(x = BlacklistRatio, fill = Sample)) +
  ggplot2::geom_histogram(bins = 100) +
  ggplot2::geom_vline(xintercept = 0.05, color = "red", linetype = "dashed") +
  ggplot2::scale_fill_manual(values = sample_color) +
  ggplot2::labs(y = "Frequency",
                x = "Blacklist ratio") +
  ggplot2::facet_wrap(~Sample,
                      scales = "free", ncol = 7) +
  my_theme() +
  ggplot2::theme(legend.position = "none")

dev.off()


# TSS vs log10 fragments --------------------------------------------------
pdf(file = here::here("data/quality_control/atac/tss_vs_log10nFrags.pdf"),
    height = 1,
    width = 6)

# Split CellColData_postfilter by sample
df <- CellColData_prefilter %>%
  BiocGenerics::as.data.frame() %>%
  dplyr::select(log10.nFrags., TSSEnrichment, Sample) %>%
  dplyr::mutate(split_by = Sample) %>%
  collapse::rsplit(~ split_by)


p7 <- Map(ggPoint_plot,
          .df = df,
          .title = names(df))

p7[2:7] <- purrr::map(p7[2:7], ~ .x + ggplot2::theme(axis.title.y = element_blank()))
patchwork::wrap_plots(p7, ncol = 7) + patchwork::plot_layout(guides = "collect") & 
  ggplot2::guides(colour = ggplot2::guide_colourbar(barwidth = 0.5, barheight = 6/2),
                  size = "none") 
dev.off()

# Filter Archr files ------------------------------------------------------
multi_islets_1 <- multi_islets[multi_islets$nFrags >= 2500 &
                                 multi_islets$TSSEnrichment >= 15 &
                                 multi_islets$BlacklistRatio <= 0.05]
# Find doublets -----------------------------------------------------------
archr_doublet_filter <- ArchR::filterDoublets(multi_islets_1)

# Get doublet barcodes
doublets_ATAC <- BiocGenerics::setdiff(BiocGenerics::as.vector(multi_islets_1$cellNames),
                                       BiocGenerics::as.vector(archr_doublet_filter$cellNames))

# Save Archr project ------------------------------------------------------
saveArchRProject(ArchRProj = multi_islets_1,
                 outputDirectory = here::here("data/archr_projects/save_multi_islets_1"),
                 threads = parallel::detectCores() - 1,
                 load = FALSE)
# Save doublets -----------------------------------------------------------
base::saveRDS(doublets_ATAC, file = here::here("data/quality_control/atac/doublet_barcodes_vector.rds"))

