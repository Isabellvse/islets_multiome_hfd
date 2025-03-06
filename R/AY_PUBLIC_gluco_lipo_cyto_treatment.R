# Description -------------------------------------------------------------
# data from this study # GSE218316
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE218316

# See how genes which are differentially expression between LFD, HFD1 and HFD3 are expressed in 
# human beta-cells treated with stressors
# Setup -------------------------------------------------------------------
base::source(here::here("R/set_up.R"))
set.seed(1000)

# set up paths ------------------------------------------------------------
# Create variable with main directory path to save results
path_main <- here::here("data/public_data/")

create_directories(c(here::here("data/public_data/"),
                     here::here("data/public_data/glu_lip_cyt_human/GSE218316/"),
                     here::here("data/public_data/glu_lip_cyt_human/GSE218316/download/"),
                     here::here("data/public_data/glu_lip_cyt_human/GSE218316/processed/")))
# download ----------------------------------------------------------------
# list of files to download

# Filtered barcodes and meta dataframe
files_dw <- list("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE218316&format=file&file=GSE218316%5Ffiltered%5Fcounts%5Fbarcodes%2Etsv%2Egz",
                 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE218316&format=file&file=GSE218316%5Ffiltered%5Fcounts%5Fdata%2Emtx%2Egz",
                 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE218316&format=file&file=GSE218316%5Ffiltered%5Fcounts%5Fgene%2Etsv%2Egz",
                 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE218316&format=file&file=GSE218316%5Fmetadata%5Fpercell%2Ecsv%2Egz")
# it would not let me download it, so I had to upload it to my pc and uploaded it
## unzip all files
# unzip all gz files
paths <- base::list.files(path = output_dir[["download"]], pattern = ".gz", recursive = TRUE, full.names = TRUE)
plyr::ldply(.data = paths, .fun = R.utils::gunzip)

#untar files
# unzip all gz files
paths <- base::list.files(path = output_dir[["download"]], pattern = ".tar", recursive = TRUE, full.names = TRUE)
plyr::ldply(c(.data = paths, .fun = untar))

# untar files
utils::untar(tarfile = here::here("data/public_data/glu_lip_cyt_human/GSE218316/download/GSE218316_RAW.tar"),
             exdir = here::here("data/public_data/glu_lip_cyt_human/GSE218316/download/raw"))
# Load --------------------------------------------------------------------
## meta data
GSE218316_metadata_percell <- read_csv("data/public_data/glu_lip_cyt_human/GSE218316/download/GSE218316_metadata_percell.csv")

## gene clusters 
gene_clusters <- base::readRDS(here::here("data/deseq2/rna_clusterbycond/files/clustering_new.rds")) %>% 
  purrr::pluck("Beta") %>% 
  purrr::pluck("cluster") %>% 
  BiocGenerics::as.data.frame() %>% 
  tibble::rownames_to_column("gene") %>% 
  dplyr::rename(cluster = ".") %>% 
  base::split(factor(.$cluster)) %>% 
  purrr::modify_depth(1, ~dplyr::pull(.,"gene"))
  
  
# Preprocess --------------------------------------------------------------
# only keep conditions you are interested in
test <- GSE218316_metadata_percell %>%
  dplyr::filter(Condition %in% c("Glucose+Palmitate", "IL1b+IFNg", "IL1b", "Glucose", "Untreated", "Palmitate", "IFNg", "IFNa")) %>%
  dplyr::rename(barcode = ...1)


# find out which donors have B something infront of them
prefix_donor <- test %>%
  dplyr::mutate(prefix = stringr::str_extract(barcode, "^.{2}"),
                sample = paste0(Donor, "_", Timepoint)) %>%
  dplyr::group_by(Donor, prefix, sample) %>%
  dplyr::summarise()

# Create seurat object ----------------------------------------------------
# Extract file paths for matrix.mtx, barcodes.tsv and features.tsv
# keep only the files from GSE183010_meta
matrix_files <- list.files(path = paste0(output_dir[["download"]], "/raw/"),
                           pattern = "matrix.mtx", recursive = TRUE, full.names = TRUE)
barcode_files <- list.files(path = paste0(output_dir[["download"]], "/raw/"),
                            pattern = "barcodes.tsv", recursive = TRUE, full.names = TRUE)
feature_files <- list.files(path = paste0(output_dir[["download"]], "/raw/"),
                            pattern = "features.tsv", recursive = TRUE, full.names = TRUE)

### Load mtx files into R ----

# get samples
sample <- stringr::str_extract(matrix_files, "[D][:digit:]{1}[_][:digit:]{2}[h]")

# Create mtx file, using gene ids (genefull - intron + exon)
mtx <- list()

for(i in 1:length(sample)) {
  new_element <- Seurat::ReadMtx(
    feature.column = 2,
    mtx = matrix_files[grep(sample[i], matrix_files)],
    cells =  barcode_files[grep(sample[i], barcode_files)],
    features = feature_files[grep(sample[i], feature_files)]
  )
  
  mtx[[sample[i]]] <- new_element
}

# replace names
base::names(mtx) <- stringi::stri_replace_all_regex(
  base::names(mtx),
  pattern = prefix_donor$sample,
  replacement = prefix_donor$prefix,
  vectorize = FALSE)

# create seurat object
seurat_list <- BiocGenerics::Map(Seurat::CreateSeuratObject,
                                 counts = mtx,
                                 project = names(mtx))

# Rename barcodes using object names as prefix
for (i in names(seurat_list)) {
  seurat_list[[i]] <- Seurat::RenameCells(seurat_list[[i]],
                                          add.cell.id = i)
}

# merge seurat object
seurat_2 <- purrr::reduce(seurat_list, merge)

# remove -1 from cell names
seurat_3 <- Seurat::RenameCells(seurat_2, new.names = stringr::str_remove(Seurat::Cells(seurat_2), "-1"))

# remove uninteresting barcodes
seurat_3 <- seurat_3 %>%
  base::subset(cells = test$barcode)

# add meta data
seurat_3 <- seurat_3 %>%
  Seurat::AddMetaData(tibble::column_to_rownames(test, "barcode"))


# Calculate UCell scores --------------------------------------------------
# Subset object
Seurat::Idents(seurat_3) <- "Celltype"
s_beta <- seurat_3 %>% 
  base::subset(idents = "Beta")

gene_cluster_h <- purrr::map(gene_clusters, mouse2human_symbol)

# Convert mouse genes to human
## Add Ucell scores
s_beta <- UCell::AddModuleScore_UCell(obj = s_beta,
                                      features= gene_cluster_h,
                                      name = "_ucell",
                                      maxRank = 1000,
                                      BPPARAM = MulticoreParam(workers = parallel::detectCores() - 1))


# Plot ucell scores -------------------------------------------------------
arrangement <- c("Untreated", "Glucose", "Palmitate", "Glucose+Palmitate", "IL1b", "IL1b+IFNg", "IFNg", "IFNa")

# remove IFNg and Il1b as there is only data from one replicate and time point
data <- s_beta@meta.data %>% 
  dplyr::select(dplyr::ends_with("ucell"), Condition, Timepoint) %>% 
  tidyr::pivot_longer(c(-Condition, -Timepoint), names_to = "ucell", values_to = "score") %>% 
  dplyr::mutate(Condition = base::factor(Condition, levels = arrangement)) %>% 
  base::split(factor(.$ucell))

p4 <- purrr::map(data, function(data){
  
  # calculate median of untreated group
  
  min_value <- base::min(data$score)
  if (min_value < 0.38) {
    breakss = c(0.1, 0.35)
  } else{
    breakss = c(0.1, 0.45)
  }
  
  med <- data%>% 
    dplyr::filter(Timepoint == "24h") %>% 
    dplyr::filter(Condition == "Untreated") %>% 
    dplyr::pull() %>% 
    stats::median()
  
  p1 <- data %>% 
    dplyr::filter(Timepoint == "24h") %>% 
    ggplot2::ggplot(aes(y = score, x = Condition)) +
    ggplot2::geom_hline(yintercept = med, color = "red") +
    ggplot2::geom_boxplot(outlier.shape = NA, lwd = 0.4) +
    ggplot2::expand_limits(y = 0) +
    ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 7), position = "left") +
    ggbreak:: scale_y_break(breaks=breakss, scales = 4) +
    my_theme() +
    ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                   axis.title.y.right = element_blank(),
                   axis.text.y.right = element_blank(),
                   axis.ticks.y.right = element_blank(),
                   axis.line.y.right = element_blank())
  
  med <- data%>% 
    dplyr::filter(Timepoint == "72h") %>% 
    dplyr::filter(Condition == "Untreated") %>% 
    dplyr::pull() %>% 
    stats::median()
  
  p2 <- data %>% 
    dplyr::filter(Timepoint == "72h") %>% 
    ggplot2::ggplot(aes(y = score, x = Condition)) +
    ggplot2::geom_hline(yintercept = med, color = "red") +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::expand_limits(y = 0) +
    ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 7), position = "left") +
    ggbreak:: scale_y_break(breaks=breakss, scales = 4) +
    my_theme() +
    ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                   axis.title.y.right = element_blank(),
                   axis.text.y.right = element_blank(),
                   axis.ticks.y.right = element_blank(),
                   axis.line.y.right = element_blank())
  
  p3 <- p1/p2
  return(p3)
})

p5 <- wrap_plots(p4, nrow = 1)
ggplot2::ggsave(filename = here::here("data/public_data/glu_lip_cyt_human/GSE218316/processed/ucell_boxplot_geneclusters.pdf"),
                plot = p5, 
                width = 7,
                height = 5)

# Differential expressed genes --------------------------------------------
# add sample info
seurat_3@meta.data <- seurat_3@meta.data %>%
  as.data.frame%>%
  dplyr::mutate(sample = paste0(Condition, "_", Timepoint, "_", Donor),
                cond_time = paste0(Condition, "_", Timepoint),
                donor_time = paste0(Timepoint, "_", Donor))

Idents(seurat_3) <- "Celltype"

s_Beta <- subset(x = seurat_3, idents = "Beta")

counts <- edgeR::Seurat2PB(s_Beta, sample = "Donor", cluster = "cond_time")

counts_2 <- counts %>%
  as.data.frame() %>%
  dplyr::rename_with(~gsub("cluster", "", .x)) %>%
  tibble::column_to_rownames("gene")


meta_data <- counts[["samples"]] %>%
  dplyr::mutate(time = stringr::str_extract(cluster, "[:digit:]{2}[h]"),
                treatment = stringr::str_extract(cluster, "[^_]+"))

rownames(meta_data) <- gsub("cluster", "", rownames(meta_data))

# Define comparisons and initialize an empty list for results
# There is only 1 replicate in IFNg, so we had to exclude this one.
comparisons <- list(
  "Glucose" = "Glucose",
  "Palmitate" = "Palmitate",
  "Glucose.Palmitate" = "Glucose+Palmitate",
  "IL1b.IFNg" = "IL1b+IFNg",
  "IFNa" = "IFNa"
)

results_list <- list()

for (comp in names(comparisons)) {
  for (time_point in c("24h", "72h")) {
    # Subset the metadata and counts for the current comparison and time point
    meta_data_subset <- subset(meta_data, time == time_point & (treatment == "Untreated" | treatment == comparisons[[comp]]))
    counts_subset <- counts_2[, rownames(meta_data_subset)]
    
    # Set reference level for treatment
    meta_data_subset$treatment <- relevel(factor(meta_data_subset$treatment), ref = "Untreated")
    
    # Create DESeq2 dataset
    dds <- DESeqDataSetFromMatrix(countData = counts_subset, colData = meta_data_subset, design = ~ sample + treatment)
    
    # Run DESeq2 analysis
    dds <- DESeq(dds)
    
    # Extract results
    res <- results(dds, contrast = c("treatment", comparisons[[comp]], "Untreated"))
    
    # Store results in the list with appropriate name
    results_list[[paste0(comp, "_", time_point)]] <- res
  }
}

# Keep significant
sig_res <- results_list %>%
  purrr::modify_depth(1, ~ BiocGenerics::as.data.frame(.) %>%
                        dplyr::filter(padj <= 0.05 & log2FoldChange >= 0.5) %>%
                        tibble::rownames_to_column("gene"))

base::saveRDS(sig_res, here::here("data/public_data/glu_lip_cyt_human/GSE218316/processed/degs_sig.rds"))

# convert from human gene symbol to mouse gene symbol
# human to mouse id -------------------------------------------------------
mouse_genes <- BiocGenerics::Map(function(df){
  vec <- df$gene
  genes <- human2mouse(humanids = vec, keytype = "SYMBOL")
  output <- as.vector(na.omit(genes$Mouse_symbol))
  return(output)
},
df = sig_res)


# discard genes that are not expressed in seurat_beta
DefaultAssay(seurat_beta) <- "RNA"

genes_mouse_keep <- BiocGenerics::Map(function(vec) {
  output <- vec[!vec %in% setdiff(vec, rownames(seurat_beta))]
  return(output)
},
vec = mouse_genes)

# save seurat
saveRDS(seurat_3, here::here("data/public_data/glu_lip_cyt_human/GSE218316/processed/seurat_object.rds"))

# plot --------------------------------------------------------------------
arrangement <- c("IFNa_24h", "IL1b.IFNg_24h", "IFNa_72h", "IL1b.IFNg_72h")

openxlsx::write.xlsx(auc_gene_df, file = here::here("data/public_data/glu_lip_cyt_human/GSE218316/processed/auc.xlsx"))
openxlsx::write.xlsx(genes_mouse_keep, file = here::here("data/public_data/glu_lip_cyt_human/GSE218316/processed/genes.xlsx"))
