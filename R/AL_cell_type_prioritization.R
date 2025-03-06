# Description -------------------------------------------------------------
# Cell type prioritization using augur
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7610525/#SD13
# https://github.com/neurorestore/Augur
# Setup -------------------------------------------------------------------
base::source(here::here("R/set_up.R"))
set.seed(1000)
create_directories(here::here("data/cell_prioritization/"))

# Load --------------------------------------------------------------------
# seurat
seurat_5 <- base::readRDS(here::here("data/seurat_objects/seurat_5.rds"))

# # RNA ---------------------------------------------------------------------
### Expression data ----
exp <- Seurat::GetAssayData(seurat_5,
                            assay = "RNA")

### Metadata ----
meta <- seurat_5@meta.data %>%
  BiocGenerics::as.data.frame() %>%
  dplyr::select(condition, manual_anno)

# # all cells, set to 12 as we want to retain stellate cells
augur <- Augur::calculate_auc(input = exp,
                              meta = meta,
                              cell_type_col = "manual_anno",
                              label_col = "condition", 
                              subsample_size = 12,
                              min_cells = 12, 
                              n_subsamples = 150,
                              n_threads = 50)

# # Save --------------------------------------------------------------------
base::saveRDS(augur, here::here("data/cell_prioritization/augur_across_all_conditions.rds"))

# Do augur randomly 100 times ---------------------------------------------

# create metadataframe, where condition is randomly sampled each time
# set iteration
n_itr <- 100

# generation n_itr random seeds
random_seeds <- sample(1:1000000, n_itr, replace=TRUE)

meta_list <- list()
for (i in seq_len(n_itr)){
  
  df <- meta
  
  # set new seed each time you choose a sample
  
  set.seed(random_seeds[i])
  print(paste0("seed used ", random_seeds[i]))
  
  new_element <- df %>%
    dplyr::mutate(condition = sample(df$condition, size = length(df$condition)))
  
  meta_list[[paste0("iteration_", i)]] <- new_element
}

set.seed(1000)

augur_iteration <- BiocGenerics::Map(function(meta, iteration){
  
  print(paste0("Currently running AUGUR ", iteration))
  
  # run augur
  augur <- Augur::calculate_auc(input = exp,
                                meta = meta,
                                cell_type_col = "manual_anno",
                                label_col = "condition",
                                subsample_size = 8,
                                n_threads = 60)
  
  # extract results
  auc <- augur$AUC
  
  # get index with max auc
  max_index <- BiocGenerics::which.max(auc$auc)
  
  # extract row with max auc
  max_row <- auc %>%
    dplyr::mutate(iteration = iteration)
  
  print(paste0("Saving AUGUR iteration ", iteration))
  saveRDS(max_row, paste0(here::here("data/cell_prioritization/augur_iteration_"), iteration, ".rds"))
  
  
},
meta = meta_list,
iteration = names(meta_list))


# Calculate Emperical p-value ---------------------------------------------
# load real augur
augur <- readRDS(here::here("data/cell_prioritization/augur_across_all_conditions.rds"))

# Extract the auc for beta cells
observed_stat <- augur$AUC$auc[1]

# get path for iterations
augur_permutation <- list.files(pattern = "augur_iteration_iteration_", recursive = TRUE, full.names = TRUE)


permuted_stats <- BiocGenerics::Map(function(path){
  output <- readRDS(here::here(path))
  return(output)
},
augur_permutation) %>%
  dplyr::bind_rows()

cell_types <- as.list(unique(unfactor(permuted_stats$cell_type)))

emperical_p_value <- BiocGenerics::Map(function(cell_types){
  
  observed <- augur$AUC %>%
    dplyr::filter(cell_type == cell_types) %>%
    dplyr::select(auc) %>%
    purrr::as_vector() %>%
    unname()
  
  print(paste0("Observed AUC for ", cell_types, " is ", round(observed, 3)))
  
  permuted <- permuted_stats %>%
    dplyr::filter(cell_type == cell_types) %>%
    dplyr::select(auc) %>%
    purrr::as_vector() %>%
    unname()
  
  emperical_p <- sum(permuted >= observed) / n_itr
  
  print(paste0("Permutated AUCs for ", cell_types, " are: "))
  print(round(permuted, 3))
  
  print(paste0("Emperical P-value for ", cell_types, " is ", emperical_p))
  print(paste0("--------------------------------------------------------"))
},
cell_types = cell_types)


dplyr::filter(cell_type == "Beta") %>%
  dplyr::select(auc) %>%
  purrr::as_vector() %>%
  unname()

emperical_p <- sum(permuted_stats >= observed_stat) / n_itr
print(emperical_p)

