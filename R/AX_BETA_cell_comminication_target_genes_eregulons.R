# Description -------------------------------------------------------------
# "Is the occurrence of any eregulon target genes overrepresented in tnf and ifnb1 target genes (found by nichenet)
# compared to what would be expected by chance?"

# Setup -------------------------------------------------------------------
base::source(here::here("R/set_up.R"))
set.seed(1000)

# Load --------------------------------------------------------------------
# eregulon
eregulon_geneset <- read.csv(here::here("data/scenicplus/beta/results/scenicplus/export/eRegulon_metadata_filtered.csv"))

# high quality nichenet_genesets
high_q <- base::readRDS(here::here("data/scenicplus/beta/results/high_quality_eregulons.rds"))

# target genes
nichenet_target <- base::readRDS(here::here("data/nichenet/beta/files/hfd1_hfd3_tnf_ifnb1_target_genes.rds"))

# universe
universe <- readRDS(here::here("data/nichenet/beta/files/nichenet_universe.rds"))

# Prepare data ------------------------------------------------------------
eregu_hfd_keep <- c("Cebpg_+", "Stat1_+", 
                    "Hivep1_+", "Irf2_+", "Irf7_+",
                    "Irf9_+", "Nfkb1_+", "Pax6_+",
                    "Pura_+", "Relb_+",
                    "Stat2_+")

## hfd induced eregulon genes
eregulon_genes <- eregulon_geneset %>%
  dplyr::mutate(Consensus_name = stringr::str_replace(Consensus_name, "_[-+]$", "")) %>%
  dplyr::filter(Consensus_name %in% eregu_hfd_keep) %>%
  base::split(f = as.factor(.$Consensus_name)) %>%
  purrr::modify_depth(1, ~ dplyr::select(., Gene) %>%
                        purrr::as_vector() %>%
                        unname() %>%
                        unique())

## Nichenet target genes - upregulated genes
nichenet_genes <- nichenet_target %>%
  dplyr::filter(direction_regulation == "up") %>%
  dplyr::mutate(split = paste0(group, "_", id)) %>%
  base::split(f = as.factor(.$split)) %>%
  purrr::modify_depth(1, ~ dplyr::select(., target) %>%
                        purrr::as_vector() %>%
                        unname() %>%
                        unique())
# Create contingency table ------------------------------------------------
contingency_table <- BiocGenerics::Map(function(genes, name_genes, pathway, universe){
  
  print(paste0("creating contingency table for ", name_genes))
  
  output <- BiocGenerics::Map(function(pathway_1, universe_1){
    a <- length(intersect(genes, pathway_1)) # genes in both genes and pathway
    b <- length(setdiff(genes, pathway_1))  # genes in genes but not in pathway
    c <- length(setdiff(pathway_1, genes))  # genes in pathway but not in genes
    d <- length(setdiff(universe_1, union(genes, pathway_1))) # genes in neither
    
    output <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE,
                     dimnames = list(c("In genes", "Not in genes"),
                                     c("In pathway", "Not in pathway")))
    return(output)
  },
  pathway_1 = pathway,
  universe_1 = rep(list(universe), length(pathway)))
  
  return(output)
},
genes = nichenet_genes,
name_genes = names(nichenet_genes),
pathway = rep(list(eregulon_genes), length(nichenet_genes)),
universe = rep(list(universe), length(nichenet_genes)))

# save contingency table
base::saveRDS(contingency_table, here::here("data/nichenet/beta/files/contingency_table.rds"))

# Significance test -------------------------------------------------------
fisher_test <- BiocGenerics::Map(function(tbl){
  output <- BiocGenerics::Map(function(tbl_con){
    output <-
      stats::fisher.test(tbl_con)
    return(output)
  },
  tbl_con = tbl)
  
  return(output)
},
tbl = contingency_table)


# Extract results and multipletesting -------------------------------------
fish_test_adj <- BiocGenerics::Map(function(tbl, nichenet_geneset){
  
  nichenet_geneset <- nichenet_geneset
  pval <- BiocGenerics::Map(function(tbl_con, geneset){
    output <- tbl_con %>%
      broom::tidy() %>%
      dplyr::mutate(geneset = geneset,
                    nichenet_geneset = stringr::str_replace(nichenet_geneset, "_[-+]$", ""))
    return(output)
  },
  tbl_con = tbl,
  geneset = names(tbl))
  
  # bind rows for all geneset results
  # from apo_a1.
  output <- bind_rows(pval)
  return(output)
},
tbl = fisher_test,
nichenet_geneset = names(fisher_test)) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(fdr =  stats::p.adjust(p.value, method = "BH"))

# Plot eregulons  ---------------------------------------------------------
nichenet_df <- fish_test_adj %>%
  dplyr::mutate(logfdr = -log10(fdr),
                condition = stringr::str_extract(nichenet_geneset, "^[^_]*")) %>%
  tidyr::pivot_wider(id_cols = c(nichenet_geneset, condition),
                     names_from = geneset,
                     values_from = logfdr) %>%
  tibble::column_to_rownames("nichenet_geneset")

nichenet_mat <- nichenet_df %>%
  dplyr::select(-condition) %>%
  as.matrix()

pdf(here::here("data/scenicplus/beta/results/hfd_hfd3_tnf_ifbn1_target_genes_ORA_eregulons_heatmap.pdf"),
    height = 3,
    width = 5.5)
# clustering using the UPGMC method
ComplexHeatmap::pheatmap(
  mat = nichenet_mat,
  color = c("white", "#fa8231"),
  breaks = c(0, -log10(0.05), max(nichenet_mat)),
  main = "Nichenet Ligand target genes",
  show_rownames = TRUE,
  cluster_cols = TRUE,
  cluster_rows = TRUE,
  heatmap_legend_param = list(title = "-log10(FDR)"),
  border_color = "#808285")

dev.off()


# Save results ------------------------------------------------------------
contingency_table_df <- BiocGenerics::Map(function(tbl, eregulon){
  
  eregulon <- eregulon
  pval <- BiocGenerics::Map(function(tbl_con, geneset){
    output <- tbl_con %>%
      BiocGenerics::as.data.frame() %>%
      tibble::rownames_to_column("category") %>%
      tidyr::pivot_wider(names_from = category, values_from = c(`In pathway`, `Not in pathway`)) %>%
      dplyr::mutate(geneset = geneset,
                    interaction = eregulon)
    return(output)
  },
  tbl_con = tbl,
  geneset = names(tbl))
  
  # bind rows for all geneset results
  # from apo_a1.
  output <- bind_rows(pval)
  return(output)
},
tbl = contingency_table,
eregulon = names(contingency_table)) %>%
  dplyr::bind_rows()


nichenet_genes_df <- nichenet_genes %>% purrr::map_df(base::data.frame, .id = "geneset") %>%
  `colnames<-`(c("interaction", "genes")) %>%
  dplyr::group_by(interaction) %>%
  summarise(genes = paste(genes, collapse=", "))

contingency_table_df <- contingency_table_df %>%
  dplyr::left_join(y = nichenet_genes_df, by = c("interaction")) %>%
  dplyr::left_join(y = fish_test_adj %>% dplyr::rename(interaction = nichenet_geneset), by = c("geneset", "interaction")) %>%
  dplyr::rename(eregulon = geneset)

saveRDS(contingency_table_df, here::here("data/scenicplus/beta/results/contingency_table_df.rds"))
