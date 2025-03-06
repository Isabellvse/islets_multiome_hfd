# Description -------------------------------------------------------------
# See if genes are positively or negatively correlated with a clinical trait is enriched in eRegulon target genes, using a GSEA analysis.
# Setup -------------------------------------------------------------------
base::source(here::here("R/set_up.R"))
set.seed(1000)

create_directories(c(here::here("data/gsea")))


# Load --------------------------------------------------------------------
# eregulon
eregulon_geneset <- read.csv(here::here("data/scenicplus/beta/results/scenicplus/export/eRegulon_metadata_filtered.csv"))

# high quality nichenet_genesets
high_q <- readRDS(here::here("data/scenicplus/beta/results/high_quality_eregulons.rds"))

# Attie correlated genes
attie_cor <- readRDS(here::here("data/attie_diabetes_database/genomic_study_f2_cohort/correlation/Cor_list_extended.rds"))

# attie annotation
attie_anno <- readxl::read_excel(here::here("data/attie_diabetes_database/genomic_study_f2_cohort/correlation/categories_2.xlsx"))

# prepare data ------------------------------------------------------------
## eregulon genes
eregulon_genes <- eregulon_geneset %>%
    dplyr::mutate(Consensus_name = stringr::str_replace(Consensus_name, "_[-+]$", "")) %>%
    dplyr::filter(Consensus_name %in% high_q) %>%
    base::split(f = as.factor(.$Consensus_name)) %>%
    purrr::modify_depth(1, ~ dplyr::select(., Gene) %>%
                            purrr::as_vector() %>%
                            unname() %>%
                            unique())

# GSEA analysis -----------------------------------------------------------
## Create a named vector of genes for each each cell, where they are sorted by correlation
attie_genes_cor_2 <- purrr::map(attie_cor, function(df){
  vec <- setNames(df$correl_coeff, df$transcript)
  vec_sort <- BiocGenerics::sort(vec, decreasing = TRUE)
  
  return(vec_sort)
})
  
## Create term2gene df with eregulons
term2gene <- eregulon_genes %>%
    purrr::imap_dfr(~ tibble::tibble(term = .y, gene = .x))

## Extract names of genesets, which are in the categories we want to test
category <- attie_anno %>%
    dplyr::filter(Category %in% c("Metabolism and Energy Regulation",
                                  "Inflammation and Immune Response"))
# subset attie_gene_cor
attie_sub <- attie_genes_cor_2[category$id]

## Run GSEA ----
# Loop to run correlation test for all the clinical traits and all the eregulon target genes
gse <- purrr::map(attie_sub, ~clusterProfiler::GSEA(geneList = .x,
                                                          TERM2GENE = term2gene,
                                                          pvalueCutoff = 1))

# put everything into a dataframe and do p-value adjustment
gse_df <- purrr::map2_df(gse, names(gse),
                          ~ BiocGenerics::as.data.frame(.x@result) %>%
                              dplyr::mutate(GeneSet = .y)) %>%
    dplyr::mutate(fdr =  stats::p.adjust(pvalue, method = "BH")) %>%
    dplyr::relocate(GeneSet, ID, fdr)


## only significant ----
gse_df_sig <- gse_df %>%
    dplyr::filter(fdr <= 0.05)

## which is the most significant
most_sig <- gse_df_sig %>%
    dplyr::arrange(fdr, desc(abs(NES))) %>%
    head(n = 10)


# nes heatmap -------------------------------------------------------------
## Prepare data ----
# Clinical traits that are significantly enriched for eRegulon target genes
sig_genesets <- gse_df_sig %>%
    dplyr::select(GeneSet) %>%
    purrr::as_vector() %>%
    unname() %>%
    unique()

# Dataframe with nes scores, NA values will be replaced with 0
nes_df <- gse_df %>%
    dplyr::select(GeneSet, ID, NES) %>%
    dplyr::filter(GeneSet %in% sig_genesets) %>%
    tidyr::pivot_wider(id_cols = GeneSet, names_from = ID, values_from = NES) %>%
    tibble::column_to_rownames("GeneSet")

# dataframe with FDR scores
fdr_df <-  gse_df %>%
    dplyr::select(GeneSet, ID, fdr) %>%
    dplyr::filter(GeneSet %in% sig_genesets) %>%
    tidyr::pivot_wider(id_cols = GeneSet, names_from = ID, values_from = fdr) %>%
    tibble::column_to_rownames("GeneSet")

# Make sure everything is equal
all(rownames(fdr_df) == rownames(nes_df))
all(colnames(fdr_df) == colnames(nes_df))
dim(nes_df)
dim(fdr_df)

# Determine significance based on FDR
significant <- fdr_df <= 0.05

# Map NES values to color based on significance and value
color_matrix <- matrix("white", nrow = nrow(nes_df), ncol = ncol(nes_df))

# Set colors based on NES values
color_matrix[nes_df > 0 & significant] <- "#A83708"  # Positive and significant
color_matrix[nes_df < 0 & significant] <- "#004B7A"   # Negative and significant
color_matrix[is.na(nes_df) == TRUE | !significant] <- "white" # NA or non-significant

# Create a discrete matrix for color mapping
discrete_mat <- matrix(0, nrow = nrow(nes_df), ncol = ncol(nes_df))
rownames(discrete_mat) <- rownames(nes_df)
colnames(discrete_mat) <- colnames(nes_df)

# Assign discrete values based on NES and significance
discrete_mat[nes_df > 0 & significant] <- 2  # Positive and significant
discrete_mat[nes_df < 0 & significant] <- 1  # Negative and significant
discrete_mat[is.na(nes_df) == TRUE | !significant] <- 3  # Zero or non-significant

nes_annotation <- discrete_mat %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("id") %>% 
  dplyr::left_join(y = attie_anno, by = "id") %>% 
  dplyr::select(id, full_name, Category) %>% 
  tibble::column_to_rownames("id") %>% 
  dplyr::mutate(Category = snakecase::to_sentence_case(Category))

# Define a custom color mapping
colors <- structure(c("#004B7A", "#A83708", "white"), names = c("1", "2", "3"))

pdf(here::here("data/gsea/GSEA_glucose_and_immune_sets.pdf"),
    height = 4,
    width = 5)
# Create a Heatmap object with the custom color mapping
ComplexHeatmap::Heatmap(
    discrete_mat,
    name = "NES",
    col = colors,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    row_title = "Genes",
    column_title = "Conditions",
    row_split = nes_annotation$Category,
    row_labels = nes_annotation$full_name, border = "black",
    row_names_gp = gpar(fontsize = 3),
    rect_gp = gpar(col = "grey", lwd = 1),
    heatmap_legend_param = list(title = "NES", at = 1:3, labels = c("Negative", "Positive", "Non-Significant/Zero"))
)
dev.off()

## HFD induced, and only those that have a significant hit
eregu_hfd_keep <- colnames(discrete_mat)[which(colnames(discrete_mat) %in% c("Cebpg_+", "Stat1_+", "Hivep1_+", "Irf2_+", "Irf7_+",
                    "Irf9_+", "Nfkb1_+", "Pax6_+",
                    "Pura_+", "Relb_+",
                    "Stat2_+"))]

geneset_hfd_sig <- gse_df %>%
    dplyr::filter(ID %in% eregu_hfd_keep & fdr <= 0.05) %>%
    dplyr::select(GeneSet) %>%
    purrr::as_vector() %>%
    unique()

nes_annotation_hfd <- nes_annotation %>%
    dplyr::filter(rownames(.) %in% geneset_hfd_sig)

# check that geneset_hfd_sig and nes_annotation_hfd are the same order
all.equal(rownames(nes_annotation_hfd), geneset_hfd_sig)

pdf(here::here("data/gsea/GSEA_glucose_and_immune_sets_hfd_sig.pdf"),
    height = 6,
    width = 5)

ComplexHeatmap::Heatmap(
    discrete_mat[geneset_hfd_sig, eregu_hfd_keep],
    name = "NES",
    col = colors,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    row_title = "Genes",
    column_title = "Conditions",
    row_split = nes_annotation_hfd $Category,
    row_labels = nes_annotation_hfd$full_name, border = "black",
    row_names_gp = gpar(fontsize = 6),
    rect_gp = gpar(col = "grey", lwd = 1),
    heatmap_legend_param = list(title = "NES", at = 1:3, labels = c("Negative", "Positive", "Non-Significant/Zero"))
)
dev.off()


# plot a single pathway ---------------------------------------------------

pdf(here::here("data/gsea/GSEA_examples_first_hit.pdf"),
    height = 6,
    width = 7)
for (i in 1:length(gse)){
    j = names(gse)[i]
    print(gse_plot(x = gse[[j]],
                   correlation_vector = attie_genes_cor_2[[j]],
                   geneSetID = 1,
                   title = paste(gse[[j]]$Description[1], "target genes enriched in", j)))
}
dev.off()


# plot for a single pathway and gene set
gse_plot(x = gse[["LDL"]], correlation_vector = attie_genes_cor[["LDL"]], geneSetID = 1, title = paste(gse[["LDL"]]$Description[1], "target genes overrepresented in", names(gse)[1]))


# save --------------------------------------------------------------------
saveRDS(attie_genes_cor_2, here::here("data/gsea/attie_gene_list_sorted.rds"))
saveRDS(term2gene, here::here("data/gsea/term2gene.rds"))
saveRDS(gse, here::here("data/gsea/gse_results_list.rds"))
saveRDS(gse_df, here::here("data/gsea/gse_results_dataframe.rds"))
saveRDS(gse_df_sig, here::here("data/gsea/gse_results_dataframe_sig.rds"))
saveRDS(most_sig, here::here("data/gsea/top_10_most_sig_NES.rds"))


