# Description -------------------------------------------------------------
# Exporting excel results from various analysis

# Setup -------------------------------------------------------------------
base::source(here::here("R/set_up.R"))
set.seed(1000)

create_directories(c(here::here("data/export/")))
seurat_5 <- readRDS(here::here("data/seurat_objects/seurat_5.rds"))
seurat_beta <-  readRDS(here::here("data/inflammatory_cluster/files/seurat_beta.rds"))

# Load and process --------------------------------------------------------
table_8 <- readRDS(here::here("data/deseq2/celltype/files/marker_genes_results.rds")) %>% 
  purrr::modify_depth(1, ~ BiocGenerics::as.data.frame(.) %>%
                        dplyr::filter(padj < 0.05 & log2FoldChange > 0) %>%
                        dplyr::select(-baseMean) %>%
                        tibble::rownames_to_column("Gene") %>%
                        dplyr::rename("Log2(FoldChange)" = log2FoldChange,
                                      "Standard Error" = lfcSE,
                                      "Wald statistic" = stat,
                                      "Wald test p-value" = pvalue,
                                      "FDR" = padj)) %>% 
  dplyr::bind_rows(.id = "Cell Type") %>% 
  dplyr::relocate("Cell Type")


## fuzzy cluster, clusters
cl <- readRDS(here::here("data/deseq2/rna_clusterbycond/files/clustering_new.rds"))
lrt_cluster <- BiocGenerics::Map(heatmap_cluster,
                                 df = cl) %>%
  purrr::modify_depth(1, ~ dplyr::rename(., "Gene" = gene,
                                         "Fuzzy Cluster ID" = membership)) %>% 
  dplyr::bind_rows(.id = "Cell Type")

## Pseudobulk Normalized counts

## diff genes across conditions
table_10 <- readRDS(here::here("data/deseq2/rna_clusterbycond/files/res_lrt_significant.rds")) %>% 
  purrr::modify_depth(1, ~ BiocGenerics::as.data.frame(.) %>%
                      dplyr::select(-baseMean, -log2FoldChange, -lfcSE, - stat) %>%
                      tibble::rownames_to_column("Gene") %>%
                      dplyr::rename("LRT" = pvalue,
                                    "FDR" = padj)) %>%
  dplyr::bind_rows(.id = "Cell Type") %>% 
  dplyr::left_join(y = lrt_cluster, by = c("Cell Type", "Gene")) %>% 
  dplyr::relocate(c("Cell Type", "Fuzzy Cluster ID"))

## pathway analysis across conditions
table_11 <- readRDS(here::here("data/deseq2/rna_clusterbycond/files/clusterprofiler_go_biolgoical_process.rds")) %>% 
  purrr::modify_depth(2, ~ BiocGenerics::as.data.frame(.)) %>% 
  purrr::modify_depth(1, ~ dplyr::bind_rows(., .id = "Fuzzy Cluster ID")) %>% 
  purrr::pluck("Beta") %>% 
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::select(-k, -n, -M, -N, -is_sig, -j_bg) %>%
  dplyr::rename("GO-term ID" = ID,
                "Gene Ratio" = GeneRatio,
                "Background Ratio" = BgRatio,
                "p-value" = pvalue,
                "FDR" = p.adjust,
                "q-value" = qvalue,
                "Genes" = geneID,
                "Jaccard Similary Index" = "j_path")

## diff genes immune reponse subcluster

table_19<- readRDS(here::here("data/inflammatory_cluster/files/deseq2.rds")) %>% 
  dplyr::bind_rows(.id = "Comparisons") %>% 
  dplyr::relocate("Comparisons") %>% 
  dplyr::rename(Gene = gene,
                "Log2(FoldChange)" = log2FoldChange,
                "Standard Error" = lfcSE,
                "Wald statistic" = stat,
                "Wald test p-value" = pvalue,
                "FDR" = padj)

immune_go_term_up <- readRDS(here::here("data/inflammatory_cluster/files/go_term_analysis_upregulated.rds")) %>% 
  purrr::modify_depth(1, ~ BiocGenerics::as.data.frame(.) %>% 
                        dplyr::mutate("Direction" = "Up")) %>%
  dplyr::bind_rows(.id = "Comparisons")

immune_go_term_down <- readRDS(here::here("data/inflammatory_cluster/files/go_term_analysis_downregulated.rds")) %>% 
  purrr::modify_depth(1, ~ BiocGenerics::as.data.frame(.) %>% 
                        dplyr::mutate("Direction" = "Down")) %>%
  dplyr::bind_rows(.id = "Comparisons")

table_20 <- dplyr::full_join(x = immune_go_term_up, y = immune_go_term_down) %>% 
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::select(-k, -n, -M, -N, -is_sig, -j_bg) %>%
  dplyr::rename("GO-term ID" = ID,
                "Gene Ratio" = GeneRatio,
                "Background Ratio" = BgRatio,
                "p-value" = pvalue,
                "FDR" = p.adjust,
                "q-value" = qvalue,
                "Genes" = geneID,
                "Jaccard Similary Index" = "j_path")

## nichenet interactions
multinichenet_output <- readRDS(here::here("data/nichenet/beta/files/beta_vs_all_multinichenet.rds"))

lfd_30 <-  multinichenet_output$prioritization_tables$group_prioritization_tbl %>% 
  dplyr::filter(group == "LFD" & fraction_expressing_ligand_receptor > 0) %>% 
  dplyr::distinct(group, sender, receiver, ligand, receptor, receiver, id, prioritization_score) %>% 
  dplyr::mutate(prioritization_rank = rank(-prioritization_score)) %>% 
  dplyr::filter(prioritization_rank <= 30)

hfd1_30 = multinichenet_output$prioritization_tables$group_prioritization_tbl %>% 
  dplyr::filter(group == "HFD1" & fraction_expressing_ligand_receptor > 0) %>% 
  dplyr::distinct(group, sender, receiver, ligand, receptor, receiver, id, prioritization_score) %>% 
  dplyr::mutate(prioritization_rank = rank(-prioritization_score)) %>% 
  dplyr::filter(prioritization_rank <= 30)

hfd3_30 = multinichenet_output$prioritization_tables$group_prioritization_tbl %>% 
  dplyr::filter(group == "HFD3" & fraction_expressing_ligand_receptor > 0) %>% 
  dplyr::distinct(group, sender, receiver, ligand, receptor, receiver, id, prioritization_score) %>% 
  dplyr::mutate(prioritization_rank = rank(-prioritization_score)) %>% 
  dplyr::filter(prioritization_rank <= 30)

table_16 <- purrr::reduce(list(lfd_30, hfd1_30, hfd3_30), dplyr::full_join) %>% 
  dplyr::rename_with(.fn = ~ snakecase::to_title_case(.)) %>% 
  dplyr::rename(ID = Id,
                Condition = Group)

## nichenet target genes
table_17 <- readRDS(here::here("data/nichenet/beta/files/hfd1_hfd3_tnf_ifnb1_target_genes.rds")) %>% 
  dplyr::rename(ID = id,
                "Pearson correlation coeffcient" = pearson,
                "Pearsons correlation P value" = pearson_pval,
                "Spearman correlation coefficient" = spearman,
                "Spearman correlation P value" = spearman_pval,
                "Prior knowledge score" = prior_score,
                "Direction of regulation" = direction_regulation,
                Condition = group
                ) %>% 
  dplyr::rename_with(.fn = ~ snakecase::to_title_case(.)) 

## meta data eregulons
table_13 <- read.csv(here::here("data/scenicplus/beta/results/scenicplus/export/eRegulon_metadata_filtered.csv")) %>%
  dplyr::select(-X)

high_q <- readRDS(here::here("data/scenicplus/beta/results/high_quality_eregulons.rds")) %>% 
  paste(collapse="|")

table_12 <- table_13 %>% dplyr::filter(base::grepl(high_q, Consensus_name)) %>% 
  dplyr::select(TF, Region_signature_name, Gene_signature_name, is_extended) %>% 
  dplyr::distinct() %>% 
  dplyr::mutate(Region_signature_name = stringr::str_extract(Region_signature_name, "(?<=\\()\\d+(?=r\\))"),
                Gene_signature_name = stringr::str_extract(Gene_signature_name, "(?<=\\()\\d+(?=g\\))")) %>%
  dplyr::rename("Transcription Factor" = TF,
                "Numer of Target Enhancers" = Region_signature_name,
                "Numer of Target Genes" = Gene_signature_name)

## eregulon correlation
pseudo_eregulon <- readRDS(here::here("data/scenicplus/beta/results/pseudo_counts_correlation.rds")) %>% 
  dplyr::select(eregulon, tf, product, condition, auc, value)

table_14 <- readRDS(here::here("data/scenicplus/beta/results/pseudo_counts_correlation_test.rds")) %>% 
  dplyr::left_join(y = pseudo_eregulon) %>% 
  dplyr::rename("Gene or Region AUC" = auc,
                "AUC score" = value,
                "TF * Activity product" = product,
                eRegulon = eregulon,
                "Upper confidence interval" = Confidence_interval_upper,
                "Lower condifence interval" = Confidence_interval_lower,
                "P value" = P_value,
                TF = tf,
                "T statistic" = "T_statistic",
                "Pearson correlation coefficient" = Correlation_coefficient)

## Attie diabetes database top 5 % of positively and negatively correlated genes
attie_genes <- readRDS(here::here("data/attie_diabetes_database/genomic_study_f2_cohort/correlation/attie_genesets_list.rds"))

## contingency table attie
table_18 <- readRDS(here::here("data/nichenet/beta/files/contingency_table_df.rds")) %>% 
  dplyr::rename_with(.fn = ~ snakecase::to_sentence_case(.)) %>% 
  dplyr::rename("Upper confidence interval" = "Conf high",
                "Lower confidence interval" = "Conf low",
                FDR = Fdr,
                eRegulon = Eregulon)

## attie anno
attie_anno <- readxl::read_excel(here::here("data/attie_diabetes_database/genomic_study_f2_cohort/correlation/categories_2.xlsx"))

## attie database GSEA
table_15<- readRDS(here::here("data/gsea/gse_results_dataframe.rds")) %>% 
  dplyr::select(-p.adjust) %>% 
  dplyr::rename("Gene set" = GeneSet,
                "eRegulon" = ID,
                "Gene set size" = setSize,
                FDR = fdr,
                "Enrichment score" = enrichmentScore,
                "Normalized enrichment score" = NES,
                "P value" = pvalue,
                "Core enrichment" = core_enrichment,
                Rank = rank,
                "Leading edge" = leading_edge)


meta_data <- seurat_5@meta.data %>%
  BiocGenerics::as.data.frame()

table_9 <- meta_data %>%
  dplyr::group_by(manual_anno) %>%
  dplyr::summarise(Number_of_nuclei_per_cluster = n()) %>%
  dplyr::left_join(meta_data %>%
                     dplyr::group_by(manual_anno, condition) %>%
                     dplyr::summarise(Number_of_nuclei_per_condition = n()), by = "manual_anno") %>%
  dplyr::left_join(meta_data %>%
                     dplyr::group_by(manual_anno, condition, orig.ident) %>%
                     dplyr::summarise(Number_of_nuclei_per_replicate = n()), by = c("manual_anno", "condition")) %>%
  dplyr::select(manual_anno, Number_of_nuclei_per_cluster, condition, Number_of_nuclei_per_condition, orig.ident, Number_of_nuclei_per_replicate)

table_24 <- readRDS(here::here("data/inflammatory_cluster/files/genes_by_pathway_perk_ire1_atf6.rds")) %>% 
  purrr::modify_depth(1, ~ as.data.frame(.) %>% 
                        dplyr::rename("Gene" = ".")) %>%  
  dplyr::bind_rows(.id = "Geneset") %>% 
  dplyr::mutate(Geneset = dplyr::case_when(Geneset == "upr" ~ "UPR",
                                            Geneset == "ins" ~ "Beta-cell identity"))

table_26 <- seurat_5@meta.data %>% BiocGenerics::as.data.frame()
table_27 <- seurat_beta@meta.data %>% BiocGenerics::as.data.frame()
# combine the tables ------------------------------------------------------
table_list <- list("Table 8" = table_8,
                   "Table 9" = table_9,
                   "Table 10" = table_10,
                   "Table 11" = table_11,
                   "Table 12" = table_12,
                   "Table 13" = table_13,
                   "Table 14" = table_14,
                   "Table 15" = table_15,
                   "Table 16" = table_16,
                   "Table 17" = table_17,
                   "Table 18" = table_18,
                   "Table 19" = table_19,
                   "Table 20" = table_20,
                   "Table 24" = table_20,
                   "Table 26" = table_20,
                   "Table 27" = table_20)

# save --------------------------------------------------------------------
openxlsx::write.xlsx(table_list, here::here("data/export/tables.xlsx"))                   

