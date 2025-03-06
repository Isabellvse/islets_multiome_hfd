# Description -------------------------------------------------------------
# Create geneset module of genes invovled in immune system, and identify beta-cell subcluster
# Setup -------------------------------------------------------------------
base::source(here::here("R/set_up.R"))
set.seed(1000)

create_directories(c(here::here("data/inflammatory_cluster/"),
                     here::here("data/inflammatory_cluster/files/"),
                     here::here("data/inflammatory_cluster/figures/")))

# Load --------------------------------------------------------------------
# seurat object
seurat_5 <- base::readRDS(here::here("data/seurat_objects/seurat_5.rds"))

# nichenet genes
nichenet_target <- base::readRDS(here::here("data/nichenet/beta/files/hfd1_hfd3_tnf_ifnb1_target_genes.rds"))

# preprocess --------------------------------------------------------------
## Nichenet target genes - upregulated genes
nichenet_genes <- nichenet_target %>%
  dplyr::filter(direction_regulation == "up") %>%
  dplyr::pull(target) %>% 
  base::unique()


nichenet_genes  <- nichenet_genes %>% stringr::str_replace(pattern = "[.]", replacement = "-") # replace "." with "-" in gene names

saveRDS(nichenet_genes, here::here("data/inflammatory_cluster/files/all_nichenet_gene.rds"))

## subset seurat ----
Idents(seurat_5) <- "manual_anno"
seurat_beta <- base::subset(seurat_5, idents = "Beta")

# Module scores -----------------------------------------------------------
seurat_beta <- UCell::AddModuleScore_UCell(obj = seurat_beta,
                                      features= list("inflammatory" = nichenet_genes),
                                      name = "_ucell",
                                      maxRank = 1000,
                                      BPPARAM = MulticoreParam(workers = parallel::detectCores() - 1))
# Find subclusters --------------------------------------------------------

x = sort(seurat_beta@meta.data[seurat_beta@meta.data$condition != "LFD","inflammatory_ucell"], decreasing = TRUE)
threshold = x[uik(y = x, x= 1:length(x))]

# ECDF Plot
ecdf_plot <- ggplot2::ggplot(seurat_beta@meta.data %>% dplyr::filter(condition %in% c("LFD", "HFD_1", "HFD_3")), 
                             aes(x = inflammatory_ucell, color = condition)) +
  ggplot2::stat_ecdf(size = 1) +
  ggplot2::geom_vline(xintercept = threshold, linetype = "dashed") +
  ggplot2::scale_y_reverse() +
  ggplot2::labs(x = "Score", y = "% of cells") +
  scale_color_manual(values = c("#004B7A", "#c71500", "#fcb21c")) +
  my_theme()

# Barplot
bar_data <- seurat_beta@meta.data %>%
  mutate(condition = factor(condition, levels = c("LFD", "HFD_1", "HFD_3"))) %>%
  group_by(condition) %>%
  summarise(proportion = sum(inflammatory_ucell >= threshold) / n())

bar_plot <- ggplot2::ggplot(bar_data, aes(x = condition, y = proportion, fill = condition)) +
  ggplot2::geom_bar(stat = "identity") +
  ggplot2::scale_fill_manual(values = c("#004B7A", "#c71500", "#fcb21c")) +
  ggplot2::ylim(0, 0.5) +
  ggplot2::labs(x = "", y = "Proportion") +
  my_theme()

# Save to PDF
ggplot2::ggsave(here::here("data/inflammatory_cluster/figures/inflammatory_beta_cells_definition.pdf"),
       plot = gridExtra::grid.arrange(ecdf_plot, bar_plot, ncol = 2),
       width = 4, height = 1)
table(seurat_beta@meta.data[seurat_beta@meta.data$inflammatory_ucell >= threshold, "condition"]) / table(seurat_beta@meta.data[, "condition"])

saveRDS(threshold, here::here("data/inflammatory_cluster/files/threshold.rds"))

# add inflammatory to seurat object ---------------------------------------
seurat_beta@meta.data <- seurat_beta@meta.data %>% 
  BiocGenerics::as.data.frame() %>% 
  dplyr::mutate(is_inflammatory_ucell = dplyr::case_when(inflammatory_ucell >= threshold ~ "high",
                                                      inflammatory_ucell < threshold ~ "low"))
# Differential expressed genes --------------------------------------------
counts <- seurat_beta %>%
  edgeR::Seurat2PB(sample = "orig.ident",
                   cluster = "is_inflammatory_ucell")
counts_2 <- counts %>%
  as.data.frame() %>%
  dplyr::rename_with(~gsub("cluster", "", .x)) %>%
  tibble::column_to_rownames("gene")

meta_data <- counts[["samples"]] %>%
  magrittr::set_rownames(base::gsub("cluster", "", BiocGenerics::rownames(.))) %>%
  dplyr::mutate(condition = stringr::str_extract(BiocGenerics::rownames(.), "LFD|HFD_1|HFD_3"),
                condition_inflammatory = base::factor(base::paste0(condition, "_", cluster)))

base::all.equal(BiocGenerics::colnames(counts_2), BiocGenerics::rownames(meta_data))

# Define the comparisons to be performed
comparisons <- list(
  "HFD_3_high_vs_HFD_3_low" = list("subset" = "HFD_3", "design" = "~ sample + cluster", "contrast" = c("cluster", "high", "low")),
  "HFD_1_high_vs_HFD_1_low" = list("subset" = "HFD_1", "design" = "~ sample + cluster", "contrast" = c("cluster", "high", "low")),
  
  "HFD_3_high_vs_HFD_1_high" = list("subset" = "high", "design" = "~ condition", "contrast" = c("condition", "HFD_3", "HFD_1")),
  
  "HFD_1_low_vs_LFD_low"    = list("subset" = "low", "design" = "~ condition", "contrast" = c("condition", "HFD_1", "LFD")),
  "HFD_3_low_vs_LFD_low"    = list("subset" = "low", "design" = "~ condition", "contrast" = c("condition", "HFD_3", "LFD")),
  "HFD_3_low_vs_HFD_1_low"  = list("subset" = "low", "design" = "~ condition", "contrast" = c("condition", "HFD_3", "HFD_1"))
)

# Loop through the comparisons and subset accordingly
results <- list()
wald_results <- list()
dds_list <- list()
dds_list_sub <- list()

for (comp_name in base::names(comparisons)) {
  comparison <- comparisons[[comp_name]]
  
  # Subset the data using dplyr
  subset_meta_data <- meta_data %>%
    dplyr::filter(
      (comparison$subset == "HFD_1" & condition == "HFD_1") |
        (comparison$subset == "HFD_3" & condition == "HFD_3") |
        (comparison$subset == "low" & cluster == "low") |
        (comparison$subset == "high" & cluster == "high")
    )
  
  subset_counts <- counts_2[, BiocGenerics::rownames(subset_meta_data)]
  
  # Create DESeq2 dataset with the specified design
  dds_subset <- DESeq2::DESeqDataSetFromMatrix(
    countData = subset_counts,
    colData = subset_meta_data,
    design = stats::as.formula(comparison$design)
  )
  
  # Normalize and run DESeq
  dds_subset <- DESeq2::DESeq(dds_subset)
  
  # Run Wald test for the specific contrast
  res <- DESeq2::results(dds_subset, contrast = comparison$contrast)
  
  # Store the result in the list
  wald_results[[comp_name]] <- res
  dds_list_sub[[comp_name]] <- dds_subset
}

results <- wald_results
dds_list <- dds_list_sub


# Pathway analysis --------------------------------------------------------
go_res_up <- list()
go_res_down <- list()

for (j in names(results)) {
  universe <- results[[j]]
  
  # Upregulated and downregulated genes
  up <- results[[j]] %>%
    as.data.frame() %>%
    filter(padj <= 0.05 & log2FoldChange > 0)
  
  down <- results[[j]] %>%
    as.data.frame() %>%
    filter(padj <= 0.05 & log2FoldChange < 0)
  
  # Skip if there are no upregulated or downregulated genes
  if (nrow(up) > 0) {
    go_up <- go_term(test_genes = up, bg_genes = universe)
    go_res_up[[j]] <- go_up
  }
  
  if (nrow(down) > 0) {
    go_down <- go_term(test_genes = down, bg_genes = universe)
    go_res_down[[j]] <- go_down
  }
}

base::saveRDS(go_res_up, file = here::here("data/inflammatory_cluster/files/go_term_analysis_upregulated.rds"))
base::saveRDS(go_res_down, file = here::here("data/inflammatory_cluster/files/go_term_analysis_downregulated.rds"))

## get top 10 pathways
## GO TERM ----
go_term_top10 <- purrr::map(list("up" = go_res_up, "down" = go_res_down), function(res) {
  
  purrr::map(res, function(res_2) {
    # Top 5 most significant pathways
    top_10_pathways <- res_2@result %>%
      tibble::rownames_to_column("go_id") %>% 
      dplyr::arrange(p.adjust) %>%
      utils::head(n = 10)  # Select top 5 most significant pathways
    
    return(top_10_pathways)
  })
}) %>% 
  purrr::modify_depth(1, ~ dplyr::bind_rows(., .id = "comparison")) %>% 
  dplyr::bind_rows(.id = "direction")

openxlsx::write.xlsx(go_term_top10, here::here("data/inflammatory_cluster/files/top10go.rds"))

# significant results -----------------------------------------------------
results_sig <- list()
res_list <- list()
for (j in names(results)) {
  res <- results[[j]] %>%
    as.data.frame() %>%
    rownames_to_column("gene")
  
  res_list[[j]] <- res
  
  res_sub <- results[[j]] %>%
    as.data.frame() %>%
    filter(padj <= 0.05) %>%
    rownames_to_column("gene")
  
  results_sig[[j]] <- res_sub
}


base::saveRDS(res_list, file = here::here("data/inflammatory_cluster/files/deseq2.rds"))
base::saveRDS(results_sig, file = here::here("data/inflammatory_cluster/files/significant_deseq2.rds"))


# Plots -------------------------------------------------------------------
## PCA ----

pca_list <- list()

for (j in names(dds_list)) {
  # Get the rlog transformed data
  rld <- get_rld(dds_list[[j]])
  
  # Get the PCA data and the percent variance explained by each PC
  pca_data <- DESeq2::plotPCA(rld, intgroup = "cluster", returnData = TRUE)
  percentVar <- base::round(100 * base::attr(pca_data, "percentVar"))  # Get % variance explained
  
  # Add sample names as a column for labeling
  pca_data$sample <- colnames(rld)
  
  # Create the PCA plot using ggplot
  pca <- ggplot2::ggplot(pca_data, aes(x = PC1, y = PC2, color = cluster)) +
    ggplot2::geom_point(size = 0.5) +
    ggplot2::geom_text(aes(label = sample), size = 1, vjust = -2, check_overlap = FALSE) +  # Add sample names as labels
    ggplot2::labs(
      title = base::paste0(j),
      x = base::paste0("PC1: ", percentVar[1], "% variance"),  # X-axis with % variance
      y = base::paste0("PC2: ", percentVar[2], "% variance")   # Y-axis with % variance
    ) +
    ggplot2::scale_color_manual(values = c("low" = "grey", "high" = "#FA8231")) +
    my_theme() +
    ggplot2::theme(legend.position = "none")
  
  # Store the PCA plot for this dataset in the list
  pca_list[[j]] <- pca
}

# Save the plots to a PDF file
pdf(here::here("data/inflammatory_cluster/figures/pca_plots.pdf"),
    width = 6,
    height = 1)

# Use patchwork to arrange and save all plots
print(patchwork::wrap_plots(pca_list) +
        plot_layout(nrow = 1))  # Adjust the layout if needed

dev.off() 

## volcano ----
# load data
volcano_list <- list()

for (j in names(results)) {
  res <- volcano_2(df = results[[j]],
                   padj_thres = 0.05,
                   logfc_thres = 0,
                   title = j)
  
  volcano_list[[j]] <- res
}

# Save the plots to a PDF file
pdf(here::here("data/inflammatory_cluster/figures/volcano_plots.pdf"),
    width = 12,  
    height = 2)  

# Use patchwork to arrange and save all plots
print(wrap_plots(volcano_list, nrow = 1))

dev.off()  # Ensure the PDF device is closed

# gene examples -----------------------------------------------------------
genes <- c("Stat1", "Stat2", "Irf9", "Irf2", "Irf1",
           "H2-K1", "H2-K2", "H2-Q4", "H2-Q6", "H2-Q7",
           "Mx1", "Mx2", "Oasl2", "Ifit1",
           "Cxcl10", "Cd274",
           "Mafa", "Nkx6-1", "Slc2a2", "Slc30a8", "Chga")

genes %in% BiocGenerics::rownames(results[["HFD_3_high_vs_HFD_3_low"]])

norm_counts <- DESeq2::counts(dds_list[["HFD_3_high_vs_HFD_3_low"]], normalized = TRUE)
  
norm_counts %>%
  BiocGenerics::as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  tidyr::pivot_longer(-gene, names_to = "condition", values_to = "exp") %>%
  dplyr::filter(gene %in% genes) %>%
  dplyr::mutate(
    Inflammatory_response = case_when(
      condition %in% c("HFD_3_R1_high", "HFD_3_R2_high") ~ "High",
      condition %in% c("HFD_3_R1_low", "HFD_3_R2_low") ~ "Low"
    ),
    Inflammatory_response = factor(Inflammatory_response, levels = c("Low", "High")),
    gene = factor(gene, levels = genes)
  ) %>%
  ggplot2::ggplot(aes(x = Inflammatory_response, y = exp, fill = Inflammatory_response)) +
  ggplot2::geom_bar(stat = "summary", fun = mean) +
  ggplot2::geom_point(color = "black") +
  ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 7)) +
  ggplot2::facet_wrap( ~ gene, scales = "free", ncol = 10) +
  ggplot2::scale_fill_manual(values = c("Low" = "grey", "High" = "#FA8231")) +
  ggplot2::labs(x = "Inflammatory Response",
                y = "Pseudo Bulk Normalized Counts") +
  my_theme() +
  ggplot2::theme(legend.position = "none",
                 axis.text.x = element_text(angle = 90, hjust = 1))

pdf(here::here("data/inflammatory_cluster/figures/genes_examples_heatmap.pdf"),
    width = 2, 
    height = 4)
## gene examples heatmap ucell ----
myCol <- colorRampPalette(c('#004B7A', 'white', '#A83708'))(100)
myBreaks <- seq(-1.5, 1.5, length.out = 100)

norm_counts %>% 
  BiocGenerics::as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  dplyr::filter(gene %in% genes) %>%
  dplyr::arrange(factor(gene, levels = genes)) %>%
  tibble::column_to_rownames("gene") %>%
  as.matrix() %>%
  pheatmap::pheatmap(
    cluster_cols = TRUE,
    cluster_rows = FALSE,
    scale = "row",
    border_color = "black",
    color = myCol,
    breaks = myBreaks)
dev.off()
## volcano specific genes ----
genes_list_up <- list("HFD_3_high_vs_HFD_3_low" = c("Nfkb1", "Stat1", "Stat2", "Irf9", "Irf2", "Irf1",
                                                    "H2-K1", "H2-K2", "H2-Q4", "H2-Q6", "H2-Q7",
                                                    "Mx1", "Mx2", "Oasl2", "Ifit1",
                                                    "Cxcl10", "Cd274"),
                      "HFD_3_low_vs_HFD_1_low" = c("Slc2a2", "Nkx6-1", "Pax6", "Glplr", "Glul",
                                                   "Foxa2", "Rfx6", "Vsnl1", "Lrp5", "Stx1a"))

genes_list_down <- list("HFD_3_high_vs_HFD_3_low" = c("Slc2a2", "Nkx6-1", "Mafa", "Slc30a8",
                                                      "Chga", "Ttr", "Syt13"),
                        "HFD_3_low_vs_HFD_1_low" = c("Hmgcr", "Mvd", "Fdft1", "Sqle", "Srebf2",
                                                     "Insig1", "Cyp51"))

volcano_list <- list()

for (j in c("HFD_3_high_vs_HFD_3_low", "HFD_3_low_vs_HFD_1_low")) {
    
    res <- volcano_genes(df = results[[j]],
                         padj_thres = 0.05,
                         logfc_thres = 0,
                         title = paste0(j),
                         genes_up = genes_list_up[[j]],
                         genes_down = genes_list_down[[j]])
    
    volcano_list[[j]] <- res
    
    print(res)
  }

pdf(here::here("data/inflammatory_cluster/figures/volcano_plots_gene_examples.pdf"),
    width = 4,  
    height = 2)  

# Use patchwork to arrange and save all plots
print(wrap_plots(volcano_list, nrow = 1))

dev.off()  # Ensure the PDF device is closed



# module score for public data --------------------------------------------
# load data
sig_public <- readRDS(here::here("data/public_data/glu_lip_cyt_human/GSE218316/processed/degs_sig.rds"))

# convert to mouse gene symbols
mouse_genes <- purrr::map(sig_public, ~{.x %>% dplyr::pull("gene") %>% human2mouse_symbol()})

# add module score
seurat_beta <- seurat_beta %>%
  UCell::AddModuleScore_UCell(features = mouse_genes,
                              name = "_ucell",
                              maxRank = 1000,
                              BPPARAM = MulticoreParam(workers = parallel::detectCores() - 1))


# plot --------------------------------------------------------------------
arrangement <- c("IFNa_24h", "IL1b.IFNg_24h", "IFNa_72h", "IL1b.IFNg_72h")

data <- seurat_beta@meta.data %>% 
  BiocGenerics::as.data.frame() %>% 
  dplyr::select(dplyr::starts_with(arrangement), condition, is_inflammatory_ucell) %>% 
  tidyr::pivot_longer(c(-condition, -is_inflammatory_ucell), names_to = "ucell", values_to = "score") %>% 
  dplyr::mutate(condition = base::factor(condition, levels = condition_levels),
                ucell = base::factor(ucell, levels = paste0(arrangement, "_ucell")))

med <- data%>% 
  dplyr::group_by(ucell) %>% 
  dplyr::filter(condition == "LFD") %>% 
  dplyr::summarise(median =  stats::median(score))

p <- data %>% 
  ggplot2::ggplot(aes(x = condition, y = score, color = is_inflammatory_ucell)) +
  ggrastr::geom_jitter_rast(width = 0.3, size = 0.2, stroke = 0, alpha = 0.8) +
  ggplot2::geom_hline(data = med, aes(yintercept = median), colour = 'red') +
  ggplot2::scale_color_manual(values = c("low" = "grey", "high" = "#FA8231")) +
  ggplot2::geom_boxplot(outlier.shape = NA,
                        width = 0.5, fill = "transparent", color = "black") +
  ggplot2::expand_limits(y = 0) +
  ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
  ggplot2::facet_wrap(~ ucell, scales = "free", nrow = 1) +
  my_theme() +
  ggplot2::theme(legend.position = "none",
                 axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))
  
ggplot2::ggsave(here::here("data/inflammatory_cluster/figures/cytokine_public_data_module_score.pdf"),
                p,
                width = 4,
                height = 1.5)  

# heatmap of eRegulons -----------------------------------------------------

# Chromvar deviations
dev_motifs_beta <- readRDS(here::here("data/scenicplus/beta/results/chromvar_deviations_motifs_beta.rds"))

# eregu_motifs
eregu_motifs <- readRDS(here::here("data/scenicplus/beta/results/eregulon_motifs.rds"))

# high quality erequlons
high_quality_eregulon <- readRDS(here::here("data/scenicplus/beta/results/high_quality_eregulons.rds"))

# Extract data ------------------------------------------------------------

# Transcription factor gene expression
tfs <- unique(eregu_motifs$tf)

# Extract normalized counts for each TF in each cell
rna_exp_tfs <- seurat_beta %>%
  Seurat::GetAssayData(assay = "RNA", slot = "data") %>%
  BiocGenerics::as.data.frame() %>%
  t() %>%
  BiocGenerics::as.data.frame() %>%
  dplyr::select(all_of(tfs)) %>%
  t() %>%
  BiocGenerics::as.data.frame()

# Deviations for each TF per cell
eregu_list <- eregu_motifs %>%
  base::split(.$tf) %>%
  purrr::map(~ dplyr::select(., motifs) %>% unlist())

# Create list of motifs associated with each TF
eregu_motifs_list <- purrr::map(eregu_list, ~ {
  unlist(lapply(., function(motif) {
    rownames(dev_motifs_beta)[grep(motif, rownames(dev_motifs_beta))]
  }))
})

# Calculate mean motif deviations per eregulon for each cell
mean_activity_eregu <- purrr::map2_dfr(
  eregu_motifs_list,
  names(eregu_motifs_list),
  ~ {
    dev_motifs_beta %>%
      BiocGenerics::as.data.frame() %>%
      tibble::rownames_to_column("motif") %>%
      dplyr::filter(motif %in% .x) %>%
      tidyr::pivot_longer(-motif, names_to = "barcode", values_to = "values") %>%
      dplyr::group_by(barcode) %>%
      dplyr::summarise(mean = mean(values), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = barcode, values_from = mean) %>%
      dplyr::mutate(tf = .y) %>%
      tibble::column_to_rownames("tf")
  }
)

# Make rownames to columns
tf_scale <- rna_exp_tfs %>%
  rownames_to_column("tf")

dev_scale <- mean_activity_eregu %>%
  rownames_to_column("tf")

# Combine data: TF expression and deviation score into 1 dataframe
merged_df <- dplyr::full_join(
  tidyr::pivot_longer(tf_scale, -tf, names_to = "cell", values_to = "tf_exp"),
  tidyr::pivot_longer(dev_scale, -tf, names_to = "cell", values_to = "tf_dev"),
  by = c("tf", "cell")
) %>%
  tidyr::drop_na() %>%
  dplyr::mutate(condition = factor(str_extract(cell, "LFD|HFD_1|HFD_3"), levels = condition_levels))


## Heatmap
# Select the appropriate column from seurat_beta@meta.data
  sub_df <- seurat_beta@meta.data %>%
    BiocGenerics::as.data.frame() %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::select(cell, is_inflammatory_ucell)
  
  # Final merged dataframe with annotations
  merged_df_anno <- merged_df %>%
    dplyr::left_join(sub_df, by = "cell")
  
  # Calculate mean expression and deviation for each TF
  mean_merged_df <- merged_df_anno %>%
    dplyr::group_by(tf, condition, is_inflammatory_ucell) %>%  
    dplyr::summarise(
      mean_exp = mean(expm1(tf_exp)),  # Reverse log transformation
      mean_dev = mean(tf_dev),
      .groups = 'drop'  # Avoids keeping grouping structure
    ) %>%
    dplyr::mutate(id = paste(condition, is_inflammatory_ucell, sep = "_")) %>%
    dplyr::filter(id != "LFD_high")  # Don't include the few cells from LFD
  
  # Define color palette and breaks
  myCol <- colorRampPalette(c('#004B7A', 'white', '#A83708'))(100)
  myBreaks <- seq(-1.5, 1.5, length.out = 100)
  
  title <- stringr::str_remove(i, "is_")
  
  # Save heatmap as PDF
  pdf(here::here("data/inflammatory_cluster/figures/tf_expression_activity_heatmap.pdf"),
      width = 2,
      height = 2)
  
  # Gene expression heatmap
  print(mean_merged_df %>%
          dplyr::ungroup() %>%
          dplyr::select(tf, id, mean_exp) %>%
          dplyr::filter(tf %in% c("Nfkb1", "Stat2", "Stat1", "Irf2", "Irf7", "Irf9")) %>%
          dplyr::arrange(factor(tf, levels = c("Nfkb1", "Stat2", "Stat1", "Irf2", "Irf7", "Irf9"))) %>%
          dplyr::arrange(factor(id, levels = c("LFD_low", "HFD_1_low", "HFD_3_low", "HFD_1_high", "HFD_3_high"))) %>%
          tidyr::pivot_wider(id_cols = tf, names_from = id, values_from = mean_exp) %>%
          tibble::column_to_rownames("tf") %>%
          pheatmap::pheatmap(main = "Gene expression",
                             cluster_cols = FALSE,
                             cluster_rows = FALSE,
                             scale = "row",
                             border_color = "black",
                             color = myCol,
                             breaks = myBreaks))
  
  # Motif activity heatmap
  print(mean_merged_df %>%
          dplyr::ungroup() %>%
          dplyr::select(tf, id, mean_dev) %>%
          dplyr::filter(tf %in% c("Nfkb1", "Stat2", "Stat1", "Irf2", "Irf7", "Irf9")) %>%
          dplyr::arrange(factor(tf, levels = c("Nfkb1", "Stat2", "Stat1", "Irf2", "Irf7", "Irf9"))) %>%
          dplyr::arrange(factor(id, levels = c("LFD_low", "HFD_1_low", "HFD_3_low", "HFD_1_high", "HFD_3_high"))) %>%
          tidyr::pivot_wider(id_cols = tf, names_from = id, values_from = mean_dev) %>%
          tibble::column_to_rownames("tf") %>%
          pheatmap::pheatmap(main = "Motif activity",
                             cluster_cols = FALSE,
                             cluster_rows = FALSE,
                             scale = "row",
                             border_color = "black",
                             color = myCol,
                             breaks = myBreaks))
  
  dev.off()

# get normalized counts my data -------------------------------------------
seurat_beta <- readRDS(here::here("data/inflammatory_cluster/files/seurat_beta.rds"))

# Pseudobulk counts per replicate per condition
counts <- seurat_beta %>%
  edgeR::Seurat2PB(sample = "orig.ident",
                   cluster = "is_inflammatory_ucell")

counts_2 <- counts %>%
  as.data.frame() %>%
  dplyr::rename_with(~gsub("cluster", "", .x)) %>%
  tibble::column_to_rownames("gene")

# Get meta data
meta_data <- counts[["samples"]] %>%
  magrittr::set_rownames(base::gsub("cluster", "", rownames(.)))

# check that colnames in counts_2 and rownames in meta_data are the same order
print(all.equal(colnames(counts_2), rownames(meta_data)))

dds_seu <-
  DESeqDataSetFromMatrix(
    countData = counts_2,
    colData =  meta_data,
    design = ~ sample + cluster)

dds_seu <- DESeq(dds_seu)

norm_counts_seurat <- DESeq2::counts(dds_seu, normalized = TRUE) %>%
  BiocGenerics::as.data.frame() %>%
  tibble::rownames_to_column("gene")

# save normalized counts
saveRDS(norm_counts_seurat, here::here("data/inflammatory_cluster/files/pseudo_norm_count.rds"))

# Gene examples UPR response ----------------------------------------------
upr_genes <- c("Atf6", "Ern1", "Eif2ak3")

p <- print(norm_counts_seurat %>%
        BiocGenerics::as.data.frame() %>%
        dplyr::select(gene, ends_with("_Low")) %>%
        tidyr::pivot_longer(-gene, names_to = "condition", values_to = "exp") %>%
        dplyr::filter(gene %in% upr_genes) %>%
        dplyr::mutate(Inflammatory_response = case_when(condition %in% c("HFD_3_R1_low", "HFD_3_R2_low") ~ "HFD_3",
                                                        condition %in% c("HFD_1_R1_low", "HFD_1_R2_low") ~ "HFD_1",
                                                        condition %in% c("LFD_R1_low", "LFD_R2_low", "LFD_R3_low") ~ "LFD"),
                      Inflammatory_response = factor(Inflammatory_response, levels = c("LFD", "HFD_1", "HFD_3"))) %>%
        ggplot2::ggplot(aes(x = Inflammatory_response, y = exp, fill = Inflammatory_response)) +
        ggplot2::geom_bar(stat="summary", fun=mean) +
        ggplot2::geom_point(color = "black", size = 0.1)+
        ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
        ggplot2::facet_wrap(~ gene, scales = "free", ncol = 10) +
        ggplot2::scale_fill_manual(values = c("LFD" = "#004B7A", "HFD_1" = "#c71500", "HFD_3" = "#fcb21c")) +
        ggplot2::labs(x = "Inflammatory Response",
                      y = "Pseudo Bulk Normalized Counts") +
        my_theme() +
        ggplot2::theme(legend.position = "none",
                       axis.text.x = element_text(angle = 90, hjust = 1)))

ggplot2::ggsave(here::here("data/inflammatory_cluster/figures/inflammatory_ucell_gene_examples_HFD1_vs_LFD_low_UPR.pdf"),
                p, width = 2, height = 1)

# Define UPR and INS gene sets ------------------------------------------------
# Source: https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/GOBP_ATF6_MEDIATED_UNFOLDED_PROTEIN_RESPONSE.html

# ATF6-mediated unfolded protein response genes
atf6 <- dplyr::bind_rows(msigdbr(species = "Mus musculus", category = "C5")) %>%
  dplyr::filter(grepl("GOBP_ATF6_MEDIATED_UNFOLDED_PROTEIN_RESPONSE", gs_name)) %>%
  dplyr::pull("gene_symbol")

# IRE1-mediated unfolded protein response genes
ire1 <- dplyr::bind_rows(msigdbr(species = "Mus musculus", category = "C5")) %>%
  dplyr::filter(grepl("GOBP_IRE1_MEDIATED_UNFOLDED_PROTEIN_RESPONSE", gs_name)) %>%
  dplyr::pull("gene_symbol")

# PERK-mediated unfolded protein response genes
perk <- dplyr::bind_rows(msigdbr(species = "Mus musculus", category = "C5")) %>%
  dplyr::filter(grepl("GOBP_PERK_MEDIATED_UNFOLDED_PROTEIN_RESPONSE", gs_name)) %>%
  dplyr::pull("gene_symbol")

# INS gene set
# Source: https://pmc.ncbi.nlm.nih.gov/articles/PMC8259463/#R87
ins <- c("Pdx1","Nkx6-1", "Mafa", "Pax6", "Insm1", "Neurod1", "Nkx2-2", "Foxo1", # identity and maturity
         "Slc2a2", # glucose sensing
         "Slc30a8", # insulin granule formation and exocytosis
         "Chga", # insulin granule cargo
         "Vsnl1", # calcium dependent insulin secretion
         "Kcnj11", # Potassium transporters
         "Foxa2", "Rfx6", "Pcsk1", # insulin expression and processing
         "Glp1r", "Hadh", # amplification of insulin secretion
         "Hnf4a", "Foxo3", "Foxo4", "Crtc2", "Nfatc1", "Nfatc2", "Syt13", "Ttr")

# Combine gene sets into a list
genes_by_pathway <- list("upr" = unique(c(atf6, ire1, perk)), "ins" = ins)

# Add UCell scores to Seurat object -------------------------------------------
seurat_beta <- seurat_beta %>% UCell::AddModuleScore_UCell(features = genes_by_pathway,
                                                           name = "_ucell",
                                                           maxRank = 1000,
                                                           BPPARAM = MulticoreParam(workers = parallel::detectCores() - 1))

# Calculate upper quantile of UCell scores for LFD condition ------------------
lfd_upper_q <- seurat_beta@meta.data %>%
  dplyr::mutate(Inflammatory_response = case_when(condition == "HFD_3" & is_inflammatory_ucell == "high" ~ "HFD_3_high",
                                                  condition == "HFD_3" & is_inflammatory_ucell == "low" ~ "HFD_3_low",
                                                  TRUE ~ condition),
                Inflammatory_response = factor(Inflammatory_response, levels = c("LFD", "HFD_1", "HFD_3_low", "HFD_3_high"))) %>%
  dplyr::filter(Inflammatory_response == "LFD") %>%
  dplyr::select(orig.ident, names(genes_by_pathway), Inflammatory_response) %>%
  tidyr::pivot_longer(c(-Inflammatory_response, -orig.ident), names_to = "gene", values_to = "score") %>%
  dplyr::group_by(gene, Inflammatory_response) %>%
  dplyr::summarize(median_score = quantile(score, na.rm = TRUE, probs = 0.75)) %>%
  dplyr::mutate(Inflammatory_response = factor(Inflammatory_response, levels = c("LFD", "HFD_1", "HFD_3_low", "HFD_3_high")))

# Prepare data for plotting ---------------------------------------------------
df <- seurat_beta@meta.data %>%
  dplyr::mututate(Inflammatory_response = case_when(condition == "HFD_3" & is_inflammatory_ucell == "high" ~ "HFD_3_high",
                                                    condition == "HFD_3" & is_inflammatory_ucell == "low" ~ "HFD_3_low",
                                                    TRUE ~ condition),
                  Inflammatory_response = factor(Inflammatory_response, levels = c("LFD", "HFD_1", "HFD_3_low", "HFD_3_high"))) %>%
  dplyr::select(names(genes_by_pathway), Inflammatory_response, is_inflammatory_ucell) %>%
  tidyr::pivot_longer(cols = c(-Inflammatory_response, -is_inflammatory_ucell), names_to = "gene", values_to = "score") %>%
  dplyr::mutate(Inflammatory_response = factor(Inflammatory_response, levels = c("LFD", "HFD_1", "HFD_3_low", "HFD_3_high")))

# Plot UCell scores -----------------------------------------------------------
p <- df %>%
  ggplot(aes(x = Inflammatory_response, y = score)) +
  ggrastr::geom_jitter_rast(color = "grey", width = 0.3, size = 0.1, stroke = 0, alpha = 0.8) +
  geom_hline(data = lfd_upper_q, aes(yintercept = median_score), colour = 'red') +
  geom_boxplot(outlier.shape = NA, width = 0.6, fill = "transparent") +
  expand_limits(y = 0) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
  labs(x = "", y = "activity score") +
  facet_wrap(~gene, scale = "free", ncol = 1) +
  my_theme() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# Save the plot to a PDF file -------------------------------------------------
ggsave(here::here("data/inflammatory_cluster/figures/upr_ins_auc_ucell_score.pdf"),
       p, width = 1, height = 3)

# Save Seurat metadata and gene pathways --------------------------------------
saveRDS(seurat_beta@meta.data %>% tibble::rownames_to_column("barcode"), here::here("data/inflammatory_cluster/files/seurat_beta_perk_ire1_atf6.rds"))
saveRDS(genes_by_pathway, here::here("data/inflammatory_cluster/files/genes_by_pathway_perk_ire1_atf6.rds"))


# get normalized counts my data -------------------------------------------
seurat_beta@meta.data <- seurat_beta@meta.data %>%
  dplyr::mutate(cond_inf = case_when(condition == "HFD_3" & is_inflammatory_ucell == "high" ~ "high",
                                     condition == "HFD_3" & is_inflammatory_ucell == "low" ~ "low",
                                     .default = ""))

# Pseudobulk counts per replicate per condition
counts <- seurat_beta %>%
  edgeR::Seurat2PB(sample = "orig.ident",
                   cluster = "cond_inf")

counts_2 <- counts %>%
  as.data.frame() %>%
  dplyr::rename_with(~gsub("cluster", "", .x)) %>%
  tibble::column_to_rownames("gene")

# Get meta data
meta_data <- counts[["samples"]] %>%
  magrittr::set_rownames(base::gsub("cluster", "", rownames(.)))

# check that colnames in counts_2 and rownames in meta_data are the same order
print(all.equal(colnames(counts_2), rownames(meta_data)))

dds_seu <-
  DESeqDataSetFromMatrix(
    countData = counts_2,
    colData =  meta_data,
    design = ~ cluster)

dds_seu <- DESeq(dds_seu)

norm_counts_seurat <- DESeq2::counts(dds_seu, normalized = TRUE) %>%
  BiocGenerics::as.data.frame() %>%
  tibble::rownames_to_column("gene")

# Plot normalized counts for selected genes -----------------------------------

genes = c("Atf6", "Ern1", "Eif2ak3", "Mafa", "Slc2a2", "Nkx6-1")
p <- norm_counts_seurat %>%
  tidyr::pivot_longer(-gene, names_to = "condition", values_to = "exp") %>%
  dplyr::filter(gene %in% genes) %>%
  dplyr::mutate(Inflammatory_response = case_when(condition %in% c("HFD_3_R1_low", "HFD_3_R2_low") ~ "HFD_3_low",
                                                  condition %in% c("HFD_3_R1_high", "HFD_3_R2_high") ~ "HFD_3_high",
                                                  condition %in% c("HFD_1_R1_", "HFD_1_R2_") ~ "HFD_1",
                                                  condition %in% c("LFD_R1_", "LFD_R2_", "LFD_R3_") ~ "LFD"),
                Inflammatory_response = factor(Inflammatory_response, levels = c("LFD", "HFD_1", "HFD_3_low", "HFD_3_high")),
                gene = factor(gene, levels = genes)) %>%
  ggplot2::ggplot(aes(x = Inflammatory_response, y = exp, fill = Inflammatory_response)) +
  ggplot2::geom_bar(stat="summary", fun=mean) +
  ggplot2::geom_point(color = "black", size = 0.2)+
  ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
  ggplot2::facet_wrap(~ gene, scales = "free_y", nrow = 2) +
  ggplot2::scale_fill_manual(values = c("LFD" = "#004B7A", "HFD_1" = "#c71500", "HFD_3_low" = "grey", "HFD_3_high" = "orange")) +
  ggplot2::labs(x = "Inflammatory Response",
                y = "Pseudo Bulk Normalized Counts") +
  my_theme() +
  ggplot2::theme(legend.position = "none",
                 axis.text.x = element_text(angle = 90, hjust = 1))

ggplot2::ggsave(here::here("data/inflammatory_cluster/figures/upr_ins_gene_examples.pdf"),
                p, width = 2, height = 2)


# Number of cells with a UPR or INS score above or below the median --------
## Get the 75th percentile of UPR and INS scores in the LFD group ----
lfd_upper_q <- seurat_beta@meta.data %>%
  dplyr::mutate(Inflammatory_response = case_when(condition == "HFD_3" & is_inflammatory_ucell == "high" ~ "HFD_3_high",
                                                  condition == "HFD_3" & is_inflammatory_ucell == "low" ~ "HFD_3_low",
                                                  .default = condition),
                Inflammatory_response = factor(Inflammatory_response, levels = c("LFD", "HFD_1", "HFD_3_low", "HFD_3_high"))) %>%
  dplyr::filter(Inflammatory_response == "LFD") %>%
  dplyr::select(orig.ident, upr, ins, Inflammatory_response) %>%
  tidyr::pivot_longer(c(-Inflammatory_response, -orig.ident), names_to = "gene", values_to = "score") %>%
  dplyr::group_by(gene, Inflammatory_response) %>%
  dplyr::summarize(uppr_quartile = stats::quantile(score, na.rm = TRUE, probs = 0.75)) %>%
  dplyr::select(-Inflammatory_response) %>%
  tibble::deframe()

## Get number of cells which have an inflammatory score above the 75th quantile of the LFD group ----
n_df <- seurat_beta@meta.data %>%
  dplyr::mutate(Inflammatory_response = case_when(condition == "HFD_3" & is_inflammatory_ucell == "high" ~ "HFD_3_high",
                                                  condition == "HFD_3" & is_inflammatory_ucell == "low" ~ "HFD_3_low",
                                                  .default = condition),
                Inflammatory_response = factor(Inflammatory_response, levels = c("LFD", "HFD_1", "HFD_3_low", "HFD_3_high"))) %>%
  dplyr::select(orig.ident, upr, ins, Inflammatory_response) %>%
  dplyr::group_by(orig.ident, Inflammatory_response) %>%
  dplyr::summarise(n_above_upr = (sum(upr > lfd_upper_q["upr_ucell"])),
                   n_above_ins = (sum(ins > lfd_upper_q["ins_ucell"])),
                   n_total = n(),
                   frac_upr = n_above_upr/n_total,
                   frac_ins = n_above_ins/n_total)

# Plot fractions of high UPR and INS scores -----------------------------------
p <- n_df %>%
  ggplot2::ggplot(aes(x = Inflammatory_response, y = frac_upr, fill = Inflammatory_response)) +
  ggplot2::geom_bar(stat="summary", fun=mean) +
  ggplot2::geom_point(size = 0.2) +
  ggplot2::labs(y = "Fraction of high UPR") +
  n_df %>%
  ggplot2::ggplot(aes(x = Inflammatory_response, y = frac_ins, fill = Inflammatory_response)) +
  ggplot2::geom_bar(stat="summary", fun=mean) +
  ggplot2::geom_point(size = 0.2) +
  ggplot2::labs(y = "Fraction of high function") &
  ggplot2::scale_fill_manual(values = c("LFD" = "#004B7A",
                                        "HFD_1" = "#c71500",
                                        "HFD_3_low" = "grey",
                                        "HFD_3_high" = "#FA8231")) &
  my_theme() &
  ggplot2::theme(legend.position = "none",
                 axis.text.x = element_text(angle = 90, hjust = 1),
                 axis.title.x = element_blank())
  
ggplot2::ggsave(here::here("data/inflammatory_cluster/figures/upr_ins_auc_ucell_fraction.pdf"),
                p, width = 1, height = 2)

pdf(here::here(paste0("data/inflammatory_cluster/figures/upr_ins_auc_ucell_fraction_taller.pdf")),
    height = 5 / 2.54,
    width = 5 / 2.54)

# Prepare data for ECDF plots -------------------------------------------------
df <- seurat_beta@meta.data %>%
  dplyr::mutate(Inflammatory_response = case_when(
    condition == "HFD_3" & is_inflammatory_ucell == "high" ~ "HFD_3_high",
    condition == "HFD_3" & is_inflammatory_ucell == "low" ~ "HFD_3_low",
    TRUE ~ condition
  ),
  Inflammatory_response = factor(Inflammatory_response, levels = c("LFD", "HFD_1", "HFD_3_low", "HFD_3_high"))) %>%
  dplyr::select(orig.ident, upr, ins, Inflammatory_response)

## ECDF plots ----
# First plot (ECDF plots)
ecdf_plot_ins <- ggplot2::ggplot(df, aes(x = ins, color = Inflammatory_response)) +
  ggplot2::stat_ecdf(size = 1) +
  ggplot2::geom_vline(xintercept = lfd_upper_q["ins"], linetype = "dashed") +
  ggplot2::labs(x = "Score", y = "% of cells") +
  ggplot2::scale_y_reverse() +
  ggplot2::scale_color_manual(values = c("LFD" = "#004B7A", "HFD_1" = "#c71500", "HFD_3_low" = "grey", "HFD_3_high" = "#FA8231")) +
  my_theme() +
  ggplot2::theme(legend.position = "none")

# Second plot (ECDF plots)
ecdf_plot_upr <- ggplot(df, aes(x = upr, color = Inflammatory_response)) +
  ggplot2::stat_ecdf(size = 1) +
  ggplot2::geom_vline(xintercept = lfd_upper_q["upr"], linetype = "dashed") +
  ggplot2::labs(x = "Score", y = "% of cells") +
  ggplot2::scale_y_reverse() +
  ggplot2::scale_color_manual(values = c("LFD" = "#004B7A", "HFD_1" = "#c71500", "HFD_3_low" = "grey", "HFD_3_high" = "#FA8231")) +
  my_theme() +
  ggplot2::theme(legend.position = "none")

ggsave(here::here("data/inflammatory_cluster/figures/ins__upr_ucell_cumulative_perc.pdf"),
       plot = gridExtra::grid.arrange(ecdf_plot_ins, ecdf_plot_upr, ncol = 2),
       width = 2, height = 1)


# Save seurat -------------------------------------------------------------
saveRDS(seurat_beta, here::here("data/inflammatory_cluster/files/seurat_beta.rds"))
