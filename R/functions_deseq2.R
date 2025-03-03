#' EdgeRs Seurat2PB code to create pseubulk counts, for different assays
#'
#' @param object a seurat object
#' @param sample column in seurat object metadata, that contain samples - e.g. orig.ident
#' @param assay the assay from which to get counts - e.g. "RNA"
#' @param cluster column in seurat object metadata, that contain cluster annotation - e.g. "seuratc_clusters"
#'
#' @return a dataframe with pseudo counts
pseudobulk <- function (object, sample, assay, cluster = "seurat_clusters")
{
    if (!is(object, "SeuratObject") & !is(object, "Seurat"))
        stop("object is not of the SeuratObject or Seurat class")
    if (!requireNamespace("SeuratObject", quietly = TRUE))
        stop("SeuratObject package required but is not installed (or can't be loaded)")

    counts <- SeuratObject::GetAssayData(object, assay = paste0(assay),
                                         slot = "counts")
    if (is.null(counts))
        stop("object doesn't contain raw RNA counts")
    meta <- object@meta.data
    if (!sample %in% names(meta))
        stop("sample information can not be found in meta.data")
    if (!cluster %in% names(meta))
        stop("cluster information can not be found in meta.data")
    sp <- meta[, sample]
    clst <- meta[, cluster]
    if (length(table(sp)) == 1)
        warning("Only 1 sample found in meta.data. Please check whether sample information is specified correctly.")
    if (length(table(clst)) == 1)
        warning("Only 1 cluster found in meta.data. Please check whether cluster information is specified correctly.")
    genes <- data.frame(gene = rownames(object[[paste0(assay)]]))
    genes <- cbind(genes, object[[paste0(assay)]][[]])
    sp_clst <- factor(paste(sp, clst, sep = "_cluster"))
    group_mat <- Matrix::sparse.model.matrix(~0 + sp_clst)
    colnames(group_mat) <- gsub("^sp_clst", "", colnames(group_mat))
    counts.pb <- counts %*% group_mat
    levels(sp_clst)
    sp.pb <- gsub("_cluster.*$", "", levels(sp_clst))
    clst.pb <- gsub("^.*_cluster", "", levels(sp_clst))
    sample.pb <- data.frame(sample = sp.pb, cluster = clst.pb)
    edgeR::DGEList(counts = as.matrix(counts.pb), samples = sample.pb,
            genes = genes)
}

#' Create pseudobulk count list
#'
#' @param object a seurat object
#' @param sample a vector of samples
#' @param cluster a vector of clusters
#' @param assay assay to get data from
#' @param sample_levels a vector of samples with factor levels
#' @param anno_levels a vector of cluster annotation factor levels
#'
#' @return a list of dataframes with pseudobulk counts divided into clusters per sample
pseudobulk_list <- function(object, sample, cluster, assay, sample_levels, anno_levels){
    output <- object %>%

        # pseudocounts
        pseudobulk(sample=paste0(sample), cluster=paste0(cluster), assay = paste0(assay)) %>%

        # remove cluster from column names
        as.data.frame() %>%
        dplyr::rename_with(~gsub("cluster", "", .x)) %>%

        # make dataframe into long format
        tidyr::pivot_longer(-gene, names_to = "sample", values_to = "count") %>%

        # create cluster which contains sample and cluster information, by splitting sample column
        tidyr::separate(sample,  into = c("sample", "cluster"), sep="_(?=[^_]+$)") %>%

        # sample and cluster as factor
        dplyr::mutate(sample = factor(sample, levels = sample_levels),
                      cluster = factor(cluster, levels = anno_levels)) %>%

        # split by cluster and create list of pseudobulk counts for each cluster
        dplyr::group_split(cluster, .keep = FALSE) %>%

        # make each dataframe wider and set gene column to rownames
        purrr::modify_depth(1, ~ tidyr::pivot_wider(., names_from = sample, values_from = count, names_sort = TRUE) %>%
                                tibble::column_to_rownames("gene"))

    # add names to each dataframe
    names(output) <- anno_levels

    return(output)
}

#' Create a list of metadata from seurat object
#'
#' @param object a seurat object
#' @param condition_levels a vector of condition factor levels
#' @param sample_levels a vector of sample factor levels
#' @param anno_levels a vector of cluster sample factor levels
#'
#' @return a list of dataframes with metadata divided into clusters
meta_list <- function(object, condition_levels, sample_levels, anno_levels){
    output <- object@meta.data %>%
        as.data.frame() %>%

        # select specific columns to keep
        dplyr::select(orig.ident, condition, manual_anno) %>%

        # get numbers of cells in each cluster + replicate
        dplyr::group_by(orig.ident, condition, manual_anno) %>%
        dplyr::summarise(cell_number = n()) %>%

        # set factor levels and create a batch column
        dplyr::mutate(orig.ident = factor(orig.ident, levels = sample_levels),
                      condition = factor(condition, levels = condition_levels),
                      manual_anno = factor(manual_anno, levels = anno_levels)) %>%
        dplyr::ungroup() %>%

        # split dataframe by cluster
        dplyr::group_split(manual_anno) %>%

        # convert orig.ident to rownames
        purrr::modify_depth(1, ~ tibble::column_to_rownames(., "orig.ident"))

    names(output) <- anno_levels
    return(output)
}

#' Filter counts and create dds object
#'
#' @param counts a dataframe of counts
#' @param meta a meta dataframe
#' @param design the design model
#'
#' @return a dds object
create_dds <- function(counts, meta, design){

    # remove counts that are 0
    counts_filt <- counts %>%
        dplyr::filter(base::rowSums(across(where(is.numeric)))>0)

    # create dds object
    output <- DESeq2::DESeqDataSetFromMatrix(countData = counts_filt,
                                             colData = meta,
                                             design = design)

    return(output)
}

#' LRT dds
#'
#' @param dds a dds object
#' @param test_use a character string of what test to use
#' @param reduced the model
#'
#' @return
#' @export
#'
#' @examples
dds_test_lrt <- function(dds, test_use, reduced){
    output <- DESeq2::DESeq(dds, test=test_use, reduced = reduced)
    return(output)
}

#' extract parwise comparisons from deseq
#'
#' @param dds a dese2 object
#' @param con a contrast
#'
#' @return a dataframe
get_results_lrt <- function(dds){
    output <- DESeq2::results(dds)
    return(output)
}


#' Wald
#'
#' @param dds a dds object
#' @param test_use a character string of what test to use
#'
#' @return
#' @export
#'
#' @examples
dds_test_wald <- function(dds, test_use) {
    output <- DESeq2::DESeq(dds, test=test_use)
    return(output)
}

#' extract parwise comparisons from deseq
#'
#' @param dds a dese2 object
#' @param con a contrast
#'
#' @return a dataframe
get_results <- function(dds, con){
    output <- DESeq2::results(dds, contrast = con)
    return(output)
}

#' Extract significant genes from dataframe
#'
#' @param df a dataframe
#' @param padj_thres p-value threshold
#' @param only_padj do you want to filter only based on adjusted p-value
#' @param logfc_thres log2fc threshold
#'
#' @return a filtered dataframe
get_sig <- function(df, padj_thres, only_padj, logfc_thres){
    if (only_padj == "yes")  {
        print(paste0("Filteing with adjust pvalue <=", padj_thres))
        output <- df %>%
            as.data.frame() %>%
            dplyr::filter(padj <= padj_thres)
    }
    else{
        print(paste0("Filteing with adjust pvalue <=", padj_thres, " and abs(log2foldchange) >= ", logfc_thres))
        output <- df %>%
            as.data.frame() %>%
            dplyr::filter(padj <= padj_thres & abs(log2FoldChange) >= logfc_thres)
    }

    return(output)
}

#' log transform data with deseq2
#'
#' @param dds a deseq2 object
#'
#' @return a dds object with log transformed counts
get_rld<- function(dds){
    output <- DESeq2::rlog(dds)
    return(output)
}

#' extract log transformed values from deseq2
#'
#' @param rld a dds object with log transofrmed counts
#'
#' @return a dataframe
extract_rld <- function(rld){
    output <- as.data.frame(assay(rld))
    return(output)
}

plot_pca <- function(rld, cluster){

    # Add replicate annotation to coldata
    colData(rld)$replicate <- BiocGenerics::rownames(SummarizedExperiment::colData(rld))

    output <- DESeq2::plotPCA(rld, intgroup = c("replicate")) +
        ggplot2::labs(title = base::paste0(cluster)) +
      ggplot2::ggtitle(base::paste0(cluster)) +
        ggplot2::scale_color_manual(values = sample_color) +
      my_theme() +
      ggplot2::coord_cartesian(clip = "off") + 
      ggplot2::theme(legend.position = "none")
    return(output)
}

get_norm_counts <- function(dds){
    norm_dds <- DESeq2::estimateSizeFactors(dds)
    output <- DESeq2::counts(norm_dds,normalized = TRUE)
    return(output)
}


#' Volcano plot with labels of top 10 genes at each point
#'
#' @param df a deseq2 results dataframe
#' @param padj_thres a value e.g. 0.05
#' @param logfc_thres a value e.g. 1.5
#' @param title title of plot e.g. "HFD_3 vs LFD"
#'
#' @return a ggplot2 plot
volcano <- function(df, padj_thres, logfc_thres, title) {

    df_2 <- df %>%
        as.data.frame() %>%
        tibble::rownames_to_column("gene")


    output <-     ggplot2::ggplot(data= df_2, aes(x=log2FoldChange, y=-log10(padj))) +
        ggplot2::geom_point() +
        ggplot2::geom_point(shape=21, color="lightgrey", fill= "lightgrey", size=2) +

        # add color to significant results
        ggplot2::geom_point(data=subset(df_2, padj <= padj_thres & log2FoldChange <= -logfc_thres), shape=21, color = "black", fill ='#004B7A', size = 2) +
        ggplot2::geom_point(data=subset(df_2, padj <= padj_thres & log2FoldChange >= logfc_thres), shape=21, color = "black", fill ='#A83708', size = 2) +

        # Add number of significant results to plot
        ggplot2::geom_text(data= df_2 %>%
                               dplyr::filter(padj <= padj_thres & log2FoldChange <= -logfc_thres) %>%
                               dplyr::summarise(n = n()),
                           aes(label = n, x = -Inf, y = Inf), color = '#1465AC',
                           hjust = -0.5, vjust = 1.5) +
        ggplot2::geom_text(data= df_2 %>%
                               dplyr::filter(padj <= padj_thres & log2FoldChange >= logfc_thres) %>%
                               dplyr::summarise(n = n()),
                           aes(label = n, x = -Inf, y = Inf), color = '#B31B21',
                           hjust = -0.5, vjust = 3) +
        # add text to top 10 significant genes
        ggrepel::geom_text_repel(data= df_2 %>%
                                     dplyr::filter(padj <= padj_thres &
                                                       log2FoldChange <= -logfc_thres) %>%
                                     dplyr::arrange(padj, log2FoldChange) %>%
                                     head(10),
                                 aes(label = gene), size=3,
                                 max.overlaps = 100) +

        ggrepel::geom_text_repel(data= df_2 %>%
                                     dplyr::filter(padj <= padj_thres &
                                                       log2FoldChange >= logfc_thres) %>%
                                     dplyr::arrange(padj, log2FoldChange) %>%
                                     head(10),
                                 aes(label = gene), size=3,
                                 max.overlaps = 100) +

        # Add lines
        ggplot2::geom_vline(xintercept = c(0, logfc_thres, -logfc_thres),
                            linetype = c(1, 2, 2),
                            color = c("#868686FF", "#868686FF", "#868686FF")) +
        ggplot2::geom_hline(yintercept = -log10(padj_thres),
                            linetype = 2,
                            color = "#868686FF") +

        # add titles
        ggplot2::xlab("log2 fold change") +
        ggplot2::ylab("-log10 adjusted p-value") +
        ggplot2::ggtitle(paste0(title)) +
        ggprism::theme_prism(border = T,
                             base_fontface = "plain",
                             base_size = 12) +
        ggplot2::coord_cartesian(clip = "off")

    return(output)
}

#' Volcano plot with labels of top 10 genes aligned om each side
#'
#' @param df a deseq2 results dataframe
#' @param padj_thres a value e.g. 0.05
#' @param logfc_thres a value e.g. 1.5
#' @param title title of plot e.g. "HFD_3 vs LFD"
#'
#' @return a ggplot2 plot
volcano_2 <- function(df, padj_thres, logfc_thres, title) {

    df_2 <- df %>%
        as.data.frame() %>%
        tibble::rownames_to_column("gene")

    top <- base::max(stats::na.omit(df_2$padj))
    bottom <- base::min(stats::na.omit(df_2$padj))

    # Subset for significant genes
    sig_down <- df_2 %>%
        dplyr::filter(padj <= padj_thres & log2FoldChange <= -logfc_thres) %>%
        dplyr::arrange(padj, log2FoldChange) %>%
        head(10)

    sig_up <- df_2 %>%
        dplyr::filter(padj <= padj_thres & log2FoldChange >= logfc_thres) %>%
        dplyr::arrange(padj, log2FoldChange) %>%
        head(10)

    # Create label positions manually
    sig_down <- sig_down %>%
        mutate(label_x = min(na.omit(df_2$log2FoldChange)) - (-1.5),
               label_y = seq(-log10(bottom), -log10(top), length.out = nrow(sig_down)))

    sig_up <- sig_up %>%
        mutate(label_x = max(na.omit(df_2$log2FoldChange)) + (-1.5),
               label_y = seq(-log10(bottom), -log10(top), length.out = nrow(sig_up)))

    output <- ggplot2::ggplot(data= df_2, aes(x=log2FoldChange, y=-log10(padj))) +
        ggrastr::geom_point_rast(shape=21,color= "lightgrey", fill = "lightgray", size=0.1) +

        # Add color to significant results
        ggrastr::geom_point_rast(data=subset(df_2, padj <= padj_thres & log2FoldChange <= -logfc_thres), shape=21, color ='#004B7A', fill = '#004B7A', size = 0.5) +
        ggrastr::geom_point_rast(data=subset(df_2, padj <= padj_thres & log2FoldChange >= logfc_thres), shape=21, color ='#A83708', fill = '#A83708', size = 0.5) +

        # Add lines
        ggplot2::geom_vline(xintercept = c(0, logfc_thres, -logfc_thres),
                            linetype = c(1, 2, 2),
                            color = c("#868686FF", "#868686FF", "#868686FF")) +
        ggplot2::geom_hline(yintercept = -log10(padj_thres),
                            linetype = 2,
                            color = "#868686FF") +

        # Add number of significant results to plot
        ggplot2::geom_text(data= df_2 %>%
                               dplyr::filter(padj <= padj_thres & log2FoldChange <= -logfc_thres) %>%
                               dplyr::summarise(n = n()),
                           aes(label = paste0("n(genes) = ", n), x = -Inf, y = Inf), color = '#1465AC',
                           hjust = -0.5, vjust = 1.5, size = 2) +
        ggplot2::geom_text(data= df_2 %>%
                               dplyr::filter(padj <= padj_thres & log2FoldChange >= logfc_thres) %>%
                               dplyr::summarise(n = n()),
                           aes(label = paste0("n(genes) = ", n), x = -Inf, y = Inf), color = '#B31B21',
                           hjust = -0.5, vjust = 3, size = 2) +

        # Add text to top 10 significant genes
        ggplot2::geom_text(data=sig_down, aes(x=label_x, y=label_y, label=gene), size = 2, hjust=1) +
        ggplot2::geom_text(data=sig_up, aes(x=label_x, y=label_y, label=gene), size = 2, hjust=0) +

        # Draw lines from labels to points
        ggplot2::geom_segment(data=sig_down, aes(x = log2FoldChange, y = -log10(padj), xend = label_x, yend = label_y), color = "#004B7A", alpha = 0.5) +
        ggplot2::geom_segment(data=sig_up, aes(x = log2FoldChange, y = -log10(padj), xend = label_x, yend = label_y), color = "#A83708", alpha = 0.5) +


        # Add titles
        ggplot2::xlab("log2 fold change") +
        ggplot2::ylab("-log10 adjusted p-value") +
        ggplot2::ggtitle(paste0(title)) +
        ggprism::theme_prism(border = T,
                             base_fontface = "plain",
                             base_size = 6) +
        ggplot2::coord_cartesian(clip = "off") +
        ggplot2::expand_limits(x = c(min(df_2$log2FoldChange) - 6, max(df_2$log2FoldChange) + 6))

    return(output)
}

#' Volcano plot with labels of specific genes
#'
#' @param df a deseq2 results dataframe
#' @param padj_thres a value e.g. 0.05
#' @param logfc_thres a value e.g. 1.5
#' @param title title of plot e.g. "HFD_3 vs LFD"
#' @param genes_up vector of genes with positive log2foldchange
#' @param genes_down vector of genes with negative log2foldchange
#'
#' @return a ggplot2 plot
volcano_genes <- function(df, padj_thres, logfc_thres, title, genes_up, genes_down) {

    df_2 <- df %>%
        as.data.frame() %>%
        tibble::rownames_to_column("gene")

    top <- base::max(stats::na.omit(df_2$padj))
    bottom <- base::min(stats::na.omit(df_2$padj))

    # Subset for significant genes
    sig_down <- df_2 %>%
        dplyr::filter(gene %in% genes_down) %>%
        dplyr::arrange(padj, log2FoldChange)

    sig_up <- df_2 %>%
        dplyr::filter(gene %in% genes_up) %>%
        dplyr::arrange(padj, log2FoldChange)

    # Create label positions manually
    sig_down <- sig_down %>%
        mutate(label_x = min(na.omit(df_2$log2FoldChange)) +2,
               label_y = seq(-log10(bottom), -log10(top), length.out = nrow(sig_down)))

    sig_up <- sig_up %>%
        mutate(label_x = max(na.omit(df_2$log2FoldChange)) + (-2),
               label_y = seq(-log10(bottom), -log10(top), length.out = nrow(sig_up)))

    output <- ggplot2::ggplot(data= df_2, aes(x=log2FoldChange, y=-log10(padj))) +
        ggrastr::geom_point_rast(shape=21,color= "lightgrey", fill = "lightgray", size=0.1) +

        # Add color to significant results
        ggrastr::geom_point_rast(data=subset(df_2, padj <= padj_thres & log2FoldChange <= -logfc_thres), shape=21, color ='#004B7A', fill = '#004B7A', size = 0.5) +
        ggrastr::geom_point_rast(data=subset(df_2, padj <= padj_thres & log2FoldChange >= logfc_thres), shape=21, color ='#A83708', fill = '#A83708', size = 0.5) +

        # Add lines
        ggplot2::geom_vline(xintercept = c(0, logfc_thres, -logfc_thres),
                            linetype = c(1, 2, 2),
                            color = c("#868686FF", "#868686FF", "#868686FF")) +
        ggplot2::geom_hline(yintercept = -log10(padj_thres),
                            linetype = 2,
                            color = "#868686FF") +

        # Add number of significant results to plot
        ggplot2::geom_text(data= df_2 %>%
                               dplyr::filter(padj <= padj_thres & log2FoldChange <= -logfc_thres) %>%
                               dplyr::summarise(n = n()),
                           aes(label = paste0("n(genes) = ", n), x = -Inf, y = Inf), color = '#1465AC',
                           hjust = -0.5, vjust = 1.5, size = 2) +
        ggplot2::geom_text(data= df_2 %>%
                               dplyr::filter(padj <= padj_thres & log2FoldChange >= logfc_thres) %>%
                               dplyr::summarise(n = n()),
                           aes(label = paste0("n(genes) = ", n), x = -Inf, y = Inf), color = '#B31B21',
                           hjust = -0.5, vjust = 3, size = 2) +

        # Add text to top 10 significant genes
        ggplot2::geom_text(data=sig_down, aes(x=label_x, y=label_y, label=gene), size = 2, hjust=1) +
        ggplot2::geom_text(data=sig_up, aes(x=label_x, y=label_y, label=gene), size = 2, hjust=0) +

        # Draw lines from labels to points
        ggplot2::geom_segment(data=sig_down, aes(x = log2FoldChange, y = -log10(padj), xend = label_x, yend = label_y), color = "#004B7A", alpha = 0.5) +
        ggplot2::geom_segment(data=sig_up, aes(x = log2FoldChange, y = -log10(padj), xend = label_x, yend = label_y), color = "#A83708", alpha = 0.5) +


        # Add titles
        ggplot2::xlab("log2 fold change") +
        ggplot2::ylab("-log10 adjusted p-value") +
        ggplot2::ggtitle(paste0(title)) +
        ggprism::theme_prism(border = T,
                             base_fontface = "plain",
                             base_size = 6) +
        ggplot2::coord_cartesian(clip = "off") +
        ggplot2::expand_limits(x = c(min(df_2$log2FoldChange) - 6, max(df_2$log2FoldChange) + 6))

    return(output)
}

volcano_archr <- function(df, FDR_thres, logfc_thres, title) {


    output <-     ggplot2::ggplot(data= df, aes(x=Log2FC, y=-log10(FDR))) +
        ggplot2::geom_point() +
        ggplot2::geom_point(shape=21, color="lightgrey", fill= "lightgrey", size=2) +

        # add color to significant results
        ggplot2::geom_point(data=subset(df, FDR <= FDR_thres & Log2FC <= -logfc_thres), shape=21, color = "black", fill ='#004B7A', size = 2) +
        ggplot2::geom_point(data=subset(df, FDR <= FDR_thres & Log2FC >= logfc_thres), shape=21, color = "black", fill ='#A83708', size = 2) +

        # Add number of significant results to plot
        ggplot2::geom_text(data= df %>%
                               dplyr::filter(FDR <= FDR_thres & Log2FC <= -logfc_thres) %>%
                               dplyr::summarise(n = n()),
                           aes(label = n, x = -Inf, y = Inf), color = '#1465AC',
                           hjust = -0.5, vjust = 1.5) +
        ggplot2::geom_text(data= df %>%
                               dplyr::filter(FDR <= FDR_thres & Log2FC >= logfc_thres) %>%
                               dplyr::summarise(n = n()),
                           aes(label = n, x = -Inf, y = Inf), color = '#B31B21',
                           hjust = -0.5, vjust = 3) +
        # add text to top 10 significant SYMBOLs
        ggrepel::geom_text_repel(data= df %>%
                                     dplyr::filter(FDR <= FDR_thres &
                                                       Log2FC <= -logfc_thres) %>%
                                     dplyr::arrange(FDR, Log2FC) %>%
                                     head(10),
                                 aes(label = SYMBOL), size=3,
                                 max.overlaps = 100) +

        ggrepel::geom_text_repel(data= df %>%
                                     dplyr::filter(FDR <= FDR_thres &
                                                       Log2FC >= logfc_thres) %>%
                                     dplyr::arrange(FDR, Log2FC) %>%
                                     head(10),
                                 aes(label = SYMBOL), size=3,
                                 max.overlaps = 100) +

        # Add lines
        ggplot2::geom_vline(xintercept = c(0, logfc_thres, -logfc_thres),
                            linetype = c(1, 2, 2),
                            color = c("#868686FF", "#868686FF", "#868686FF")) +
        ggplot2::geom_hline(yintercept = -log10(FDR_thres),
                            linetype = 2,
                            color = "#868686FF") +

        # add titles
        ggplot2::xlab("log2 fold change") +
        ggplot2::ylab("-log10 adjusted p-value") +
        ggplot2::ggtitle(paste0(title)) +
        ggprism::theme_prism(border = T,
                             base_fontface = "plain",
                             base_size = 12) +
        ggplot2::coord_cartesian(clip = "off")

    return(output)
}


ma <- function(df, padj_thres, fc_thres, title){
    output <- df %>%
        ggpubr::ggmaplot(
            main = paste0(title),
            fdr = padj_thres,
            fc = fc_thres,
            size = 0.4,
            palette = c("#B31B21", "#1465AC", "darkgray"),
            genenames = as.vector(df$gene),
            legend = "top",
            top = 10,
            font.label = c("bold", 11),
            font.legend = "bold",
            font.main = "bold") +
        ggprism::theme_prism(border = T,
                             base_fontface = "plain",
                             base_size = 12) +
        ggplot2::coord_cartesian(clip = "off")
    return(output)
}


# determine m1 and number of clusters
mfuzzy_params <- function(genes, df, cluster) {

    print(paste0("Mean norm counts of ", cluster))

    df_2 <- df %>%
        as.data.frame() %>%
        tibble::rownames_to_column("gene") %>%
        dplyr::filter(gene %in% genes) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(LFD = mean(c(LFD_R1, LFD_R2, LFD_R3)),
                      HFD_1 = mean(c(HFD_1_R1, HFD_1_R2)),
                      HFD_3 = mean(c(HFD_3_R1, HFD_3_R2)))  %>%
        tibble::column_to_rownames("gene") %>%
        dplyr::select(LFD, HFD_1, HFD_3) %>%
        as.matrix()

    # Create expression set object
    exprSet <- Biobase::ExpressionSet(assayData=df_2)

    # filter genes with low variance (0 genes)
    tmp=Mfuzz::filter.std(exprSet, min.std=0, visu=FALSE)

    # z-scores
    exprSet.s=Mfuzz::standardise(tmp)

    # Estimate fuzzifier
    # For fuzzyfier m, we would like to choose a value which prevents clustering of random
    # data. Note, that fuzzy clustering can be tuned in such manner, that random data is not
    # clustered. This is a clear advantage to hard clustering (such as k-means), which commonly
    # detects clusters even in random data

    # this function gives the minimum fuzzifier value which prevents clustering of randomized data (3.619766)
    m1=Mfuzz::mestimate(exprSet.s)
    print(paste0("Fuzzyfier is ", m1))

    # evaluate optimal number of clusters by calculating minimum distance to centroids for 1:20 clusters.
    # when the optimal number of clusters is found, the minimunm distance to centroids will dcrease more slowly
      
    temp <- Mfuzz::Dmin(exprSet.s, m=m1, crange=base::seq(2,20,1), repeats = 15, visu = FALSE)

        # calculate differences in dmin
    diff <- temp[-length(temp)] - temp[-1]
    diff <- c(diff, NA)
    
    dmin_df <- base::data.frame(
      cluster_1 = base::seq(2,20,1),
      cluster_2 = (1 + base::seq(2,20,1)),
      dmin_1 = temp,
      dmin_2 = c(temp[2:20]),
      diff = diff
    ) 

    # create plot
    plot <- dmin_df %>% ggplot2::ggplot() +
        geom_point(aes(x = cluster_1, y = dmin_1), color = "red") +
        geom_point(aes(x = cluster_1, y = dmin_2), color = "blue") +
        geom_text(aes(label = cluster_1, x = cluster_1, y = (0.1 + dmin_1))) +
        geom_text(aes(label = cluster_2, x = cluster_1, y = (dmin_2 -0.1))) +
        labs( y = "minimum centroid distance",
              x = "Cluster",
              title = paste0(cluster),
              subtitle = paste0("fuzzyfier = ", m1)) +
        geom_segment(aes(x = cluster_1,
                         xend = cluster_1,
                         y = dmin_2,
                         yend = dmin_1),
                     col = 'black',
                     linetype = "dashed") +
      my_theme()
    output <- list("stand_expr" = exprSet.s,
                   "m1" = m1,
                   "min_c_dis" = dmin_df,
                   "plot" = plot)

    return(output)

}



# plot dmins
dmin_plot <- function(plot){
    output <- plot[["plot"]]
    return(output)
}

# cluster
mfuzz_cluster <- function(exp, n){
    output <- mfuzz(exp[["stand_expr"]],c=n,m=exp[["m1"]])
    return(output)
}

# plot cluster
plot_cluster <- function(norm, genes, c, title){

    # mean normalized expression z-scores
    df <- norm %>%
        as.data.frame() %>%
        tibble::rownames_to_column("gene") %>%
        dplyr::filter(gene %in% genes) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(LFD = mean(c(LFD_R1, LFD_R2, LFD_R3)),
                      HFD_1 = mean(c(HFD_1_R1, HFD_1_R2)),
                      HFD_3 = mean(c(HFD_3_R1, HFD_3_R2)))  %>%
        tibble::column_to_rownames("gene") %>%
        dplyr::select(LFD, HFD_1, HFD_3) %>%
        t() %>%
        scale() %>%
        t() %>%
        as.data.frame() %>%
        tibble::rownames_to_column("gene") %>%
        tidyr::pivot_longer(-gene, names_to = "condition", values_to = "exp")

    # Extract cluster information
    cl_2 <- c$cluster  %>%
        as.data.frame() %>%
        tibble::rownames_to_column("gene")
    colnames(cl_2) <- c("gene", "membership")

    # number of genes in each cluster
    label <- cl_2 %>%
        dplyr::group_by(membership) %>%
        dplyr::summarise(n = n())

    # concanate gene and membership columns
    cl_2 <- cl_2 %>%
        dplyr::mutate(mem = paste0(gene, "_", membership)) %>%
        dplyr::select(mem)

    # combine cluster membership assignment and membership score for that cluster as well as scaled gene expression across condition
    output <- c$membership  %>%
        as.data.frame() %>%
        tibble::rownames_to_column("gene") %>%
        tidyr::pivot_longer(-gene, names_to = "membership", values_to = "score") %>%
        dplyr::mutate(mem = paste0(gene, "_", membership)) %>%
        dplyr::right_join(y = cl_2, by = "mem") %>%
        dplyr::select(-mem) %>%
        dplyr::left_join(y = df, by = "gene") %>%
        dplyr::arrange(membership, desc(score)) %>%
        dplyr::mutate(condition = factor(condition, levels = condition_levels)) %>%
        ggplot2::ggplot(aes(y = exp, x = condition, color = score)) +
        ggplot2::geom_line(aes(group=gene)) +
        geom_hline(yintercept=0,
                   color = "black", linewidth=0.5) +
        ggplot2::labs(
            y = "Z-scores (Mean pseudobulk expression)",
            x = "Samples",
            title = paste0(title),
            color = "Membership\n Score") +
        ggplot2::stat_summary(fun.data = mean_se, colour="black",
                              aes(group=1),
                              geom="line",
                              lwd=1,
                              lty=2) +
        ggplot2::facet_wrap(~membership, ncol = 1) +
        ggplot2::scale_colour_gradientn(colours = rev(RColorBrewer::brewer.pal(n = 11, name = "Spectral"))) +
        ggprism::theme_prism(border = T,
                             base_fontface = "plain",
                             base_size = 12) +
        ggplot2::coord_cartesian(clip = "off") +
        ggplot2::theme(axis.title.x = element_blank(),
                       legend.title = element_text()) +
        ggplot2::geom_text(
            data    = label,
            mapping = aes(x = -Inf, y = -Inf, label = n),
            hjust   = -0.1,
            vjust   = -1,
            col = "black")
    return(output)

}

# plot cluster
plot_cluster_colors <- function(norm, genes, c, title, color){

    # mean normalized expression z-scores
    df <- norm %>%
        as.data.frame() %>%
        tibble::rownames_to_column("gene") %>%
        dplyr::filter(gene %in% genes) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(LFD = mean(c(LFD_R1, LFD_R2, LFD_R3)),
                      HFD_1 = mean(c(HFD_1_R1, HFD_1_R2)),
                      HFD_3 = mean(c(HFD_3_R1, HFD_3_R2)))  %>%
        tibble::column_to_rownames("gene") %>%
        dplyr::select(LFD, HFD_1, HFD_3) %>%
        t() %>%
        scale() %>%
        t() %>%
        as.data.frame() %>%
        tibble::rownames_to_column("gene") %>%
        tidyr::pivot_longer(-gene, names_to = "condition", values_to = "exp")

    # Extract cluster information
    cl_2 <- c$cluster  %>%
        as.data.frame() %>%
        tibble::rownames_to_column("gene")
    colnames(cl_2) <- c("gene", "membership")

    # number of genes in each cluster
    label <- cl_2 %>%
        dplyr::group_by(membership) %>%
        dplyr::summarise(n = n())

    # concanate gene and membership columns
    cl_2 <- cl_2 %>%
        dplyr::mutate(mem = paste0(gene, "_", membership)) %>%
        dplyr::select(mem)

    # combine cluster membership assignment and membership score for that cluster as well as scaled gene expression across condition
    output <- c$membership  %>%
        as.data.frame() %>%
        tibble::rownames_to_column("gene") %>%
        tidyr::pivot_longer(-gene, names_to = "membership", values_to = "score") %>%
        dplyr::mutate(mem = paste0(gene, "_", membership)) %>%
        dplyr::right_join(y = cl_2, by = "mem") %>%
        dplyr::select(-mem) %>%
        dplyr::left_join(y = df, by = "gene") %>%
        dplyr::arrange(membership, desc(score)) %>%
        dplyr::mutate(condition = factor(condition, levels = condition_levels)) %>%
        ggplot2::ggplot(aes(y = exp, x = condition, color = membership)) +
        ggplot2::geom_line(aes(group=gene)) +
        geom_hline(yintercept=0,
                   color = "black", linewidth=0.5) +
        ggplot2::labs(
            y = "Z-scores (Mean pseudobulk expression)",
            x = "Samples",
            title = paste0(title),
            color = "Membership\n Score") +
        ggplot2::stat_summary(fun.data = mean_se, colour="black",
                              aes(group=1),
                              geom="line",
                              lwd=1,
                              lty=2) +
        ggplot2::facet_wrap(~membership, ncol = 1) +
        ggplot2::scale_colour_manual(values = color) +
        ggprism::theme_prism(border = T,
                             base_fontface = "plain",
                             base_size = 12) +
        ggplot2::coord_cartesian(clip = "off") +
        ggplot2::theme(axis.title.x = element_blank(),
                       legend.title = element_text()) +
        ggplot2::geom_text(
            data    = label,
            mapping = aes(x = -Inf, y = -Inf, label = n),
            hjust   = -0.1,
            vjust   = -1,
            col = "black")
    return(output)

}
## Enrichment analysis

genes_enrichr <- function(c){
    cl_2 <- c$cluster  %>%
        as.data.frame() %>%
        tibble::rownames_to_column("gene")

    colnames(cl_2) <- c("gene", "membership")

    output <- cl_2 %>%
        dplyr::group_by(membership) %>%
        dplyr::group_split() %>%
        purrr::modify_depth(1, ~ dplyr::select(., gene) %>%
                                purrr::as_vector() %>%
                                unname() %>%
                                unique())
    names(output) <- c(1:length(output))
    return(output)
}


plot_enrich <- function(df, n_term, title){
    # calcualte new p-value
    df_2 <- df %>%
        dplyr::mutate(gene_count = as.numeric(sub("/\\d+$", "", as.character(Overlap))),
                      gene_path_count =  as.numeric(sub("^\\d+/", "", as.character(Overlap))),
                      ratio = gene_count/gene_path_count,
                      gg_term = stringr::str_wrap(Term, width = 11),
                      gg_term = factor(gg_term, levels = rev(gg_term))) %>%
        dplyr::filter(gene_path_count > 2 & gene_path_count < 100)

    df_2$new_pval <- p.adjust(df_2$P.value, method = "bonferroni")

    output <- df_2 %>%
        dplyr::arrange(new_pval) %>%
        head(n = n_term) %>%
        ggplot2::ggplot(aes(x = Combined.Score, y = gg_term, fill = new_pval)) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::geom_text(aes(label = Overlap), hjust = 0) +
        ggplot2::scale_fill_gradient(low = "#A83708", high = "#004B7A") +
        ggplot2::guides(fill = guide_colorbar(title = "Adjusted P-value", reverse = TRUE)) +
        ggplot2::labs(title = paste0(title),
                      x = "Combined score",
                      fill = "Adjusted p-value") +
        ggprism::theme_prism(border = T,
                             base_fontface = "plain",
                             base_size = 12) +
        ggplot2::theme(axis.title.y = element_blank(),
                       legend.title = element_text())

    return(output)
}

genes_enrich <- function(df, n_term){
    df_2 <- df %>%
        dplyr::mutate(gene_count = as.numeric(sub("/\\d+$", "", as.character(Overlap))),
                      gene_path_count =  as.numeric(sub("^\\d+/", "", as.character(Overlap))),
                      ratio = gene_count/gene_path_count,
                      gg_term = stringr::str_wrap(Term, width = 11),
                      gg_term = factor(gg_term, levels = rev(gg_term))) %>%
        dplyr::filter(gene_path_count >2 & gene_path_count < 100)

    df_2$new_pval <- p.adjust(df_2$P.value, method = "bonferroni")

    output <- df_2 %>%
        dplyr::mutate(sig = case_when(new_pval > 0.05 ~ "NS",
                                      new_pval <=0.05 ~ "S"),
                      sig = factor(sig, levels = c("S", "NS"))) %>%
        dplyr::arrange(new_pval) %>%
        head(n = n_term)

    return(output)
}

heatmap_cluster <- function(df){
    output <- df$cluster  %>%
        as.data.frame() %>%
        tibble::rownames_to_column("gene")

    colnames(output) <- c("gene", "membership")

    return(output)

}

heat_row_anno <- function(df, norm_counts){
    df_norm <- norm_counts %>%
        as.data.frame() %>%
        tibble::rownames_to_column("gene")

    output <- df %>%
        as.data.frame() %>%
        dplyr::select(membership, gene) %>%
        distinct() %>%
        dplyr::left_join(y = df_norm, by = "gene") %>%
        tibble::column_to_rownames("gene") %>%
        dplyr::select(membership)
    return(output)
}

gg_fuzzy_cluster_heatmap <- function(df, heat_color, breaks, title, gene, row_anno, col_anno, anno_color){
    output <- df  %>%
        as.data.frame() %>%
        tibble::rownames_to_column("genes") %>%
        dplyr::filter(genes %in% gene) %>%
        tibble::column_to_rownames("genes")  %>%
        t() %>%
        scale() %>%
        t() %>%
        ComplexHeatmap::pheatmap(
            color = heat_color,
            main = paste0(title),
            breaks = breaks,
            show_rownames = F,
            cluster_cols = FALSE,
            cluster_rows = FALSE,
            annotation_color = anno_color,
            annotation_col = col_anno,
            column_split = col_anno$condition,
            row_split = row_anno$membership,
            annotation_row = row_anno,
            heatmap_legend_param = list(title = "z-score (Norm counts)"))
    return(output)
}

heat_of_genes <- function(df, heat_color, breaks, title, gene, col_anno, anno_color, row_anno){
    df_2 <- df  %>%
        as.data.frame() %>%
        tibble::rownames_to_column("genes") %>%
        dplyr::filter(genes %in% gene) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(LFD = mean(c(LFD_R1, LFD_R2, LFD_R3)),
                      HFD_1 = mean(c(HFD_1_R1, HFD_1_R2)),
                      HFD_3 = mean(c(HFD_3_R1, HFD_3_R2)))  %>%
        tibble::column_to_rownames("genes") %>%
        dplyr::select(LFD, HFD_1, HFD_3) %>%
        t() %>%
        scale() %>%
        t()

    df_3 <-
        df_2 %>%
        as.data.frame() %>%
        tibble::rownames_to_column("gene")

    row_anno_2 <- row_anno[match(df_3$gene, row_anno$gene),]

    rownames(row_anno_2) <- NULL

    row_anno_3 <- row_anno_2 %>%
        tibble::column_to_rownames("gene")



        output <- df_2 %>%
        ComplexHeatmap::pheatmap(
            color = heat_color,
            main = paste0(title),
            breaks = breaks,
            show_rownames = T,
            cluster_cols = FALSE,
            cluster_rows = T,
            row_title_gp = gpar(fontsize = 8),
            row_title_rot = 0,
            annotation_colors = anno_color,
            annotation_row = row_anno_3,
            row_split = row_anno_3$class,
            heatmap_legend_param = list(title = "z-score (Norm counts)"))
    return(output)
}


#' Get genes present in each gene cluster
#'
#' @param df a dataframe with a column of gene names called "gene" and a column of clusters called "membership"
#'
#' @return a list of gene names in a vector
cluster_genes <- function(df){
    output <- df %>%
        dplyr::group_by(membership) %>%
        dplyr::group_split() %>%
        purrr::modify_depth(1, ~ dplyr::select(., gene) %>%
                                purrr::as_vector() %>%
                                unname() %>%
                                unique())
}


num_genes <- function(df, cluster){
    output <- df %>%
        dplyr::summarize(n = n()) %>%
        dplyr::mutate(cluster = paste0(cluster))
    return(output)
}

gg_num_genes <- function(df, cell_type_levels){
  output <- df %>% 
    dplyr::mutate(cluster = factor(cluster, levels = cell_type_levels)) %>%
    ggplot2::ggplot(aes(x = cluster, y = n, fill = cluster)) +
    ggplot2::geom_bar(position = position_dodge(), stat = "identity") +
    ggplot2::scale_fill_manual(values = cluster_anno) +
    ggplot2::labs(y = "# DEGs") +
    my_theme() +
    ggplot2::theme(legend.position = "none",
                   axis.title.x = element_blank())
  
  return(output)
  
}

gg_augur_score <- function(df, cell_type_levels){
    output <- df %>% 
        dplyr::mutate(cell_type = factor(cell_type, levels = cell_type_levels)) %>%
      ggplot2::ggplot(aes(x = cell_type, y = auc, fill = cell_type)) +
        ggplot2::geom_bar(position = position_dodge(), stat = "identity") +
        ggplot2::scale_fill_manual(values = cluster_anno) +
        ggplot2::labs(y = "Prioritization score \n(snRNA-seq)") +
      my_theme() +
        ggplot2::theme(legend.position = "none",
                       axis.title.x = element_blank())

    return(output)

}

#' Perfrom GO-term analysis with clusterprofiler for a wald test
#'
#' @param test_genes - vector of genes to be tested
#' @param test_cell - name of gene set to be tested e.g. "gene_set_1
#' @param bg_genes - background genes (all genes that were tested) - dataframe where genes are rownames
#' @param bg_cell - name of background genes e.g. "beta cell"
#'
#' @return a cluster profiler result
go_term_wald <- function(test_genes, test_cell, bg_genes, bg_cell){

    print(paste0("testing:", test_cell))
    # Get genes for each cell type

    # add entrez IDs to genes

    genes <- BiocGenerics::Map(function(df){
        if (length(rownames(df) > 0)){
            df$entrez = AnnotationDbi::mapIds(org.Mm.eg.db::org.Mm.eg.db,
                                              keys=rownames(df),
                                              column="ENTREZID",
                                              keytype="SYMBOL",
                                              multiVals="first")

            output <- df
        } else {
            output <- "NA"
        }

        return(output)
    },
    df = test_genes)

    # Add entrez IDs to background genes

    bg_genes <- BiocGenerics::Map(function(df){
        df$entrez = AnnotationDbi::mapIds(org.Mm.eg.db::org.Mm.eg.db,
                                          keys=rownames(df),
                                          column="ENTREZID",
                                          keytype="SYMBOL",
                                          multiVals="first")

        return(df)
    },
    df = bg_genes)

    print(paste0("background:", bg_cell))

    # Create universe - this should be any gene that could have been positive
    # And I guess that should be all genes that were tested

    my_universe <- BiocGenerics::Map(function(df){
        output <- df %>%
            as.data.frame() %>%
            filter(!is.na(entrez)) %>%
            pull(entrez) %>%
            unlist() %>%
            unname() %>%
            unique()

        return(output)
    },
    df = bg_genes)



    gene_list <- genes
    # Pathway analysis with clusterprofiler + jaccard index

    # Jaccard
    # k is the overlap between your genes-of-interest and the geneset
    # n is the number of all unique genes-of-interest

    # BgRatio=M/N

    # M is the number of genes within each geneset
    # N is the number of all unique genes across all genesets (universe)
    output <- BiocGenerics::Map(function(df, universe){
        if (is.character(df) == TRUE) {
            output <- "no genes"
        } else {
        output <- df %>%
            dplyr::select(entrez) %>%
            purrr::as_vector() %>%
            unname() %>%
            unique() %>%
            clusterProfiler::enrichGO(OrgDb = "org.Mm.eg.db",
                                      pAdjustMethod = "fdr",
                                      ont = "BP",
                                      pvalueCutoff  = 0.2,
                                      qvalueCutoff  = 0.2,
                                      readable      = TRUE,
                                      universe = universe
            ) %>%
            dplyr::mutate(k = as.numeric(sub("/\\d+$", "", as.character(GeneRatio))),
                          n = as.numeric(sub("^\\d+/", "", as.character(GeneRatio))),
                          M = as.numeric(sub("/\\d+$", "", as.character(BgRatio))),
                          N = as.numeric(sub("^\\d+/", "", as.character(BgRatio))),
                          j_path = k/((n + M) - k),
                          j_bg = M/((N + M) - M),
                          is_sig = case_when(p.adjust <= 0.0001 ~ "****",
                                             p.adjust <= 0.001 ~ "****",
                                             p.adjust <= 0.01 ~ "**",
                                             p.adjust <= 0.05 ~ "*",
                                             p.adjust > 0.05 ~ "ns"))
        }

        return(output)

    },
    df = gene_list,
    universe = my_universe)


    return(output)
}

gg_go_term_wald <- function(res, way) {
    plot <- BiocGenerics::Map(function(res_2, cluster) {
        if (is.character(res_2) == TRUE) {
            x <- rnorm(20)
            y <- rnorm(20, 1, 0.5)
            df <- data.frame(x, y)

            plot <-
                ggplot2::ggplot(df, aes(x, y)) +
                ggplot2::geom_blank() +
                ggprism::theme_prism(
                    border = T,
                    base_fontface = "plain",
                    base_size = 12
                ) + ggtitle("no data")
        } else {
            # Top 5 most significant pathways
            df <- res_2@result %>%
                dplyr::arrange(p.adjust) %>%
                head(n = 5)

            # Plot
            if (way == "up") {
                output <- df %>%
                    ggplot2::ggplot(aes(
                        y = fct_reorder(Description, j_path),
                        x = j_path
                    )) +
                    ggplot2::geom_bar(
                        stat = "identity",
                        position = position_dodge(),
                        fill = "#B31B21"
                    ) +
                    ggplot2::geom_text(
                        aes(label = is_sig),
                        data = df,
                        color = "black",
                        size = 5
                    ) +
                    ggplot2::scale_y_discrete(position = "right") +
                    xlab("Jaccard similarity index") +
                    ylab(NULL) +
                    ggtitle(paste0(cluster, " Up")) +
                    ggprism::theme_prism(
                        border = T,
                        base_fontface = "plain",
                        base_size = 12
                    ) +
                    ggplot2::theme(legend.title = element_text())
            } else {
                output <- df %>%
                    ggplot2::ggplot(aes(
                        y = fct_reorder(Description, j_path),
                        x = j_path
                    )) +
                    ggplot2::geom_bar(
                        stat = "identity",
                        position = position_dodge(),
                        fill = "#1465AC"
                    ) +
                    ggplot2::geom_text(
                        aes(label = is_sig),
                        data = df,
                        color = "black",
                        size = 5
                    ) +
                    ggplot2::scale_x_reverse() +
                    xlab("Jaccard similarity index") +
                    ylab(NULL) +
                    ggtitle(paste0(cluster, " Down")) +
                    ggprism::theme_prism(
                        border = T,
                        base_fontface = "plain",
                        base_size = 12
                    ) +
                    ggplot2::theme(legend.title = element_text())
            }

            return(output)
        }
    },
    res_2 = res,
    cluster = names(res))

    output <- plot %>% patchwork::wrap_plots(ncol = 1)

    return(output)
}

