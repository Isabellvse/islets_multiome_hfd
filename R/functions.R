#' Create Directories if They Do Not Exist
#'
#' This function checks if each directory in the provided list exists. If any directory
#' does not exist, it creates the directory (along with any necessary parent directories).
#' It prints a message for each directory indicating whether it was created or already exists.
#'
#' @param output_dir A character vector or list of directory paths to be checked and created.
#'
#' @return NULL This function performs a side-effect (creating directories) and does not return anything.
#' 
#' @export
#'
#' @examples
#' # Define a list of directories
#' output_dirs <- list("data/quality_control/rna", "data/plots", "data/other_output")
#' 
#' # Call the function to create the directories
#' create_directories(output_dirs)
create_directories <- function(output_dir) {
  purrr::walk(output_dir, ~{
    if (!dir.exists(.x)) {
      dir.create(.x, recursive = TRUE)  # create the directory if it doesn't exist
      print(paste0(.x, " has been created!"))
    } else {
      print(paste0(.x, " already exists!"))
    }
  })
}




#' Process Doublets in Seurat Objects
#'
#' This function processes Seurat objects to identify and adjust for doublets using the DoubletFinder package.
#'
#' @param s_obj A Seurat object to be processed.
#' @param sample_name A character string representing the name of the sample.
#' @param df_freq A data frame containing the sample names and their corresponding doublet frequencies.
#'
#' @return A Seurat object with doublets identified and adjusted.
#' @export
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' seurat_list_doublets <- pmap(list(seurat_list, names(seurat_list), list(df_freq)), process_doublets)
#' }
process_doublets <- function(s_obj, sample_name, df_freq) {
  # pK Identification (no ground truth)
  base::print(paste0("pK Identification (no ground truth)"))
  sweep.res.list <- DoubletFinder::paramSweep_v3(
    s_obj,
    PCs = 1:20,
    num.cores = parallel::detectCores() - 1
  )
  
  sweep.stats <- DoubletFinder::summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- DoubletFinder::find.pK(sweep.stats)
  
  # Optimal pK
  pK <- bcmvn %>%
    dplyr::filter(BCmetric == max(BCmetric)) %>%
    dplyr::select(pK)
  
  opti_pK <- base::as.numeric(as.character(pK[[1]]))
  base::print(base::paste0("optimal pK = ", opti_pK))
  
  # Homotypic doublet proportion estimate
  annotations <- s_obj@meta.data$seurat_clusters
  homotypic.prop <- DoubletFinder::modelHomotypic(annotations)
  
  # Get doublet frequency
  freq <- dplyr::filter(df_freq, sample == sample_name)$freq
  base::print(base::paste0("frequency used:", freq))
  
  # Adjust number of doublets, based on homotypic doublets
  nExp_poi <- base::round(freq * base::nrow(s_obj@meta.data))
  nExp_poi.adj <- base::round(nExp_poi * (1 - homotypic.prop))
  
  # Plot optimal pK and number of expected doublets
  plot <- ggplot2::ggplot(bcmvn, aes(pK, BCmetric, group = 1)) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::labs(subtitle = base::paste0("pK: ", opti_pK, " Number of expected doublets: ", nExp_poi.adj)) +
    my_theme()
  
  # Save plot
  ggplot2::ggsave(plot, filename = base::paste0(here::here("data/quality_control/rna/"), sample_name, "_DoubletFinder_optimal_pK.pdf"), width = 6, height = 6)
  
  # Run DoubletFinder
  output <- DoubletFinder::doubletFinder_v3(
    s_obj,
    PCs = 1:20,
    pN = 0.25,
    pK = opti_pK,
    nExp = nExp_poi.adj,
    reuse.pANN = FALSE,
    sct = FALSE
  )
  return(output)
}

my_theme <- function() {
  ggplot2::theme_classic(base_size = 6) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line(color = "grey85"),
      panel.grid.minor = ggplot2::element_line(color = "grey95")
    )
}

my_theme_void <- function(remove_strip_text = TRUE) {
  base_theme <- ggplot2::theme_void(base_size = 6)
  
  if (remove_strip_text) {
    base_theme <- base_theme + ggplot2::theme(strip.text = ggplot2::element_blank())
  }
  
  base_theme
}


#' Function to plot several of archrs ggPoint plots.
#'
#' @param a dataframe / tibble containing the colunms log10.nFrags. .df$TSSEnrichment
#' @param .title The title of your plot "some title"
#'
#' @return a ggplot
ggPoint_plot <- function(.df, .title){
  a <- BiocGenerics::as.vector(.df$log10.nFrags.)
  b <- BiocGenerics::as.vector(.df$TSSEnrichment)
  output <- ArchR::ggPoint(
    x = a,
    y = b,
    title = .title,
    size = 0.1,
    colorDensity = TRUE,
    continuousSet = "sambaNight",
    xlabel = "Log10 Unique Fragments",
    ylabel = "TSS Enrichment",
    xlim = c(log10(500), quantile(a, probs = 0.99)),
    ylim = c(0, quantile(b, probs = 0.99))
  ) + geom_hline(yintercept = 15, lty = "dashed") + geom_vline(xintercept = 3.3, lty = "dashed") +
    my_theme()
  return(output)
}

#' Convert ArchR barcode syntax (LFD_R1#BARCODE-1) to seurat (RNA syntax (LFD_R1_BARCODE)
#'
#' @param vec a character vector containing barcodes in Archr syntax
#'
#' @return a character vector with barcodes in RNA syntax (Seurat)

atac_to_rna_syntax <- function(vec){
  output <-
    vec %>% stringi::stri_replace_last_regex("#", "_") %>%
    stringr::str_remove("-1")
  return(output)
}

## function adapted from: http://www.github.com/madsen-lab/JOINTLY_reproducibility
ilisi_comp <- function(x, metadata, dims, batch_var) {
  # Ensure dims is correctly defined as a vector
  dims <- dims
  
  # Compute LISI
  space <- x
  LISIS <- lisi::compute_lisi(space[,1:dims], metadata, c(batch_var))
  
  # Calculate global iLISI
  global_iLISI <- (median(LISIS[,1]) - 1) / (length(unique(metadata[,batch_var])) - 1)
  
  # Return global iLISI
  list(global_iLISI = global_iLISI)
}

## function adapted from: http://www.github.com/madsen-lab/JOINTLY_reproducibility
evaluateEmbedding <- function(x, metadata, dims, batch_var, label_var, cl.min = 1, cl.max = 50) {
  # Find which number of dimensions to use
  dims <- dims
  
  # Setup to capture clustering metrics
  global_aris <- c()
  global_nmis <- c()
  global_vs <- c()
  global_completeness <- c()
  global_homogeneity <- c()
  dataset_aris <- c()
  dataset_nmis <- c()
  dataset_vs <- c()
  dataset_completeness <- c()
  dataset_homogeneity <- c()
  cluster_list <- list()
  
  # Evaluate clustering metrics
  dend <- HGC::HGC.dendrogram(G = HGC::SNN.Construction(x[,1:dims]))
  for (cl in seq(cl.min, cl.max, 1)) {
    # Cluster
    if (cl == 1) {
      clusters <- rep(1, nrow(metadata))
      names(clusters) <- rownames(metadata)
    } else {
      clusters <- cutree(dend, k = cl)
      names(clusters) <- rownames(metadata)
    }
    cluster_list[[length(cluster_list)+1]] <- clusters
    
    # Capture global metrics
    global_aris <- c(global_aris, aricode::ARI(clusters, factor(metadata[,label_var])))
    global_nmis <- c(global_nmis, aricode::NMI(clusters, factor(metadata[,label_var])))
    global_vs <- c(global_vs, clevr::v_measure(clusters, factor(metadata[,label_var])))
    global_completeness <- c(global_completeness, clevr::completeness(clusters, factor(metadata[,label_var])))
    global_homogeneity <- c(global_homogeneity, clevr::homogeneity(clusters, factor(metadata[,label_var])))
    
    # Setup to capture dataset metrics
    ds_aris <- c()
    ds_nmis <- c()
    ds_vs <- c()
    ds_completeness <- c()
    ds_homogeneity <- c()
    
    # Calculate dataset metrics
    for (ds in unique(metadata[,batch_var])) {
      md <- metadata[ metadata[,batch_var] == ds,]
      clusters.ds <- clusters[ names(clusters) %in% rownames(md) ]
      clusters.ds <- clusters.ds[ match(rownames(md), names(clusters.ds))]
      ds_aris <- c(ds_aris, aricode::ARI(clusters.ds, factor(md[,label_var])))
      ds_nmis <- c(ds_nmis, aricode::NMI(clusters.ds, factor(md[,label_var])))
      ds_vs <- c(ds_vs, clevr::v_measure(clusters.ds, factor(md[,label_var])))
      ds_completeness <- c(ds_completeness, clevr::completeness(clusters.ds, factor(md[,label_var])))
      ds_homogeneity <- c(ds_homogeneity, clevr::homogeneity(clusters.ds, factor(md[,label_var])))
    }
    
    # Capture dataset metrics
    dataset_aris <- c(dataset_aris, min(ds_aris))
    dataset_nmis <- c(dataset_nmis, min(ds_nmis))
    dataset_vs <- c(dataset_vs, min(ds_vs))
    dataset_completeness <- c(dataset_completeness, min(ds_completeness))
    dataset_homogeneity <- c(dataset_homogeneity, min(ds_homogeneity))
  }
  
  ## cLISI and iLISI
  space <- x
  LISIS <- lisi::compute_lisi(space[,1:dims], metadata, c(batch_var, label_var))
  global_cLISI <- (length(unique(metadata[,label_var]))-median(LISIS[,2])) / (length(unique(metadata[,label_var])) - 1)
  global_iLISI <- (median(LISIS[,1])-1) / (length(unique(metadata[,batch_var])) - 1)
  dataset_cLISI <- c()
  dataset_iLISI <- c()
  for (ds in unique(metadata[,batch_var])) {
    dataset_cLISI <- c(dataset_cLISI, (length(unique(metadata[ metadata[,batch_var] == ds,label_var]))-median(LISIS[ rownames(LISIS) %in% rownames(metadata[ metadata[,batch_var] == ds,]),2])) / (length(unique(metadata[ metadata[,batch_var] == ds,label_var])) - 1))
    dataset_iLISI <- c(dataset_iLISI, (median(LISIS[ rownames(LISIS) %in% rownames(metadata[ metadata[,batch_var] == ds,]),1])-1) / (length(unique(metadata[,batch_var])) - 1))
  }
  
  ## Label and batch ASW
  # Subsample if too large
  if (dim(space)[1] >= 40000) {
    smps <- c()
    pr_ds_cells <- ceiling(table(metadata[,batch_var]) / (sum(table(metadata[,batch_var])) / 40000))
    for (ds in unique(metadata[,batch_var])) {
      pr_ct_cells <- ceiling(table(metadata[ metadata[,batch_var] == ds, label_var]) / (sum(table(metadata[ metadata[,batch_var] == ds, label_var]))  / pr_ds_cells[which(names(pr_ds_cells) == ds)]))
      ds.idx <- which(metadata[,batch_var] == ds)
      for (ct in unique(metadata[ metadata[,batch_var] == ds,label_var])) {
        ct.idx <- which(metadata[ metadata[,batch_var] == ds,label_var] == ct)
        selected <- sample(ct.idx, size = pr_ct_cells[names(pr_ct_cells) == ct], replace = FALSE)
        smps <- c(smps, ds.idx[selected])
      }
    }
    smps <- unique(smps)
    md.subset <- metadata[ rownames(metadata) %in% rownames(metadata)[smps],]
  } else {
    md.subset <- metadata
  }
  
  # Distance matrix
  space <- x
  space <- space[ rownames(space) %in% rownames(md.subset),]
  space <- space[ match(rownames(md.subset), rownames(space)),]
  dist.mat <- Rfast::Dist(space[,1:dims])
  
  # Batch ASW
  sil <- cluster::silhouette(as.numeric(factor(md.subset[,batch_var], labels = seq(length(unique(md.subset[,batch_var]))))), dmatrix = dist.mat, do.col.sort = FALSE)
  sil <- sil[,"sil_width"]
  sil <- abs(sil)
  avg.sil <- c()
  for (cl in unique(md.subset[,label_var])) {
    avg.sil <- c(avg.sil, sum(1 - sil[ which(md.subset[,label_var] == cl)]) / length( which(md.subset[,label_var] == cl)))
  }
  global.basw <- mean(avg.sil)
  dataset.basw <- c()
  for (ds in unique(md.subset[,batch_var])) {
    avg.sil <- c()
    for (cl in unique(md.subset[md.subset[,batch_var] == ds,label_var])) {
      avg.sil <- c(avg.sil, sum(1 - sil[ which(md.subset[,label_var] == cl & md.subset[,batch_var] == ds)]) / length( which(md.subset[,label_var] == cl & md.subset[,batch_var] == ds)))
    }
    dataset.basw <- c(dataset.basw, mean(avg.sil))
  }
  
  # Cell type ASW
  sil <- cluster::silhouette(as.numeric(factor(md.subset[,label_var], labels = seq(length(unique(md.subset[,label_var]))))), dmatrix = dist.mat, do.col.sort = FALSE)
  sil <- sil[,"sil_width"]
  global.casw <- mean((sil + 1) / 2)
  dataset.casw <- c()
  for (ds in unique(md.subset[,batch_var])) {
    dataset.casw <- c(dataset.casw, mean((sil[ which(md.subset[,batch_var] == ds)] + 1) / 2))
  }
  
  # Format results
  asw <- data.frame(celltype = dataset.casw, batch = dataset.basw)
  lisi <- data.frame(celltype = dataset_cLISI, batch = dataset_iLISI)
  clustering_metrics <- data.frame(globalARI = global_aris, globalNMI = global_nmis, globalVmeasure = global_vs, globalHomogeneity = global_homogeneity, globalCompleteness = global_completeness, datasetARI = dataset_aris, datasetNMI = dataset_nmis, datasetVmeasure = dataset_vs, datasetHomogeneity = dataset_homogeneity, datasetCompleteness = dataset_completeness)
  summary <- data.frame(globalARI = max(global_aris), datasetARI = max(dataset_aris), globalNMI = max(global_nmis), datasetNMI = max(dataset_nmis), globalVmeasure = max(global_vs), datasetVmeasure = max(dataset_vs), globalHomogeneity = max(global_homogeneity), datasetHomogeneity = max(dataset_homogeneity), globalCompleteness = max(global_completeness), datasetCompletness = max(dataset_completeness), globalcLISI = global_cLISI, datasetcLISI = min(dataset_cLISI), globaliLISI = global_iLISI, datasetiLISI = min(dataset_iLISI), globalcASW = global.casw, datasetcASW = min(dataset.casw), globalbASW = global.basw, datasetbASW = min(dataset.basw))
  
  # Return
  list(summary = summary, cluster_metrics = clustering_metrics, ASW = asw, LISI = lisi, clusters = cluster_list)
}

rna_to_atac_syntax <- function(vec){
  output <- vec %>%
    stringi::stri_replace_last_regex("[_]", "#") %>%
    stringi::stri_join("-1")
  return(output)
}

