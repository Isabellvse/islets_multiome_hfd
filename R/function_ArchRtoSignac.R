## Had to rewrite the ArchR2Signac function from the ArchRtoSignac package. (19/07/22)

#' Convert Archr project to signac (seurat) object
#'
#' @param ArchRProject = an ArchR project
#' @param refversion = The assembly release and versions of UCSC genome reference
#' @param samples = List of all the samples from the ArchRProject
#' @param fragments_dir = path to fragment files
#' @param pm = peak matrix
#' @param output_dir = output directory
#' @param annotation = Gene annotation
#'
#' @return a chromatin assay object
#' @export
#'
#' @examples
ArchR2Signac_multiome <- function (ArchRProject, refversion, samples = NULL, fragments_dir = NULL,
                                   pm, output_dir = "/outs/", annotation)
{
    if (is.null(samples)) {
        samples <- unique(proj@cellColData$Sample)
    }
    print("In Progress:")
    print("Prepare Seurat list for each sample")
    seurat_list <- lapply(samples, function(cur_sample) {
        print(cur_sample)
        cur_fragments <- paste0(fragments_dir, cur_sample, output_dir,
                                "atac_fragments.tsv.gz") ### changed from 'fragments.tsv.gz'
        cur_pm <- pm[, grepl(paste0(cur_sample, "#"), colnames(pm))]
        cur_meta <- ArchRProject@cellColData %>% as.data.frame %>%
            subset(Sample == cur_sample)
        colnames(cur_pm) <- do.call(rbind, str_split(colnames(cur_pm),
                                                     "#"))[, 2]
        rownames(cur_meta) <- do.call(rbind, str_split(rownames(cur_meta),
                                                       "#"))[, 2]
        print(dim(cur_pm))
        cur_chromatin <- Signac::CreateChromatinAssay(counts = cur_pm,
                                                      sep = c("-", "-"), fragments = cur_fragments, ranges = ArchRProject@peakSet,
                                                      genome = refversion, annotation = annotation)
        cur_atac <- Seurat::CreateSeuratObject(cur_chromatin,
                                               assay = "peaks", meta.data = cur_meta, )
    })
    print("In Progress:")
    print("Merge Seurat list")
    SeuratObject <- merge(x = seurat_list[[1]], y = seurat_list[2:length(seurat_list)],
                          add.cell.ids = samples)
    print("Return SeuratObject")
    SeuratObject
}



## Convert archr tile matrix to chromatin assay
ArchR2Signac_tiles <- function (ArchRProject, refversion, samples = NULL, fragments_dir = NULL,
                                pm, output_dir = "/outs/", annotation)
{
    if (is.null(samples)) {
        samples <- unique(proj@cellColData$Sample)
    }
    print("In Progress:")
    print("Prepare Seurat list for each sample")
    seurat_list <- lapply(samples, function(cur_sample) {
        print(cur_sample)
        cur_fragments <- paste0(fragments_dir, cur_sample, output_dir,
                                "atac_fragments.tsv.gz") ### changed from 'fragments.tsv.gz'
        cur_pm <- pm[, grepl(paste0(cur_sample, "#"), colnames(pm))]
        cur_meta <- ArchRProject@cellColData %>% as.data.frame %>%
            subset(Sample == cur_sample)
        colnames(cur_pm) <- do.call(rbind, str_split(colnames(cur_pm),
                                                     "#"))[, 2]
        rownames(cur_meta) <- do.call(rbind, str_split(rownames(cur_meta),
                                                       "#"))[, 2]
        print(dim(cur_pm))
        cur_chromatin <- Signac::CreateChromatinAssay(counts = cur_pm,
                                                      sep = c("-", "-"), fragments = cur_fragments,
                                                      genome = refversion, annotation = annotation)
        cur_atac <- Seurat::CreateSeuratObject(cur_chromatin,
                                               assay = "tiles", meta.data = cur_meta, )
    })
    print("In Progress:")
    print("Merge Seurat list")
    SeuratObject <- merge(x = seurat_list[[1]], y = seurat_list[2:length(seurat_list)],
                          add.cell.ids = samples)
    print("Return SeuratObject")
    SeuratObject
}

#' Extract tile matrix from archr project, with rownames chr-start-end
#'
#' @param archr_proj = an ArchR project
#' @param use_matrix = name of matrix to extract: "TileMatrix"
#' @param bin_size = size of tile matrix bins, e.g. 5000
#'
#' @return a matrix object
#' Author: Isabell V. S. Ernst
create_tile_matrix <- function(archr_proj,
                               use_matrix,
                               bin_size){
    tiles <- ArchR::getMatrixFromProject(
        ArchRProj = archr_proj,
        useMatrix = use_matrix,
        useSeqnames = NULL,
        verbose = TRUE,
        binarize = TRUE,
        threads = 1,
        logFile = createLogFile("getMatrixFromProject"))
    # get matrix
    tiles_matrix <- SummarizedExperiment::assays(tiles)$TileMatrix
    # Row data
    tiles_row <- SummarizedExperiment::rowData(tiles)
    # Create rownames chr-start-end, from: https://github.com/GreenleafLab/ArchR/issues/330
    tiles_row <- tiles_row %>%
        as.data.frame() %>%
        dplyr::mutate(end = (start + bin_size) - 1,
                      rown = paste0(seqnames, "-", start, "-", end))
    # add new rownames to tile matrix
    rownames(tiles_matrix) <- tiles_row$rown
    # return
    return(tiles_matrix)
}
