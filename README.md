# Unveiling Early Responses of Mouse β-Cells to High Fat Diet feeding through Single-Cell Multiome Analysis

To investigate early adaptive changes to short-term HFD, we fed eleven-week-old male C57BL/6JBomTac mice with a HFD (60% energy from fat) or purified low-fat diet (LFD) (10% energy from fat) as a control, for one and three weeks.

To capture the early genomic changes within islets of Langerhans following short-term HFD exposure, we utilized a single-nucleus method known as chromium single-cell multiome ATAC + Gene expression from 10X Genomics (snMultiome) using nuclei from isolated islets. snMultiome enables simultaneous capture of gene expression data and chromatin accessibility data from the same single nuclei through RNA-sequencing (snRNA-seq) and transposase-accessible chromatin with sequencing (snATAC-seq), respectively.

![Illustration of work](https://github.com/Isabellvse/islets_multiome/blob/main/illustrations/github_image.png)

###### *A) Schematic overview of the experimental setup: 8-week-old male C57BL/6JBomTac mice underwent a three-week acclimatization period during which they were fed a low-fat-diet (LFD). Subsequently, mice were divided into groups receiving either LFD (orange, as a control) or high-fat-diet (blue, HFD) for one or three weeks.* 
###### *B) Macronutrient composition of LFD (7 % energy from sucrose 10 % energy from fat) and HFD (7 % energy from sucrose, 60 % energy from fat) presented in % Kcal.* 
###### C) *Illustration of proposed mechanism: Following one week of HFD (HFD 1 wk) feeding UPR-mediated stress and transcriptional downregulation of β-cell identity genes in most β-cells, as well as induction of an inflammatory response in a small subset of β-cells (orange cells). After three weeks of HFD (HFD 3 wk), the UPR stress is resolved, but the subset of β-cells with an inflammatory response is increased. In β-cells without an inflammatory response, β-cell identity gene transcription is partially re-stored, but not in β-cells with an inflammatory response. Read arrows indicate upregulation and blue arrows indicate downregulation.*


# Brief description of folder

The following folders contain:

-   `preprocessing/`: This folder contains scripts for alignment of snRNA and snATAC-seq data
-   `R/`: This folder contains scripts for quality control and analysis performed in R and python
-   `libraries/`: Different versions of packages are used in this analysis (indicated in each script). This file contains session info for each library.
-   `illustrations/`: Illustrations used for this readme file

# Workflow

![Overview of analysis workflow](https://github.com/Isabellvse/islets_multiome_hfd/blob/main/illustrations/workflow.png)

###### *Illustration of workflow, and what scripts are used at each step. Illustration was created with [Excalidraw](https://excalidraw.com/)*


# Description of scripts

### Alignment in terminal

-   `preprocessing/`: *Package versions: libraries/alignment.txt*
    -   `1_Download_files.sh`: Download files used for alignment
    -   `2_Generating_genome.sh`: Generate reference genome for STARsolo
    -   `3_fastqc.sh`: Quality control of fastq files
    -   `4_Alignment_RNA.sh`: Alignment of snRNA data with STARsolo both introns and exons
        -   `alignment_to_gene.sh`: Alignment of snRNA data with STARsolo just exons
    -   `5_cellranger_arc_processing.sh`: Alignment of snATAC-seq data with cellranger arc

### Quality control and find clusters in R

-   `R/`: *Package versions: libraries/r_library.txt*
    -   `A_RNA_quality_control.R`: RNA data quality control -- remove empty droplets and other quality control using Validrops, create seurat object, identify doublets using DoubletFinder.
    -   `AB_ATAC_quality_control.R`: ATAC data quality control -- create ArchR project quality control and doublet identify doublets using ArchR
    -   `AC_RNA_ATAC_quality_control.R`: Identify barcodes which pass quality control in both assays, and also cells marked as doublets both by DoubletFinder and ArchR.
    -   `AD_RNA_polyhormal_doublet_cells.R`: Using high quality nuclei from script 3 (passed quality control in both assays, not marked as doublet in both assays), cluster nuclei based on RNA data, and integrate data using JOINTLY, and remove clusters which have polyhormone expression of *Ins1*, *Sst*, *Ppy* and *Gcg*. Now the RNA data is integrated, and we only have high quality nuclei.
    -   `AE_ATAC_keep_good_quality_nuclei_and_cluster.R`: For the ATAC data, keep only nuclei which are of high quality (identified in 4_RNA_polyhormal_doublet_cells.R)
    -   `AF_ATAC_integration_tiles.R`: For the ATAC data, integrate tile data using LIGER
    -   `AG_ATAC_determine_cluster.R`: For the ATAC data, find clusters using the integrated data found in AH_ATAC_integration_tiles.R, and annotate clusters using marker gene activity scores and RNA gene expression.
    -   `AI_ATAC_find_peaks.R`: For ATAC data, using the clusters identified in 7_ATAC_determine_cluster.R, find reproducible peaks.
    -   `AJ_ATAC_integration_peaks.R`: For ATAC data, same workflow as 7_ATAC_determine_cluster.R, but integrate data uing peaks identified in 8_ATAC_find_peaks.R.
    -   `AK_RNA_ATAC_combine_peaks_and_genes.R`: Combine RNA and ATAC data in a seurat object, cluster ATAC-data, and perform weighted nearest neighbor analysis (RNA+ATAC clustering)
    -   `AL_RNA_ATAC_cell_type_annotation_manual.R`: Annotate cells using marker genes expression and activity, and also make some plots of cell type distributions.

### Celltype prioritization and differential gene expression analysis (RNA) in R

-   `R/`: *Package versions: libraries/r_library.txt*
    -   `AL_cell_type_prioritization.R`: Explore which cell types is most peturbed by diet using Augur
    -   `AM_differential_expression_analysis_pseudobulk.R`: Differential gene expression analysis across conditions using DESeq2 and pseudobulk counts

### Cell-Cell communication analysis (RNA) in R

-   `R/`: *Package versions: libraries/r_nichenet_library.txt*
    -   `AN_cell_cell_communication.R`: cell-cell communication analysis using MultinichenetR

### Identifying enhancer driven eRegulons with SCENIC+ (RNA and ATAC)

#### Preparing data in R

-   `R/`: *Package versions: libraries/r_library.txt*
    -   `AO_BETA_preparedata_for_scenicplus.R`: Create AnnData object and extract fragments from seurat object as input for SCENIC+

#### Run SCENIC+ analysis in phyton

-   `R/`: *Package versions: libraries/scenic_plus.txt*
    -   `AO_BETA_scenic_plus.sh`: Run the entire SCENIC+ pipeline
    -   `AQ_BETA_scenicplus_export_data.sh`: Export data from the scenicplus object, such as eRegulon meta data, and transcription factor motifs associated with each eRegulon

#### eRegulon characterization

-   `R/`: *Package versions: libraries/r_library.txt*
    -   `AR_BETA_get_scenicplus_motif_pwm.R`: Genereate motif position weight metrics (pwm) list from pwms used in scenicplus
    -   `AS_BETA_chromvar_analysis.R`: Calculate motif activity using ChromVar using motifs from scenicplus, retrived in script 19
    -   `AT_BETA_eregulon_filter.R`: Filter high quality eRegulons based on the product of transcription factor gene expression and motif activity correlation to AUC score of eRegulon target regions and target genes
    -   `AU_ATAC_BETA_eregulon_heatmaps.R`: Visualize transcription factor expression, motif activity and AUC scores for high quality eRegulons in a heatmap and dotplot
    -   `AV_prepare_geneset_attie_diabetes_database.R`: Prepare genesets from the Attie Lab diabetes database
    -   `AW_BETA_eregulon_attie_genetset_GSEA.R`: Perform GSEA of genes from attie diabetes database and eREgulon target genes
    -   `AX_BETA_cell_communication_target_genes_and_eRegulon_genes.R`: Overrepresentation analysis (ORA) with eRegulon target genes, and putative tnf and ifnb1 target genes (found in script 14_cell_cell_communication.R)

#### β-cell subpopulation

-   `R/`: *Package versions: libraries/r_library.txt*
    -   `AY_PUBLIC_gluco_lipo_cyto_treatment`: Prepare public scRNA-seq data from human islets treated with various stressors for analysis, and calculate gene module activity scores of gene clusters identified in 16_differential_expression_analysis_pseudobulk.R within these human β-cells.
    -   `AZ_BETA_inflammatory_subcluster.R`: Generate an inflammatory gene set of putative tnf and ifnb1 target genes (found in script 14_cell_cell_communication.R), and catagorize β-cells into a high inflammatory subpopulation or low inflammatory subpopulation based on the gene module activity score of the inflammatory geneset with the β-cell, and find differentially expressed genes between subtypes and more. 
    -   `BA_PUBLIC_inflammatory_signal_cohort_1` and `BB_PUBLIC_inflammatory_signal_cohort_2`: For public data, calculate gene module activity scores with the same genes as in AX_BETA_inflammatory_subcluster.R, in human beta-cells from non-diabetes, prediabetic, type 1 diabetic and type 2 diabeteic individuals. 

#### Other scripts

-   `R/`:
    -   `BC_exporting_tables.R`: code to export and set up various result tables
    -   `BD_cexporting_seurat_objects.R`: Save seurat object as h5ad
    -   `Package_load.R`: to load packages
    -   `misc.R`: colors and levels used for plotting
    -   `functions.R`: functions for tasks
    -   `functions_deseq2.R`: functions used for DEseq2 analysis and plotting
    -   `function_ArchRtoSignac.R`: function used to transfer data from ArchR to seurat

##### Documents
`meta.xlsx`: meta data for each sample
