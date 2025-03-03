# Set library path
.libPaths("/work/Home/islets_multiome/nichenet_library/")

# source files
library(multinichenetr)
library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(SingleCellExperiment)
base::source(here::here("R/functions_multinichenetr.R"))
base::source(here::here("R/functions.R"))
base::source(here::here("R/functions_deseq2.R"))
base::source(here::here("R/misc.R"))
