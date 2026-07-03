install_if_missing <- function(packages) {
  missing <- packages[!packages %in% rownames(installed.packages())]
  
  if (length(missing) > 0) {
    install.packages(missing)
  } else {
    cat("All packages already installed!\n")
  }
}

install_if_missing(c(
  "rprojroot", "patchwork", "inflection", "gridExtra", "lmerTest", "lme4", 
  "emmeans", "car", "ggrastr", "khroma", "ggpubr", "spatstat", "coin", 
  "magick", "future", "furrr", "concaveman", "jsonlite", "rstatix", "gridExtra", "ggbreak", "glmmTMB", "broom.mixed"
))

library(furrr)
library(dplyr)
library(ggplot2)
library(purrr)
library(stringr)
library(patchwork)
library(lmerTest)
library(lme4)
library(ggrastr)
library(spatstat)
library(concaveman)
library(jsonlite)
library(rstatix)

source(here::here("R/functions_revisions.R"))