# Description -------------------------------------------------------------
# Get geneset from islets that correlate with various clinical phenotypes from obese F2 mice.
# This data is from the ATTIE LAB DIABETES DATABASE: http://diabetes.wisc.edu/index.php

# Gene expression from islets of langerhans were correlated with a clinical trial with a
# minimum correlation of 0.1, each result was copied manually into an excel document
# correlation databases can be found here: http://diabetes.wisc.edu/correl_f2.php

# *Description of the correlation search tool from the website:
# "Allows for gene-gene or gene-clinical trait correlation analysis among all F2 mice in each of 6 tissues
# (islet, liver, adipose, muscle, hypothalamus, kidney).
# Can also search for genes within a particular tissue that correlate with specific clinical traits.
# Correlated genes lists can be downloaded or easily transferred to DAVID for Gene Ontology enrichment analysis."

# *Description of F2 mice from the website
# " The Genetic study of F2 cohort database derives from the construction of
# ~500 B6:BTBR F2 mice that were all sacrificed at the same age (10 weeks) and were all obese,
# thereby focusing the study on genetic differences between the parental strains.
# All F2 mice were genotyped at >2,000 informative SNPs,
# allowing us to link genotypic variation to phenotypic variation.
# Gene expression was profiled in 6 tissues (islet, liver, adipose, hypothalamus, gastrocnemius and kidney) in every F2 mouse.
# In addition, >100 diabetes-related clinical phenotypes were measured in all mice (e.g., plasma insulin).
# Tools are provided to determine genetic linkage for expression or clinical traits;
# eQTL and cQTL. Correlation analysis among genes or between genes and clinical phenotypes is also available. "

# supplementary table 9
# Setup -------------------------------------------------------------------
base::source(here::here("R/set_up.R"))
set.seed(1000)

# Load --------------
path = here::here("data/attie_diabetes_database/genomic_study_f2_cohort/correlation/Correlations_extended.xlsx")

attie_cor <- path %>%
  readxl::excel_sheets() %>%
  rlang::set_names() %>%
  purrr::map(readxl::read_excel,
             path = path,
             col_names = c("tissue", "transcript", "correl_coeff"),
             col_types = c("text", "text", "numeric"),
             skip = 2)
# Save -----
saveRDS(attie_cor, here::here("data/attie_diabetes_database/genomic_study_f2_cohort/correlation/Cor_list_extended.rds"))
