# Description -------------------------------------------------------------
# Identify high inflammatory beta-cells in data from humans without Diabetes, type 1 diabetes or type 2 diabetes
# Setup -------------------------------------------------------------------
base::source(here::here("R/set_up.R"))
set.seed(1000)

create_directories(c(here::here("data/inflammatory_cohort_2"),
                     here::here("data/inflammatory_cohort_2/files/"),
                     here::here("data/inflammatory_cohort_2/figures/")))

Diabetes_levels <- c("ND", "T2D", "T1D")
# Load --------------------------------------------------------------------
all_nichenet_gene <- base::readRDS(here::here("data/inflammatory_cluster/files/all_nichenet_gene.rds"))
all_immune <- all_nichenet_gene %>%
  stringr::str_replace(pattern = "[.]", replacement = "-") # replace "." with "-" in gene names

## Convert mouse genes to human genes
immune_gene_human <- mouse2human_symbol(all_immune)

# hpap
hpap <- readRDS(here::here("data/public_data/hpap_islet_scRNAseq.rds"))
donor_info <- readxl::read_excel(path = here::here("data/public_data/donor_information.xlsx"), sheet = "donor")

# Preprocess --------------------------------------------------------------
# subset hpap data
Idents(hpap) <- "Cell Type Grouped"
hpap_beta <- subset(hpap, idents = "Beta")

# add Diabetes scale
metadata <- donor_info %>%
  dplyr::select(donor_ID, hba1c, clinical_diagnosis, disease_duration)
metadata$Diabetes <- "ND"
metadata[ metadata$clinical_diagnosis %in% c("T1DM", "T1DM (recent DKA)", "Recent T1DM Unsuspected", "T1DM or MODY, undetermined"), "Diabetes"] <- "T1D"
metadata[ metadata$clinical_diagnosis %in% c("T2DM", "T2DM polycystic ovaries"), "Diabetes"] <- "T2D"
metadata[ metadata$clinical_diagnosis %in% c("T2DM (? prediabetic)"), "Diabetes"] <- "Pre"
metadata[ metadata$clinical_diagnosis %in% c("T2DM LADA?"), "Diabetes"] <- "LADA_T2D"
metadata[ metadata$clinical_diagnosis %in% c("T2DM Gastric bypass"), "Diabetes"] <- "Gastic_Bypass_T2D"
metadata <- metadata %>% dplyr::rename(Library = donor_ID)

# add metadata to seurat object
hpap_beta@meta.data <- hpap_beta@meta.data %>%
  BiocGenerics::as.data.frame() %>%
  tibble::rownames_to_column("cell") %>%
  dplyr::left_join(y = metadata, by = "Library") %>%
  tibble::column_to_rownames("cell")

# Module score ------------------------------------------------------------
hpap_beta  <- hpap_beta  %>%
  UCell::AddModuleScore_UCell(features = base::list("tnf_ifnb1" = immune_gene_human),
                              name = "_ucell",
                              maxRank = 1000,
                              BPPARAM = MulticoreParam(workers = parallel::detectCores() - 1))
# Plot distribution of module scores --------------------------------------
hpap_beta@meta.data <- hpap_beta@meta.data %>%
  dplyr::rename(donor = Library)


p <- hpap_beta@meta.data %>% 
  dplyr::group_by(donor, Diabetes) %>%
  dplyr::summarize(mean_score = base::mean(tnf_ifnb1_ucell)) %>% 
  ggplot2::ggplot(aes(x = Diabetes, y = mean_score)) + 
  ggplot2::geom_point(size = 0.1) +
  ggplot2::geom_boxplot(outlier.shape = NA,
                        width = 0.5, fill = "transparent", color = "black") +
  ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
  ggplot2::expand_limits(y = 0) +
  ggplot2::labs(y = "Mean UCell activity score") +
  my_theme()

ggsave(here::here("data/inflammatory_cohort_2/figures/tnf_ifnb1_ucell_score.pdf"),
       p,
       width = 1, height = 1)


# Statistical test --------------------------------------------------------
## Mean module score ----
test <- hpap_beta@meta.data %>%
  dplyr::group_by(donor, Diabetes) %>%
  dplyr::summarize(n_total = n(),
                   mean_score = base::mean(tnf_ifnb1_ucell)) %>%
  dplyr::ungroup() %>%
  rstatix::wilcox_test(mean_score ~ Diabetes)

eff <-hpap_beta@meta.data %>%
  dplyr::group_by(donor, Diabetes) %>%
  dplyr::summarize(n_total = n(),
                   mean_score = base::mean(tnf_ifnb1_ucell)) %>%
  dplyr::ungroup() %>%
  rstatix::wilcox_effsize(mean_score ~ Diabetes)

res <- dplyr::full_join(test, eff)

openxlsx::write.xlsx(res, here::here("data/inflammatory_cohort_2/files/wilcox_test_inflammatory_genes_module_score.xlsx"))

# Number of cells with an inflammatory score above contorl median ---------
# Find the 75th percentile of mean inflammatory score in nondiabetic beta-cells
control_upper_q <- hpap_beta@meta.data %>%
  dplyr::filter(Diabetes == "ND") %>%
  dplyr::select(donor, tnf_ifnb1_ucell) %>%
  dplyr::summarize(uppr_quartile = stats::quantile(tnf_ifnb1_ucell, na.rm = TRUE, probs = 0.75))

# Calculate percentage of beta-cell above this threshold
n_df <- hpap_beta@meta.data %>%
  dplyr::select(donor, tnf_ifnb1_ucell, Diabetes) %>%
  dplyr::group_by(donor, Diabetes) %>%
  dplyr::summarise(n_above = base::sum(tnf_ifnb1_ucell  > control_upper_q[["uppr_quartile"]]),
                   n_total = n(),
                   frac_inf = (n_above/n_total)*100)
# plot it
p <- n_df %>%
  ggplot2::ggplot(aes(x = Diabetes, y = frac_inf)) +
  ggplot2::geom_boxplot(outlier.shape = NA) +
  ggplot2::geom_point(color = "black", size = 0.1) +
  ggplot2::expand_limits(y = 0) +
  ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 7)) +
  ggplot2::labs(y = "% above 75th percentile\n(ND UCell score)",
                title = "Cohort 2") +
  my_theme() +
  ggplot2::theme(legend.position = "none")

ggsave(here::here("data/inflammatory_cohort_2/figures/fraction_high_tnf_ifnb1_ucell_score.pdf"),
       p,
       width = 1, height = 1)

# Statitical test
test <- n_df %>%
  dplyr::select(Diabetes, frac_inf) %>%
  dplyr::ungroup() %>%
  rstatix::wilcox_test(frac_inf ~ Diabetes)
# effect size
eff <- n_df %>%
  dplyr::ungroup() %>%
  rstatix::wilcox_effsize(frac_inf ~ Diabetes)
# combine
res <- dplyr::full_join(test, eff)
openxlsx::write.xlsx(res, here::here("data/inflammatory_cohort_2/files/wilcox_test_fraction_of_high_inflammatory_ucell.xlsx"))

# Save seurat object ------------------------------------------------------
base::saveRDS(hpap_beta, here::here("data/inflammatory_cohort_2/files/hpap_beta.rds"))
