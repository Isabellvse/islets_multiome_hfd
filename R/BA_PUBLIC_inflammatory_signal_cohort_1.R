# Description -------------------------------------------------------------
# Identify high inflammatory beta-cells in data from humans without disease, predisease or type 2 disease
# Setup -------------------------------------------------------------------
base::source(here::here("R/set_up.R"))
set.seed(1000)

create_directories(c(here::here("data/inflammatory_cohort_1"),
                     here::here("data/inflammatory_cohort_1/files/"),
                     here::here("data/inflammatory_cohort_1/figures/")))

disease_levels <- c("ND", "Pre", "T2D")
# Load --------------------------------------------------------------------
# Download ----------------------------------------------------------------
# https://cellxgene.cziscience.com/collections/58e85c2f-d52e-4c19-8393-b854b84d516e
# in terminal: https://cellxgene.cziscience.com/collections/58e85c2f-d52e-4c19-8393-b854b84d516e
# I downlaoded the "Reintegrated transcriptomes of Beta pancreatic islet cells" 99,029 cells
seu_beta <- base::readRDS(here::here("data/inflammatory_cohort_1/files/motakis_beta.rds"))
seu_beta@meta.data <- seu_beta@meta.data %>% 
  BiocGenerics::as.data.frame() %>% 
  dplyr::mutate(disease = dplyr::case_when(disease == "normal" ~ "nd",
                                           disease == "prediabetes syndrome" ~ "pre",
                                           disease == "type 2 diabetes mellitus" ~ "t2d"),
                disease = base::factor(disease, levels = c("nd", "pre", "t2d")))

all_nichenet_gene <-base:: readRDS(here::here("data/inflammatory_cluster/files/all_nichenet_gene.rds"))
all_immune <- all_nichenet_gene %>%
  stringr::str_replace(pattern = "[.]", replacement = "-") # replace "." with "-" in gene names

## Convert mouse genes to human genes
immune_gene_human <- mouse2human_ensembl(all_immune)

# Module score ------------------------------------------------------------
seu_beta  <- seu_beta  %>%
  UCell::AddModuleScore_UCell(features = base::list("tnf_ifnb1" = immune_gene_human),
                              name = "_ucell",
                              maxRank = 1000,
                              BPPARAM = MulticoreParam(workers = parallel::detectCores() - 1))

# Plot distribution of module scores --------------------------------------
p <- seu_beta@meta.data %>% 
  dplyr::group_by(donor_id, disease) %>%
  dplyr::summarize(mean_score = base::mean(tnf_ifnb1_ucell)) %>% 
  ggplot2::ggplot(aes(x = disease, y = mean_score)) + 
  ggplot2::geom_point(size = 0.1) +
  ggplot2::geom_boxplot(outlier.shape = NA,
                        width = 0.5, fill = "transparent", color = "black") +
  ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
  ggplot2::expand_limits(y = 0) +
  ggplot2::labs(y = "Mean UCell activity score") +
  my_theme()

ggsave(here::here("data/inflammatory_cohort_1/figures/tnf_ifnb1_ucell_score.pdf"),
       p,
       width = 1, height = 1)


# Statistical test --------------------------------------------------------
## Mean module score ----
test <- seu_beta@meta.data %>%
  dplyr::group_by(donor_id, disease) %>%
  dplyr::summarize(n_total = n(),
                   mean_score = base::mean(tnf_ifnb1_ucell)) %>%
  dplyr::ungroup() %>%
  rstatix::wilcox_test(mean_score ~ disease)

eff <-seu_beta@meta.data %>%
  dplyr::group_by(donor_id, disease) %>%
  dplyr::summarize(n_total = n(),
                   mean_score = base::mean(tnf_ifnb1_ucell)) %>%
  dplyr::ungroup() %>%
  rstatix::wilcox_effsize(mean_score ~ disease)

res <- dplyr::full_join(test, eff)

openxlsx::write.xlsx(res, here::here("data/inflammatory_cohort_1/files/wilcox_test_inflammatory_genes_module_score.xlsx"))

# Number of cells with an inflammatory score above contorl median ---------
# Find the 75th percentile of mean inflammatory score in nondiabetic beta-cells
control_upper_q <- seu_beta@meta.data %>%
  dplyr::filter(disease == "nd") %>%
  dplyr::select(donor_id, tnf_ifnb1_ucell) %>%
  dplyr::summarize(uppr_quartile = stats::quantile(tnf_ifnb1_ucell, na.rm = TRUE, probs = 0.75))

# Calculate percentage of beta-cell above this threshold
n_df <- seu_beta@meta.data %>%
  dplyr::select(donor_id, tnf_ifnb1_ucell, disease) %>%
  dplyr::group_by(donor_id, disease) %>%
  dplyr::summarise(n_above = base::sum(tnf_ifnb1_ucell  > control_upper_q[["uppr_quartile"]]),
                   n_total = n(),
                   frac_inf = (n_above/n_total)*100)
# plot it
p <- n_df %>%
  ggplot2::ggplot(aes(x = disease, y = frac_inf)) +
  ggplot2::geom_boxplot(outlier.shape = NA) +
  ggplot2::geom_point(color = "black", size = 0.1) +
  ggplot2::expand_limits(y = 0) +
  ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 7)) +
  ggplot2::labs(y = "% above 75th percentile\n(ND UCell score)",
                title = "Cohort 1") +
  my_theme() +
  ggplot2::theme(legend.position = "none")

ggsave(here::here("data/inflammatory_cohort_1/figures/fraction_high_tnf_ifnb1_ucell_score.pdf"),
       p,
       width = 1, height = 1)

# Statitical test
test <- n_df %>%
  dplyr::select(disease, frac_inf) %>%
  dplyr::ungroup() %>%
  rstatix::wilcox_test(frac_inf ~ disease)
# effect size
eff <- n_df %>%
  dplyr::ungroup() %>%
  rstatix::wilcox_effsize(frac_inf ~ disease)
# combine
res <- dplyr::full_join(test, eff)
openxlsx::write.xlsx(res, here::here("data/inflammatory_cohort_1/files/wilcox_test_fraction_of_high_inflammatory_ucell.xlsx"))

# Save seurat object ------------------------------------------------------
base::saveRDS(seu_beta, here::here("data/inflammatory_cohort_1/files/seu_beta.rds"))
