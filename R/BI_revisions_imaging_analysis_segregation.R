# Description -------------------------------------------------------------
# TEst if STAT1 high and STAT1 low cells are spatially segregated

# Setup -------------------------------------------------------------------
source(here::here("islets_multiome_hfd/R/functions.R"))
#source(here::here("islets_multiome_hfd/R/prepare data.R"))
create_directories(here::here("data/revisions/percentile/plots"))
create_directories(here::here("data/revisions/percentile/files"))
set.seed(100)

# Load --------------------------------------------------------------------
cell <- vroom::vroom(here::here("data/revisions/percentile/files/cell_data_status.csv")) |> 
  dplyr::mutate(diet = factor(diet, diet_lvl))
islet <- vroom::vroom(here::here("data/revisions/spatial_lfd_hfd/islet_data.csv")) |> 
  dplyr::mutate(diet = factor(diet, diet_lvl))
window_list <- readRDS(here::here("data/revisions/spatial_lfd_hfd/window_list.rds"))

# Preprocess --------------------------------------------------------------
# Create list of status
per_cell_or <- cell |> 
  dplyr::select(unique_islet, unique_cell, mouse_id,
                centroid_x_px, centroid_y_px, distance_from_edge_px, 
                diet, mouse_id, ends_with("or") & contains("stat1_status"), stat1_nucleus_mean, stat1_peri_mean) |> 
  tidyr::pivot_longer(ends_with("_or"), names_to = "percentiles", values_to = "stat1_status") |> 
  split(~percentiles)

# Calculate percentage of high STAT1 cells
per_islet <- per_cell_or |>
  purrr::list_rbind(names_to = "percentiles") |>
  dplyr::group_by(percentiles, diet, unique_islet, stat1_status) |>
  dplyr::summarise(
    cells = n_distinct(unique_cell),
    .groups = "drop"
  ) |>
  tidyr::pivot_wider(
    names_from = stat1_status,
    values_from = cells,
    values_fill = 0
  ) |>
  dplyr::rename(cells_high = high, cells_low = low) |> 
  dplyr::group_by(percentiles, diet, unique_islet) |>
  dplyr::mutate(
    total_cells = cells_high + cells_low,
    pct_high = round((cells_high / total_cells) * 100, 2)
  ) |>
  dplyr::ungroup() |>
  dplyr::select(percentiles, diet, unique_islet, total_cells, pct_high, cells_high, cells_low)


# Distribution of pct high islets -----------------------------------------
pdf(here::here("data/revisions/percentile/plots/stat1_high_islet_hfd_histogram.pdf"), width = 4, height = 1.5)
per_islet |> 
  dplyr::filter(diet == "hfd") |> 
  dplyr::mutate(
    percentiles = stringr::str_remove(percentiles, "stat1_status_"),
    percentiles = stringr::str_remove(percentiles, "_or")
  ) |> 
  ggplot2::ggplot(ggplot2::aes(x = pct_high)) +
  ggplot2::geom_histogram(fill = "#fcb21c") +
  ggplot2::facet_wrap(~percentiles, scales = "free_y", nrow = 1) +
  ggplot2::labs(x = "% of STAT1 high cells per islet") +
  ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n =1)) +
  ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n =5)) +
  my_theme()
dev.off()
# Segregation test --------------------------------------------------------
# Are High and low STAT1 cells in HFD islets spatially segregated or random?

plan(multisession, workers = 30)

results <- per_cell_or |> 
  purrr::modify_depth(1, \(df){
    df |> 
      dplyr::filter(diet == "hfd") |> 
      split(~unique_islet) |> 
      furrr::future_imap(\(data, name){
        segregation_test(islet_data = data, islet_name = name)
      }, 
      .options = furrr_options(seed = TRUE)) |> 
      purrr::compact()
  })

plan(sequential)

# Extract results ---------------------------------------------------------
# Get sigma used
seg_sigma <- results  |> 
  purrr::modify_depth(1, \(vec){vec |> 
      purrr::map("sigma")}) |> 
  purrr::modify_depth(2, \(vec){data.frame(sigma_x = vec[1], 
                                           sigma_y = vec[2], 
                                           row.names = NULL)}) |> 
  purrr::modify_depth(1, \(df){df |> 
      purrr::list_rbind(names_to = "unique_islet")})|> 
  purrr::list_rbind(names_to = "percentiles")

seg_results <- results  |> 
  purrr::modify_depth(1, \(df){df |> 
      purrr::map("results") |> 
      purrr::list_rbind(names_to = "unique_islet") |> 
      dplyr::mutate(fdr = stats::p.adjust(p.value, "fdr"))}) |> 
  purrr::list_rbind(names_to = "percentiles") |> 
  dplyr::left_join(seg_sigma) |> 
  dplyr::left_join(per_islet)

# save --------------------------------------------------------------------
saveRDS(results, here::here("data/revisions/percentile/files/segregation_test.rds"))
vroom::vroom_write(seg_results, here::here("data/revisions/percentile/files/segregation_test_results.csv"))

# Load again --------------------------------------------------------------
results <- readRDS(here::here("data/revisions/percentile/files/segregation_test.rds"))
seg_results <- vroom::vroom(here::here("data/revisions/percentile/files/segregation_test_results.csv"))

# Estimates the spatially-varying probability -----------------------------
vr_results <- 
  seg_results |> 
  split(~percentiles) |> 
  purrr::map(~split(.x, ~unique_islet)) |> 
  purrr::map(purrr::map, ~with(.x, c(sigma_x, sigma_y))) |> 
  purrr::imap(\(islets, percentile) {
    purrr::imap(islets, \(sigma, islet) {
      per_cell_or[[percentile]] |>
        dplyr::filter(unique_islet == islet) |>
        relrisk_test(islet, sigma)
    })
  })

# Save --------------------------------------------------------------------
saveRDS(vr_results, here::here("data/revisions/percentile/files/relative_riskt.rds"))
vr_results <- readRDS(here::here("data/revisions/percentile/files/relative_riskt.rds"))

# QC segregation test -----------------------------------------------------
## Are significant islets imbalanced by high/low ratio of cells
pdf(here::here("data/revisions/percentile/plots/segregation_test_high_low_ratio.pdf"), width = 8, height = 3)
seg_results |>
  dplyr::mutate(
    ratio_high_low = cells_high / cells_low,
    significant = fdr < 0.05
  ) |>
  ggplot2::ggplot(ggplot2::aes(x = ratio_high_low, fill = significant)) +
  ggplot2::geom_histogram(alpha = 0.7) +
  ggplot2::facet_wrap(~percentiles, nrow = 1) +
  ggplot2::scale_x_log10() +
  ggplot2::scale_fill_manual(values = c("blue", "red")) +
  ggplot2::labs(title = "Ratio of high / low STAT1 in islets") +
  my_theme()
dev.off()

# Summary of segregation test ---------------------------------------------
seg_summary <- seg_results |>
  dplyr::mutate(
    significant = fdr < 0.05,
    percentiles = stringr::str_remove(percentiles, "stat1_status_"),
    percentiles = stringr::str_remove(percentiles, "_or")
  ) |>
  dplyr::group_by(percentiles, significant) |>
  dplyr::summarise(
    n_islets = n(),
    mean_statistic = round(mean(statistic), 2),
    mean_pvalue = round(mean(p.value), 4),
    mean_fdr = round(mean(fdr), 4),
    mean_high = round(mean(cells_high), 1),
    mean_low = round(mean(cells_low), 1),
    mean_total = round(mean(total_cells), 1),
    mean_pct_high = round(mean(pct_high), 1),
    median_statistic = round(median(statistic), 2),
    median_pvalue = round(median(p.value), 4),
    median_fdr = round(median(fdr), 4),
    median_high = round(median(cells_high), 1),
    median_low = round(median(cells_low), 1),
    median_total = round(median(total_cells), 1),
    median_pct_high = round(median(pct_high), 1),
    .groups = "drop"
  )

# Save --------------------------------------------------------------------
vroom::vroom_write(seg_summary, here::here("data/revisions/percentile/files/segregation_test_summary.csv"))

# Plot results (islets)------------------------------------------------------------
# Create lookup table for significance
is_sig_lookup <- seg_results |>
  dplyr::mutate(is_sig = fdr < 0.05) |> 
  dplyr::select(percentiles, unique_islet, is_sig)

# Plot based on if results are significant or not
plot_list <- purrr::imap(vr_results, \(percentile_islets, percentile) {
  purrr::imap(percentile_islets, \(vr_obj, islet) {
    is_sig <- is_sig_lookup |>
      dplyr::filter(percentiles == percentile, unique_islet == islet) |>
      dplyr::pull(is_sig)
    
    plot_relrisk_ggplot_void_sig(vr_obj$vr, is_significant = is_sig)
  })
})

# Plot based on if results are significant or not with legend
plot_list_legend <- purrr::imap(vr_results, \(percentile_islets, percentile) {
  purrr::imap(percentile_islets, \(vr_obj, islet) {
    is_sig <- is_sig_lookup |>
      dplyr::filter(percentiles == percentile, unique_islet == islet) |>
      dplyr::pull(is_sig)
    
    plot_relrisk_ggplot_void_sig_legend(vr_obj$vr, is_significant = is_sig)
  })
})


plot_list_pp <- purrr::imap(vr_results, \(percentile_islets, percentile) {
  purrr::map(percentile_islets, \(vr_obj) {
    plot_pp_ggplot_void(pp_obj = vr_obj$pp)
  })
})


# Save plots --------------------------------------------------------------
## pp
for (name in names(plot_list_pp)) {
  pdf(here::here(paste0("data/revisions/percentile/plots/pp_plots_high_low_", name, ".pdf")))
  print(plot_ordered_by_stat(plot_list_pp, seg_results, name))
  dev.off()
}

# relrisk
for (name in names(plot_list)) {
  pdf(here::here(paste0("data/revisions/percentile/plots/relative_risk_plots_high_low_", name, ".pdf")))
  print(plot_ordered_by_stat(plot_list, seg_results, name))
  dev.off()
}
# relrisk with legend
pdf(here::here("data/revisions/percentile/plots/relative_risk_plots_high_low_legend.pdf"))
plot_list_legend[["stat1_status_p95_or"]][["hfd_0617_2_Seq0131_1"]]
dev.off()



# Overall descriptives ----------------------------------------------------
# Of islets that has more than 1 identity, how many islets are not randomly distributed?
pdf(here::here("data/revisions/percentile/plots/pct_segregated_islets.pdf"), 
    width = 1.5, height = 2.5)
seg_results |> 
  dplyr::mutate(seg = dplyr::case_when(fdr < 0.05 ~ "segregated",
                                       .default = "random"),
                mouse_id = str_extract(unique_islet, "\\d{4}")) |> 
  dplyr::group_by(mouse_id, percentiles, seg) |>
  dplyr::summarise(
    islets = n_distinct(unique_islet),
    .groups = "drop"
  ) |>
  tidyr::pivot_wider(
    names_from = seg,
    values_from = islets,
    values_fill = 0
  ) |>
  dplyr::group_by(mouse_id, percentiles) |>
  dplyr::mutate(
    total_islets = random + segregated,
    pct_seg = round((segregated / total_islets) * 100, 2),
    percentiles = stringr::str_remove(percentiles, "stat1_status_"),
    percentiles = stringr::str_remove(percentiles, "_or")
  ) |>
  dplyr::ungroup() |> 
  ggplot2::ggplot(ggplot2::aes(x = percentiles, y = pct_seg)) +
  ggplot2::geom_bar(stat = "summary", fun = "mean", fill = "#fcb21c") +
  ggplot2::stat_summary(fun.data = mean_sdl, geom = "errorbar", width = 0.2,
                        linewidth = 0.5) +
  ggplot2::geom_point(size = 2) +
  ggplot2::labs(x = "Percentile threshold",
                y = "Mean (%)",
                title = "Segregated islet across\npercentile thresholds",
                subtitle = "Per mouse") +
  my_theme()
dev.off()

# how many of them ARE randomly distributed?
pdf(here::here("data/revisions/percentile/plots/pct_nonsegregated_islets.pdf"), 
    width = 1.5, height = 1.5)
seg_results |> 
  dplyr::mutate(seg = dplyr::case_when(fdr < 0.05 ~ "segregated",
                                       .default = "random"),
                mouse_id = str_extract(unique_islet, "\\d{4}")) |> 
  dplyr::group_by(mouse_id, percentiles, seg) |>
  dplyr::summarise(
    islets = n_distinct(unique_islet),
    .groups = "drop"
  ) |>
  tidyr::pivot_wider(
    names_from = seg,
    values_from = islets,
    values_fill = 0
  ) |>
  dplyr::group_by(mouse_id, percentiles) |>
  dplyr::mutate(
    total_islets = random + segregated,
    pct_ran = round((random/ total_islets) * 100, 2),
    percentiles = stringr::str_remove(percentiles, "stat1_status_"),
    percentiles = stringr::str_remove(percentiles, "_or")
  ) |>
  dplyr::ungroup() |> 
  ggplot2::ggplot(ggplot2::aes(x = percentiles, y = pct_ran)) +
  ggplot2::geom_bar(stat = "summary", fun = "mean", fill = "#fcb21c") +
  ggplot2::stat_summary(fun.data = mean_sdl, geom = "errorbar", width = 0.2,
                        linewidth = 0.5) +
  ggplot2::geom_point(size = 1, alpha = 0.5) +
  ggplot2::labs(x = "Percentile threshold",
                y = "Mean (%)") +
  my_theme()
dev.off()


# Plot per islet ----------------------------------------------------------
seg_results_wide <- seg_results |> 
  dplyr::select(unique_islet, percentiles, fdr) |> 
  tidyr::pivot_wider(
    names_from = percentiles,
    values_from = fdr,
    values_fill = 100) |> 
  tibble::column_to_rownames("unique_islet")

row_dist <- dist(seg_results_wide, method = "euclidean")
row_hclust <- hclust(row_dist, method = "complete")
row_order <- rownames(seg_results_wide)[row_hclust$order]

pdf(here::here("data/revisions/percentile/plots/segregated_islet_heatmap.pdf"), width = 1.8, height = 3)
seg_results |> 
  dplyr::select(unique_islet, percentiles, fdr) |> 
  dplyr::mutate(
    percentiles = stringr::str_remove(percentiles, "stat1_status_"),
    percentiles = stringr::str_remove(percentiles, "_or")
  ) |> 
  tidyr::pivot_wider(
    names_from = percentiles,
    values_from = fdr,
    values_fill = NA) |> 
  tidyr::pivot_longer(-c(unique_islet), names_to = "percentiles",
                      values_to = "fdr") |> 
  dplyr::mutate(
    unique_islet = factor(unique_islet, levels = row_order),
    significance = dplyr::case_when(
      is.na(fdr) ~ "NA",
      fdr < 0.05 ~ "significant",
      .default = "non-significant"
    )) |> 
  ggplot2::ggplot(ggplot2::aes(x = percentiles, y = unique_islet, fill = significance)) +
  ggplot2::geom_tile() +
  geom_tile(color = NA) +
  ggplot2::scale_fill_manual(
    values = c("significant" = "#A83708", "non-significant" = "white", "NA" = "lightgrey"),
    name = "FDR"
  ) +
  my_theme() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank()
  )
dev.off()
