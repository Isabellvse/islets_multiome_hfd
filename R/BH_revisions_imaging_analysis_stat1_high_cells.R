# Description -------------------------------------------------------------
# Heterogenity of cells that have high STAT1 expression

# Setup -------------------------------------------------------------------
source(here::here("islets_multiome_hfd/R/set_up_revisions.R"))
create_directories(here::here("islets_multiome_hfd/data/revisions/percentile/plots"))
create_directories(here::here("islets_multiome_hfd/data/revisions/percentile/files"))

# Load --------------------------------------------------------------------
cell <- vroom::vroom(here::here("islets_multiome_hfd/data/revisions/spatial_lfd_hfd/cell_data/revisions.csv")) |> 
  dplyr::mutate(diet = factor(diet, diet_lvl))
islet <- vroom::vroom(here::here("islets_multiome_hfd/data/revisions/spatial_lfd_hfd/islet_data/revisions.csv")) |> 
  dplyr::mutate(diet = factor(diet, diet_lvl))

logic_types <- c("and", "only_peri", "only_nuc", "or")
percentiles <- c("p70", "p75", "p80", "p85", "p90", "p95")

# Define STAT1 high cells -------------------------------------------------
percentile_data/revisions <- cell |> 
  dplyr::select(stat1_nucleus_mean, stat1_peri_mean) |> 
  tidyr::pivot_longer(everything(), names_to = "name") |>
  dplyr::group_by(name) |>
  dplyr::summarise(dplyr::across(value, list(
    p70 = ~quantile(., 0.7, na.rm = TRUE),
    p75 = ~quantile(., 0.75, na.rm = TRUE),
    p80 = ~quantile(., 0.8, na.rm = TRUE),
    p85 = ~quantile(., 0.85, na.rm = TRUE),
    p90 = ~quantile(., 0.90, na.rm = TRUE),
    p95 = ~quantile(., 0.95, na.rm = TRUE)
  )), .groups = "drop") |>
  tidyr::pivot_longer(cols = -name, names_to = "percentile", values_to = "value") |>
  dplyr::mutate(percentile = stringr::str_remove(percentile, "value_"))

# Plot distribution with percentile lines
pdf(here::here("islets_multiome_hfd/data/revisions/percentile/plots/stat1_intensity_distribution_cell.pdf"), 
    width = 4, height = 2.5)
cell |> 
  dplyr::select(stat1_nucleus_mean, stat1_peri_mean) |>
  tidyr::pivot_longer(everything(), names_to = "name") |> 
  ggplot2::ggplot(ggplot2::aes(x = value)) +
  ggplot2::geom_density(alpha = 0.5, fill = "lightblue") +
  ggplot2::geom_vline(data/revisions = percentile_data/revisions,
                      ggplot2::aes(xintercept = value, color = percentile),
                      linetype = 2, size = 0.5) +
  ggplot2::facet_wrap(~name, scales = "free_x") +
  ggplot2::labs(x = "Mean STAT1 intensity (A.U.)",
                title = "Distribution of STAT1 intensity",
                subtitle = "Per cell",
                color = "Percentile") +
  my_theme()
dev.off()

## Define STAT1 status based on percentiles ----
# Create threshold columns and add status classifications
thresholds <- percentile_data/revisions |> 
  dplyr::mutate(name = stringr::str_remove(name, "stat1_"),
                name = stringr::str_remove(name, "_mean")) |> 
  tidyr::pivot_wider(names_from = c(name, percentile), values_from = value)

cell <- cell |> dplyr::bind_cols(thresholds)


# Define STAT high cells ---------------------------------------------------

# Define logic functions
logic_conditions <- list(
  "and" = function(nuc, peri, nuc_threshold, peri_threshold) {
    (nuc > nuc_threshold) & (peri > peri_threshold)
  },
  "only_peri" = function(nuc, peri, nuc_threshold, peri_threshold) {
    (peri > peri_threshold) & (nuc <= nuc_threshold)
  },
  "only_nuc" = function(nuc, peri, nuc_threshold, peri_threshold) {
    (nuc > nuc_threshold) & (peri <= peri_threshold)
  },
  "or" = function(nuc, peri, nuc_threshold, peri_threshold) {
    (nuc > nuc_threshold) | (peri > peri_threshold)
  }
)

# Apply mutations for each percentile and logic type
for (perc in percentiles) {
  nuc_col <- paste0("nucleus_", perc)
  peri_col <- paste0("peri_", perc)
  
  for (logic_type in names(logic_conditions)) {
    status_col <- paste0("stat1_status_", perc, "_", logic_type)
    condition_fn <- logic_conditions[[logic_type]]
    
    cell <- cell |>
      dplyr::mutate(
        !!status_col := dplyr::if_else(
          condition_fn(
            stat1_nucleus_mean,
            stat1_peri_mean,
            .data/revisions[[nuc_col]],
            .data/revisions[[peri_col]]
          ),
          "high", "low"
        )
      )
  }
}

# STAT1 status per islet --------------------------------------------------
per_islet <- cell |>
  dplyr::select(unique_islet, diet, mouse_id, slide_name, unique_cell, 
                matches("^stat1_status_")) |>
  tidyr::pivot_longer(
    cols = matches("^stat1_status_"),
    names_to = "percentile",
    values_to = "status"
  ) |>
  dplyr::group_by(percentile, status, diet, mouse_id, slide_name, unique_islet) |>
  dplyr::summarise(
    cells = n_distinct(unique_cell),
    .groups = "drop"
  ) |>
  dplyr::group_by(diet, percentile, mouse_id, slide_name, unique_islet) |>
  dplyr::mutate(
    total_cells = sum(cells),
    percentage = round((cells / total_cells) * 100, 2)
  ) |>
  dplyr::ungroup() |> 
  dplyr::mutate(percentile = stringr::str_remove(percentile, "stat1_status_")) |> 
  tidyr::pivot_wider(names_from = c(status, percentile), values_from = c(percentage, cells), values_fill = 0)




# Create plots for all combinations
plot_islet <- purrr::map(percentiles, ~ {
  perc <- .x
  
  purrr::map(logic_types, ~ {
    logic <- .x
    y_var_name <- paste0("percentage_high_", perc, "_", logic)
    
    analyze_and_plot_violin(
      list(
        y_var = y_var_name,
        subtitle = "Per islet",
        title = glue::glue("Distributions of STAT1 high cells\n({perc}) - {toupper(logic)}"),
        y_axis = "STAT1 high cells (%)"
      ),
      data/revisions = per_islet
    )
  })
})


long <- map(plot_islet, ~wrap_plots(.x, ncol = 1))

pdf(here::here("islets_multiome_hfd/data/revisions/percentile/plots/pct_stat1_high_per_islet.pdf"), 
    width = 15, height = 15)
wrap_plots(long, nrow = 1)
dev.off()


# STAT1 status per mouse --------------------------------------------------
per_mouse <- cell |>
  dplyr::select(diet, mouse_id, unique_cell, matches("^stat1_status_")) |>
  tidyr::pivot_longer(
    cols = matches("^stat1_status_"),
    names_to = "percentile",
    values_to = "status"
  ) |>
  dplyr::group_by(percentile, status, diet, mouse_id) |>
  dplyr::summarise(
    cells = n_distinct(unique_cell),
    .groups = "drop"
  ) |>
  dplyr::group_by(diet, percentile, mouse_id) |>
  dplyr::mutate(
    total_cells = sum(cells),
    percentage = round((cells / total_cells) * 100, 2)
  ) |>
  dplyr::ungroup() |> 
  dplyr::mutate(percentile = stringr::str_remove(percentile, "stat1_status_")) |> 
  tidyr::pivot_wider(names_from = c(status, percentile), values_from = c(percentage, cells), values_fill = 0)


plot_mouse<- purrr::map(percentiles, ~ {
  perc <- .x
  
  purrr::map(logic_types, ~ {
    logic <- .x
    y_var_name <- paste0("percentage_high_", perc, "_", logic)
    
    analyze_and_plot_wilcox_bar(
      list(
        y_var = y_var_name,
        subtitle = "Per mouse",
        title = glue::glue("Distributions of STAT1 high cells\n({perc}) - {toupper(logic)}"),
        y_axis = "STAT1 high cells (%)"
      ),
      data/revisions = per_mouse
    )
  })
})

long <- map(plot_mouse, ~wrap_plots(.x, ncol = 1))
pdf(here::here("islets_multiome_hfd/data/revisions/percentile/plots/pct_stat1_high_per_mouse.pdf"), 
    width = 15, height = 15)
wrap_plots(long, nrow = 1)
dev.off()


# WE will use OR

# FC plot over thresholds -------------------------------------------------
per_mouse_or_wilcox <- per_mouse |> 
  dplyr::select(diet, mouse_id, total_cells, ends_with("or") & contains("percentage_high")) |> 
  dplyr::mutate(diet = factor(diet, levels = c("hfd", "lfd"))) |> 
  tidyr::pivot_longer(ends_with("_or")) |> 
  (\(df) base::split(df, factor(df$name)))() |> 
  purrr::map(~rstatix::wilcox_test(.x, value ~ diet, detailed = T)) |> 
  purrr::list_rbind(names_to = "percentiles")

per_mouse_or_eff <- per_mouse |> 
  dplyr::select(diet, mouse_id, total_cells, ends_with("or") & contains("percentage_high")) |> 
  dplyr::mutate(diet = factor(diet, levels = c("hfd", "lfd"))) |> 
  tidyr::pivot_longer(ends_with("_or")) |> 
  (\(df) base::split(df, factor(df$name)))() |> 
  purrr::map(~rstatix::wilcox_effsize(.x, value ~ diet, detailed = T)) |> 
  purrr::list_rbind(names_to = "percentiles")

per_mouse_or_mean <- per_mouse |> 
  dplyr::select(diet, mouse_id, total_cells, ends_with("or") & contains("percentage_high")) |> 
  dplyr::mutate(diet = factor(diet, levels = c("hfd", "lfd"))) |> 
  tidyr::pivot_longer(ends_with("_or")) |> 
  dplyr::group_by(diet, percentiles = name) |> 
  dplyr::summarise(mean = mean(value),
                   sd = sd(value), .groups = "drop") |> 
  tidyr::pivot_wider(names_from = diet, values_from = c(mean, sd))

per_mouse_or <- per_mouse_or_wilcox |> 
  full_join(per_mouse_or_eff) |> 
  left_join(per_mouse_or_mean) |> 
  dplyr::mutate(star = dplyr::case_when(p < 0.05 ~ "*",
                                        .default = ""))


# Plot mean values --------------------------------------------------------
pdf(here::here("islets_multiome_hfd/data/revisions/percentile/plots/high_or_mean_point.pdf"), 
    width = 2.5, height = 2.5)
bind_rows(
  per_mouse_or |> 
    select(percentiles, mean = mean_hfd, sd = sd_hfd, star) |>
    mutate(group = "hfd"),
  per_mouse_or |> 
    select(percentiles, mean = mean_lfd, sd = sd_lfd) |>
    mutate(group = "lfd")
) |>
  dplyr::mutate(percentiles = stringr::str_extract(percentiles, "p[:digit:]{2}"), 
                y_pos = mean + 10) |> 
  
  ggplot2::ggplot(ggplot2::aes(x = percentiles, y = mean, color = group, group = group)) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = mean - sd, ymax = mean + sd),
    width = 0.2,
    linewidth = 0.5
  ) +
  ggplot2::geom_line() +
  ggplot2::geom_point(size = 2) +
  ggplot2::expand_limits(y = 0) +
  ggplot2::geom_text(ggplot2::aes(y = y_pos, label = star), size = 5, color = "black") +
  ggplot2::scale_color_manual(values = condition_color, name = "Diet") +
  my_theme() +
  ggplot2::labs(
    x = "Percentile Threshold",
    y = "Mean (%)",
    title = "STAT high cells percentages\nacross percentile thresholds",
    subtitle = "Per Mouse"
  ) 
dev.off()


# islet size and pct high cells -------------------------------------------
per_cell_or <- cell |> 
  dplyr::select(unique_islet, unique_cell, mouse_id,
                centroid_x_px, centroid_y_px, distance_from_edge_px, 
                diet, mouse_id, ends_with("or") & contains("stat1_status"), stat1_nucleus_mean, stat1_peri_mean) |> 
  tidyr::pivot_longer(ends_with("_or"), names_to = "percentiles", values_to = "stat1_status") |> 
  split(~percentiles)

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

p1 <- per_islet |> 
  dplyr::filter(diet == "hfd") |> 
  dplyr::left_join(y = islet |> dplyr::select(unique_islet, diet, slide_name, insulin_area_um2)) |> 
  ggplot2::ggplot(ggplot2::aes(x = pct_high, y = insulin_area_um2)) +
  ggrastr::geom_point_rast(scale = 0.5) +
  ggplot2::geom_smooth(method = "lm") +
  ggplot2::labs(title = "% STAT1 positive cells and Area",
                subtitle = "HFD islets",
                y = "Insulin area (um2)",
                x = "STAT1 high cells (%)") +
  ggpubr::stat_cor(aes(label = paste(after_stat(rr.label), p.label, sep = "~ `,`~")), color = "brown", label.y.npc = 1) +
  ggpubr::stat_regline_equation(color = "brown", label.y.npc = 0.8) +
  ggplot2::facet_wrap(~percentiles, nrow = 1) +
  my_theme()

p2 <- per_islet |> 
  dplyr::filter(diet == "hfd") |> 
  dplyr::left_join(y = islet |> dplyr::select(unique_islet, diet, slide_name, insulin_area_um2)) |> 
  ggplot2::ggplot(ggplot2::aes(x = pct_high, y = insulin_area_um2)) +
  ggrastr::geom_point_rast(scale = 0.5) + 
  ggplot2::geom_smooth(method = "lm") +
  ggplot2::labs(title = "% STAT1 positive cells and Area",
                subtitle = "HFD islets",
                y = "Total number of cells",
                x = "STAT1 high cells (%)") +
  ggpubr::stat_cor(aes(label = paste(after_stat(rr.label), p.label, sep = "~ `,`~")), color = "brown", label.y.npc = 1) +
  ggpubr::stat_regline_equation(color = "brown", label.y.npc = 0.8) +
  ggplot2::facet_wrap(~percentiles, nrow = 1) +
  my_theme()
pdf(here::here("islets_multiome_hfd/data/revisions/percentile/plots/area_hfd_high_stat1_correlation.pdf"), 
    width = 5, height = 3)
p1/p2
dev.off()
# Save --------------------------------------------------------------------
vroom::vroom_write(per_mouse_or, here::here("islets_multiome_hfd/data/revisions/percentile/files/wilcox_test_per_mouse_or.csv"))
vroom::vroom_write(cell, here::here("islets_multiome_hfd/data/revisions/percentile/files/cell_data/revisions_status.csv"))

cell <- vroom::vroom(here::here("islets_multiome_hfd/data/revisions/percentile/files/cell_data/revisions_status.csv"))
cell %>% dplyr::group_by(unique_islet, stat1_status_p95_or) %>% 
  tally() %>% 
  tidyr::pivot_wider(names_from = stat1_status_p95_or, values_from = n) %>% 
  dplyr::arrange(desc(high))
