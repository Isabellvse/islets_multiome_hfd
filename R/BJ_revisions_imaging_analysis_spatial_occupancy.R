# Description -------------------------------------------------------------
# Do STAT1 high cells preferentially appear in the periphery of the islet?
# But defined as logodss

# Setup -------------------------------------------------------------------
source(here::here("islets_multiome_hfd/R/set_up_revisions.R"))
create_directories(here::here("data/resivions/percentile/plots"))
create_directories(here::here("data/resivions/percentile/files"))
set.seed(1000)

# Load --------------------------------------------------------------------
cell <- vroom::vroom(here::here("data/resivions/percentile/files/cell_data_status.csv")) |> 
  dplyr::mutate(diet = factor(diet, diet_lvl))
islet <- vroom::vroom(here::here("data/resivions/spatial_lfd_hfd/islet_data.csv")) |> 
  dplyr::mutate(diet = factor(diet, diet_lvl))
window_list <- readRDS(here::here("data/resivions/spatial_lfd_hfd/window_list.rds"))

# Prepare data ------------------------------------------------------------
per_cell_or <- cell |> 
  dplyr::filter(diet == "hfd") |> 
  dplyr::select(unique_islet, unique_cell, mouse_id, slide_name,
                centroid_x_px, centroid_y_px, distance_edge_px_scale_islet,distance_from_edge_px, 
                diet, mouse_id, ends_with("or") & contains("stat1_status"), stat1_nucleus_mean, stat1_peri_mean) |> 
  tidyr::pivot_longer(ends_with("_or"), names_to = "percentiles", values_to = "stat1_status") |> 
  base::split(~percentiles) |> 
  purrr::map(\(df){
    df |> 
      dplyr::group_by(unique_islet) |>
      dplyr::filter(n_distinct(stat1_status) == 2) |> 
      dplyr::mutate(
        distance_q30 = quantile(distance_from_edge_px, 0.30, na.rm = TRUE),
        location = dplyr::case_when(
          distance_from_edge_px <= distance_q30 ~ "edge",
          .default = "core"
        )
      ) |>
      dplyr::ungroup()
  })


# Plot --------------------------------------------------------------------
figlocation <- per_cell_or[["stat1_status_p95_or"]] |>
  ggplot(aes(x = centroid_x_px, y = centroid_y_px, color = location)) +
  scale_color_manual(values = c("edge" = "#B2182B", "core" = "#2166AC")) +
  ggrastr::geom_point_rast(size = 0.5) +
  facet_wrap(~unique_islet, scales = "free", ncol = 40, labeller = labeller(
    unique_islet = function(x)
      ggplot2::label_wrap_gen(width = 12)(
        stringr::str_to_sentence(gsub("_", " ", x))
      ))) +
  labs(
    title = "Location of cells",
    subtitle = "Each point is a cell",
    x = "X coordinate (pixels)",
    y = "Y coordinate (pixels)"
  ) +
  theme_void() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    aspect.ratio = 1,
    strip.text = element_text(size = 5),
    legend.position = "right"
  )

ggsave(plot = figlocation, filename = here::here("data/resivions/percentile/plots/spatial_localization_core_edge_p95.png"), 
       dpi = 300,
       width = 20, height = 10)

# Run per test ------------------------------------------------------------
plan(multisession, workers = 63)  
perm_test <- names(per_cell_or) |>
  purrr::set_names() |> 
  purrr::map_df(\(percentile) {
    per_cell_or[[percentile]] |> 
      perm_func()
  }, .id = "percentile")
vroom::vroom_write(perm_test, here::here("data/resivions/percentile/files/core_periphety_permutation_test.csv"))
plan(sequential)

pdf(here::here("data/resivions/percentile/plots/core_edge_test.pdf"), 
    width = 1.2, height = 1.5)
perm_test |> 
  dplyr::filter(term == "locationedge") |> 
  dplyr::mutate(
    OR = exp(estimate),
    OR_low = exp(conf.low),
    OR_high = exp(conf.high),
    percentiles = stringr::str_extract(percentile, "p[:digit:]{2}"),
    star = dplyr::case_when(p_value_2tailed < 0.05 ~ "*", .default = ""),
    y_pos = OR
  ) |> 
  ggplot2::ggplot(ggplot2::aes(x = percentiles, y = OR)) +
  ggplot2::geom_bar(stat = "identity", fill = "#B2182B") +
  # CRITICAL FIX: Use geom_errorbar for vertical intervals on a standard bar chart
  ggplot2::geom_errorbar(ggplot2::aes(ymin = OR_low, ymax = OR_high), 
                         width = 0.2, color = "#2B2B2B", linewidth = 0.5) +
  ggplot2::geom_text(ggplot2::aes(y = y_pos, label = star),
                     size = 5, color = "black", vjust = 0) +
  ggplot2::labs(y = "Odds ratio (Edge vs Core)",
                x = "Percentile threshold") +
  my_theme()
dev.off()
