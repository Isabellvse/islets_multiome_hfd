# Description -------------------------------------------------------------
# Combine all cell and islet data, remove islets with less than 10 cells
# Setup -------------------------------------------------------------------
source(here::here("islets_multiome_hfd/R/set_up_revisions.R"))

# Load --------------------------------------------------------------------
# Cell data ----
cell <- list.files(
  here::here("islets_multiome_hfd/data/revisions/spatial_lfd_hfd/files"),
  pattern = "_cell\\.csv$",
  recursive = TRUE,
  full.names = TRUE
) |>
  map(vroom::vroom, show_col_types = FALSE, num_threads = 30) |>
  list_rbind() |>
  mutate(
    diet = str_extract(slide_name, "^[^_]+"),
    mouse_id = str_extract(slide_name, "\\d{4}"),
    slide_number = str_extract(slide_name, "[^_]$"),
    unique_islet = paste(slide_name, image_name, islet_id, sep = "_"),
    unique_cell = paste(slide_name, image_name, islet_id, cell_id, sep = "_")
  )

# Islet data ----
islet <- list.files(
  here::here("islets_multiome_hfd/data/revisions/spatial_lfd_hfd/files"),
  pattern = "_islet\\.csv$",
  recursive = TRUE,
  full.names = TRUE
) |>
  map(vroom::vroom, show_col_types = FALSE, num_threads = 30) |>
  list_rbind() |>
  mutate(
    diet = str_extract(slide_name, "^[^_]+"),
    mouse_id = str_extract(slide_name, "\\d{4}"),
    slide_number = str_extract(slide_name, "[^_]$"),
    unique_islet = paste(slide_name, image_name, islet_id, sep = "_")
  )

# Insulin mask ----
insulin <- data.frame("path" = list.files(
  here::here("islets_multiome_hfd/data/revisions/spatial_lfd_hfd/insulin_masks"),
  pattern = "_insulin_masks\\.json$",
  recursive = TRUE,
  full.names = TRUE)) |> 
  dplyr::mutate(
    # Get folder name (handles "hfd_0615_2" and "hfd_0615_2 v2")
    slide_name = basename(dirname(path)),
    # Get image name from filename (Seq0010, Seq0007, etc.)
    image_name = stringr::str_extract(basename(path), "^[^_]+")
  )

# Load all insulin masks
insulin <- data.frame("path" = list.files(
  here::here("islets_multiome_hfd/data/revisions/spatial_lfd_hfd/insulin_masks"),
  pattern = "_insulin_masks\\.json$",
  recursive = TRUE,
  full.names = TRUE)) |>
  dplyr::mutate(
    slide_name = basename(dirname(path)),
    image_name = stringr::str_extract(basename(path), "^[^_]+")
  ) |>
  dplyr::mutate(
    # Read JSON for each path
    masks_list = purrr::map(path, jsonlite::read_json),
    # Extract islet IDs from the JSON (they're the keys)
    islet_ids = purrr::map(masks_list, names)
  ) |>
  # Unnest to get one row per islet
  tidyr::unnest(islet_ids) |>
  dplyr::mutate(
    islet_id = as.numeric(islet_ids),
    # Create unique identifier
    unique_islet = paste(slide_name, image_name, islet_id, sep = "_"),
    # Create window for each islet
    window = purrr::map2(path, islet_id, create_window_from_json),
  ) |>
  dplyr::select(unique_islet, slide_name, image_name, islet_id, 
                window, path)

# Filter islets -----------------------------------------------------------
# Filter islets that are less than 10 cells
islet_filter <- cell |> 
  dplyr::group_by(unique_islet) |> 
  tally() |> 
  dplyr::filter(n > 10) |> 
  dplyr::pull(unique_islet)

# Filter
cell <- cell |> 
  dplyr::filter(unique_islet %in% islet_filter)

islet <- islet |> 
  dplyr::filter(unique_islet %in% islet_filter)

insulin <- insulin |> 
  dplyr::filter(unique_islet %in% islet_filter)


# Scale cell data ---------------------------------------------------------
cell <- cell |>
  # Scale ACROSS all cells
  dplyr::mutate(
    stat1_peri_scale_all = (stat1_peri_mean - min(stat1_peri_mean, na.rm = TRUE)) /
      (
        max(stat1_peri_mean, na.rm = TRUE) - min(stat1_peri_mean, na.rm = TRUE)
      ),
    stat1_nuc_scale_all = (stat1_nucleus_mean - min(stat1_nucleus_mean, na.rm = TRUE)) /
      (
        max(stat1_nucleus_mean, na.rm = TRUE) - min(stat1_nucleus_mean, na.rm = TRUE)
      )
  ) |>
  # Scale WITHIN each islet
  dplyr::group_by(unique_islet) |>
  dplyr::mutate(
    stat1_peri_scale_islet = (stat1_peri_mean - min(stat1_peri_mean, na.rm = TRUE)) /
      (
        max(stat1_peri_mean, na.rm = TRUE) - min(stat1_peri_mean, na.rm = TRUE)
      ),
    stat1_nuc_scale_islet = (stat1_nucleus_mean - min(stat1_nucleus_mean, na.rm = TRUE)) /
      (
        max(stat1_nucleus_mean, na.rm = TRUE) - min(stat1_nucleus_mean, na.rm = TRUE)
      ),
    distance_edge_px_scale_islet = (
      distance_from_edge_px - min(distance_from_edge_px, na.rm = TRUE)
    ) /
      (
        max(distance_from_edge_px, na.rm = TRUE) - min(distance_from_edge_px, na.rm = TRUE)
      )
  ) |>
  dplyr::ungroup()


# Create windows ----------------------------------------------------------
# Create named list of windows
window_list <- stats::setNames(insulin$window, insulin$unique_islet)

# Convert wndows to dataframe
windows_df <- purrr::map2_df(
  window_list, 
  names(window_list),
  ~window_to_df(.x, .y)
)

# Save --------------------------------------------------------------------
vroom::vroom_write(cell, here::here("islets_multiome_hfd/data/revisions/spatial_lfd_hfd/cell_data.csv"))
vroom::vroom_write(islet, here::here("islets_multiome_hfd/data/revisions/spatial_lfd_hfd/islet_data.csv"))
vroom::vroom_write(windows_df, here::here("islets_multiome_hfd/data/revisions/spatial_lfd_hfd/window_data.csv"))
saveRDS(window_list, here::here("islets_multiome_hfd/data/revisions/spatial_lfd_hfd/window_list.rds"))
