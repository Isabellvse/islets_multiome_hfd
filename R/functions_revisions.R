# Functions ---------------------------------------------------------------
#' Create Directories if They Do Not Exist
#'
#' This function checks if each directory in the provided list exists. If any directory
#' does not exist, it creates the directory (along with any necessary parent directories).
#' It prints a message for each directory indicating whether it was created or already exists.
#'
#' @param output_dir A character vector or list of directory paths to be checked and created.
#'
#' @return NULL This function performs a side-effect (creating directories) and does not return anything.
#' 
#' @export
#'
#' @examples
#' # Define a list of directories
#' output_dirs <- list("data/quality_control/rna", "data/plots", "data/other_output")
#' 
#' # Call the function to create the directories
#' create_directories(output_dirs)
create_directories <- function(output_dir) {
  purrr::walk(output_dir, ~{
    if (!dir.exists(.x)) {
      dir.create(.x, recursive = TRUE)  # create the directory if it doesn't exist
      print(paste0(.x, " has been created!"))
    } else {
      print(paste0(.x, " already exists!"))
    }
  })
}

# Load insulin mask -------------------------------------------------------
#' Helper function to create owin from JSON mask data
#'
#' @param json_path path to json file
#' @param islet_id character string, islet id
#'
#' @returns retuns An object of the class owin (observation window in the two-dimensional plane)
create_window_from_json <- function(json_path, islet_id) {
  # Read the JSON file
  masks_list <- jsonlite::read_json(json_path)
  
  # Get the specific islet mask
  mask_data <- masks_list[[as.character(islet_id)]]
  
  if (is.null(mask_data)) {
    return(NULL)
  }
  
  # Create owin object
  spatstat.geom::owin(poly = list(
    x = unlist(mask_data$x),
    y = unlist(mask_data$y)
  ))
}

#' Convert owin into a dataframe - to be used for plotting
#'
#' @param window owin object
#' @param unique_id A character string, islet id 
#'
#' @returns A data frame 
window_to_df <- function(window, unique_id = NULL) {
  if (is.null(window)) return(NULL)
  
  # Extract polygon coordinates from window
  poly <- window$bdry[[1]]
  
  data.frame(
    unique_id = unique_id,
    x = poly$x,
    y = poly$y,
    stringsAsFactors = FALSE
  )
}


# Segregation test --------------------------------------------------------

#' Create a point pattern object on x and y positions
#'
#' @param islet_data a dataframe
#' @param islet_name a character string
#' @param marks marks of each point, e.g. "high" or "low"
#'
#' @returns A ppp object
make_ppp <- function(islet_data, islet_name, marks) {
  spatstat.geom::ppp(
    x = islet_data$centroid_x_px,
    y = islet_data$centroid_y_px,
    window = window_list[[islet_name]], 
    marks = marks
  )
}

#' Function for segregation test and permutation analysis
#'
#' @param islet_data A data frame
#' @param islet_name Name of islet
#' @param nsim Number of permutations to perform (default = 999)
#'
#' @returns A list 
segregation_test <- function(islet_data, islet_name, nsim = 999){
  if (length(unique(islet_data$stat1_status)) < 2) return(NULL)
  
  pp_all <- make_ppp(islet_data, islet_name, factor(islet_data$stat1_status))
  sigma <- spatstat.explore::bw.scott(pp_all)
  
  list(
    pp = pp_all,
    sigma = sigma,
    results = spatstat.explore::segregation.test(pp_all, sigma = sigma, nsim = nsim) |>
      broom::tidy() |>
      dplyr::mutate(nsim = nsim)
  )
}


#' Perform Estimate of Spatially-Varying Relative Risk on islets
#'
#' @param islet_data a dataframe
#' @param islet_name islet name
#' @param sigma what bandwith to use
#'
#' @returns A list
relrisk_test <- function(islet_data, islet_name, sigma){
  pp_all <- make_ppp(islet_data, islet_name, factor(islet_data$stat1_status, levels = c("low", "high")))
  
  list(
    pp = pp_all,
    vr = spatstat.explore::relrisk(pp_all, sigma =sigma, casecontrol = TRUE)
  )
}

plot_relrisk_ggplot <- function(relrisk_im, title = "") {
  as.data.frame(relrisk_im) |>
    ggplot2::ggplot(ggplot2::aes(x = x, y = y, fill = value)) +
    ggplot2::geom_raster() +
    ggplot2::scale_fill_viridis_c() +
    ggplot2::coord_equal() +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      legend.position = "none"
    ) +
    ggplot2::labs(title = title, fill = "Relative Risk")
}

plot_relrisk_ggplot_void <- function(relrisk_im, title = "") {
  as.data.frame(relrisk_im) |>
    ggplot2::ggplot(ggplot2::aes(x = x, y = y, fill = value)) +
    ggplot2::geom_raster() +
    ggplot2::scale_fill_viridis_c() +
    ggplot2::coord_equal() +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "none")
}

# relative risk ratio color limits between 0 and 1
plot_relrisk_ggplot_void_sig <- function(relrisk_im, is_significant = FALSE, title = "") {
  color_scale <- if (is_significant) {
    ggplot2::scale_fill_viridis_c(option = "plasma", limits = c(0, 1))
  } else {
    ggplot2::scale_fill_viridis_c(option = "viridis", limits = c(0, 1))
  }
  
  as.data.frame(relrisk_im) |>
    ggplot2::ggplot(ggplot2::aes(x = x, y = y, fill = value)) +
    ggplot2::geom_raster() +
    color_scale +
    ggplot2::coord_equal() +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "none")
}

plot_relrisk_ggplot_void_sig_legend <- function(relrisk_im, is_significant = FALSE, title = "") {
  color_scale <- if (is_significant) {
    ggplot2::scale_fill_viridis_c(option = "plasma", limits = c(0, 1))
  } else {
    ggplot2::scale_fill_viridis_c(option = "viridis", limits = c(0, 1))
  }
  
  as.data.frame(relrisk_im) |>
    ggplot2::ggplot(ggplot2::aes(x = x, y = y, fill = value)) +
    ggplot2::geom_raster() +
    color_scale +
    ggplot2::coord_equal() +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "right")
}

plot_pp_ggplot_void <- function(pp_obj, title = "") {
  points_df <- data.frame(
    x = pp_obj$x,
    y = pp_obj$y,
    marks = pp_obj$marks
  )
  
  # Get window as polygon (handles different window types)
  window_poly <- pp_obj$window |> as.data.frame()
  
  color_scheme <- c("low" = "blue", "high" = "red")
  
  ggplot2::ggplot() +
    ggrastr::geom_point_rast(data = points_df |> dplyr::filter(marks == "low"), 
                             ggplot2::aes(x = x, y = y, color = marks), scale = 0.2) +
    ggrastr::geom_point_rast(data = points_df |> dplyr::filter(marks == "high"), 
                             ggplot2::aes(x = x, y = y, color = marks), scale = 0.2) +
    ggplot2::geom_polygon(data = window_poly, ggplot2::aes(x = x, y = y), 
                          fill = NA, color = "black", linewidth = 0.1) +
    ggplot2::scale_color_manual(values = color_scheme) +
    ggplot2::coord_equal() +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "none")
}

plot_ordered_by_stat <- function(plot_list, seg_results, percentile_name, ncol = 30) {
  ordered_islets <- seg_results |>
    dplyr::filter(percentiles == percentile_name) |>
    dplyr::mutate(is_sig = fdr < 0.05) |> 
    dplyr::group_by(is_sig) |> 
    dplyr::arrange(fdr, desc(statistic)) |>
    dplyr::pull(unique_islet)
  
  patchwork::wrap_plots(plot_list[[percentile_name]][ordered_islets], ncol = ncol)
}


# spatial occupancy -------------------------------------------------------
#' Prepare data spatial analysis
#'
#' @param df A data frane
#'
#' @returns A data frame
prep_data <- function(df) {
  prep_df <- df |> 
    # Group by islet and structural zone
    dplyr::group_by(unique_islet, location, mouse_id, slide_name) |> 
    # Calculate the successes (high) and the trials (total cells in that zone)
    dplyr::summarise(
      n_cells = sum(stat1_status == "high"),
      total = n(),
      .groups = "drop_last"
    ) |> 
    dplyr::ungroup()
  
  return(prep_df)
}

#' Define model for spatial analysis
#'
#' @param prep_df A data frame (from prepare data)
#'
#' @returns A data frame
model_fun <- function(prep_df) {
  # Fits the proportion model while controlling for islet-to-islet variation
  model <- glmmTMB::glmmTMB(
    cbind(n_cells, total - n_cells) ~ location + (1|unique_islet),
    family = binomial,
    data = prep_df
  )
  
  # Tidy the fixed effects
  model_results <- broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE)
  return(model_results)
}

#' Perform permutation test on spatial analysis data
#'
#' @param df A data frame 
#' @param model_results A data frame
#' @param n_perm Number of permutations (default = 999)
#'
#' @returns A data frame
perm_model <- function(df, model_results, n_perm = 999){
  observed_intercept <- model_results |> dplyr::filter(term == "(Intercept)") |> dplyr::pull(estimate)
  observed_periphery <- model_results |> dplyr::filter(stringr::str_detect(term, "locationedge")) |> dplyr::pull(estimate)
  
  # Parallelize the permutation loop across your cluster
  perm_coefs <- furrr::future_map(1:n_perm, \(i) {
    
    # 1. Permute location labels within each islet at raw cell level
    df_perm <- df |>
      dplyr::group_by(unique_islet) |>
      dplyr::mutate(stat1_status = sample(stat1_status)) |>
      dplyr::ungroup()
    
    # 2. Aggregate the shuffled data structures
    prep_df <- prep_data(df_perm)
    
    # 3. Fit the null model
    model_perm <- glmmTMB::glmmTMB(
      cbind(n_cells, total - n_cells) ~ location + (1|unique_islet),
      family = binomial,
      data = prep_df
    )
    
    # 4. Extract both coefficients safely
    perm_tidy <- broom.mixed::tidy(model_perm, effects = "fixed")
    
    est_intercept <- perm_tidy |> dplyr::filter(term == "(Intercept)") |> dplyr::pull(estimate)
    est_periphery <- perm_tidy |> dplyr::filter(stringr::str_detect(term, "locationedge")) |> dplyr::pull(estimate)
    
    # Handle rare non-convergence edge cases by outputting NA_real_ if coefficients are missing
    tibble::tibble(
      intercept = if(length(est_intercept) > 0) est_intercept else NA_real_,
      periphery = if(length(est_periphery) > 0) est_periphery else NA_real_
    )
  }, .progress = TRUE, .options = furrr::furrr_options(seed = TRUE, scheduling = 2)) |>
    purrr::list_rbind() |> 
    tidyr::drop_na() # Cleanly discard any failed/unconverged permutation steps
  
  # two tailed p-values (Checks for extreme clustering in EITHER direction)
  perm_p_2tail_periphery <- (sum(abs(perm_coefs$periphery) >= abs(observed_periphery)) + 1) / (nrow(perm_coefs) + 1)
  perm_p_2tail_intercept <- (sum(abs(perm_coefs$intercept) >= abs(observed_intercept)) + 1) / (nrow(perm_coefs) + 1)
  
  # one tailed p-values (Checks if edge enrichment is greater than the null)
  if (observed_periphery >= 0) {
    perm_p_1tail_periphery <- (sum(perm_coefs$periphery >= observed_periphery) + 1) / (nrow(perm_coefs) + 1)
  } else {
    perm_p_1tail_periphery <- (sum(perm_coefs$periphery <= observed_periphery) + 1) / (nrow(perm_coefs) + 1)
  }
  
  if (observed_intercept >= 0) {
    perm_p_1tail_intercept <- (sum(perm_coefs$intercept >= observed_intercept) + 1) / (nrow(perm_coefs) + 1)
  } else {
    perm_p_1tail_intercept <- (sum(perm_coefs$intercept <= observed_intercept) + 1) / (nrow(perm_coefs) + 1)
  }
  
  # Add both sets of p-values to results 
  results <- model_results |>
    dplyr::mutate(
      p_value_2tailed = ifelse(term == "(Intercept)", perm_p_2tail_intercept, perm_p_2tail_periphery),
      p_value_1tailed = ifelse(term == "(Intercept)", perm_p_1tail_intercept, perm_p_1tail_periphery)
    )
  
  return(results)
}

#' Master function spatial analysis
#'
#' @param df A data frame
#'
#' @returns A data frame
perm_func <- function(df) {
  prep_df       <- prep_data(df)
  model_results <- model_fun(prep_df = prep_df)
  results       <- perm_model(df = df, model_results = model_results)
  return(results)
}

# Plotting ----------------------------------------------------------------
plot_islet <- function(islet_name, cell, status_col, percentile = 0.8, save_plot = TRUE) {
  
  # Filter data for this islet
  islet_data <- cell |> 
    dplyr::filter(unique_islet == islet_name)
  
  # Calculate quantiles for ALL cells (for p1)
  q_nuc_all <- quantile(cell$stat1_nuc_all, percentile, na.rm = TRUE)
  q_peri_all <- quantile(cell$stat1_peri_all, percentile, na.rm = TRUE)
  
  # Calculate quantiles for THIS islet only (for p2)
  q_nuc_islet <- quantile(islet_data$stat1_nuc_islet, percentile, na.rm = TRUE)
  q_peri_islet <- quantile(islet_data$stat1_peri_islet, percentile, na.rm = TRUE)
  
  # Plot 1: Scaled across all cells
  p1 <- islet_data |> 
    ggplot2::ggplot(ggplot2::aes(x = centroid_x_px, y = centroid_y_px)) +
    ggplot2::geom_point(ggplot2::aes(color = stat1_peri_all), size = 3) +
    ggplot2::geom_point(ggplot2::aes(fill = stat1_nuc_all), shape = 21, size = 1) +
    ggplot2::scale_fill_gradientn(
      colors = rev(roma(256)),
      oob = scales::squish,
      limits = c(0, q_nuc_all),
      name = "Nucleus"
    ) +
    ggplot2::scale_color_gradientn(
      colors = rev(roma(256)),
      oob = scales::squish,
      limits = c(0, q_peri_all),
      name = "Perinuclear\nspace"
    ) +
    ggplot2::scale_y_reverse() +
    ggplot2::labs(x = "Centroid x-axis (Pixels)",
                  y = "Centroid y-axis (Pixels)",
                  title = "STAT1 intensity",
                  subtitle = "Min-Max scaled across all cells") +
    my_theme()
  
  # Plot 2: Scaled within this islet
  p2 <- islet_data |> 
    ggplot2::ggplot(ggplot2::aes(x = centroid_x_px, y = centroid_y_px)) +
    ggplot2::geom_point(ggplot2::aes(color = stat1_peri_islet), size = 3) +
    ggplot2::geom_point(ggplot2::aes(fill = stat1_nuc_islet), shape = 21, size = 1) +
    ggplot2::scale_fill_gradientn(
      colors = rev(roma(256)),
      oob = scales::squish,
      limits = c(0, q_nuc_islet),  # Using islet-specific quantile
      name = "Nucleus"
    ) +
    ggplot2::scale_color_gradientn(
      colors = rev(roma(256)),
      oob = scales::squish,
      limits = c(0, q_peri_islet),  # Usingislet-specific quantile
      name = "Perinuclear\nspace"
    ) +
    ggplot2::scale_y_reverse() +
    ggplot2::labs(x = "Centroid x-axis (Pixels)",
                  y = "Centroid y-axis (Pixels)",
                  title = "STAT1 intensity",
                  subtitle = "Min-Max scaled within islet") +
    my_theme()
  
  # Plot 3: Status
  p3 <- islet_data |> 
    ggplot2::ggplot(ggplot2::aes(x = centroid_x_px, y = centroid_y_px)) +
    ggplot2::geom_point(ggplot2::aes(fill = stat1_status), size = 3, shape = 21) +
    ggplot2::scale_fill_manual(values = c("low" = "blue", "high" = "red")) +
    ggplot2::scale_y_reverse() +
    my_theme()
  
  # Combine plots
  combined <- p1 + p2 + p3
  # Optionally save
  if (save_plot) {
    ggplot2::ggsave(
      here::here(glue::glue("multiome_revisions/data/percentile/plots/{islet_name}.pdf")),
      plot = combined,
      width = 12,
      height = 4
    )
  }
  
  combined_black <- combined &
    theme(
      plot.background = element_rect(fill = "black", color = NA),
      panel.background = element_rect(fill = "black", color = NA),
      text = element_text(color = "white"),
      axis.text = element_text(color = "white"),
      axis.title = element_text(color = "white"),
      legend.background = element_rect(fill = "black", color = NA),
      legend.text = element_text(color = "white"),
      legend.title = element_text(color = "white"),
      panel.grid = element_line(color = "gray30"),
      axis.line = element_line(color = "white")
    )
  
  return(combined_black)
}



roma <-khroma::color("roma")
diet_lvl <- c("lfd", "hfd")
my_theme <- function() {
  ggplot2::theme_classic(base_size = 7) +
    ggplot2::theme(strip.background = element_blank())
}

condition_color <- c(
  "lfd" = "#004B7A",
  "hfd" = "#fcb21c"
)

stat1_status_color <- c("high" = "red",
                        "low" = "blue")
