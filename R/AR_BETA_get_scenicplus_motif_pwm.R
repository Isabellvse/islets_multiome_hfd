# Description -------------------------------------------------------------
# Genereate motif pwm list from pwms used in scenicplus

# Setup -------------------------------------------------------------------
base::source(here::here("R/set_up.R"))
set.seed(1000)

# Download ----------------------------------------------------------------
# in terminal
#wget https://resources.aertslab.org/cistarget/motif_collections/v10nr_clust_public/v10nr_clust_public.zip
#unzip v10nr_clust_public.zip

# Import files ------------------------------------------------------------
## motif annotation for mouse ----
motif_anno <- RcisTarget::importAnnotations(here::here("data/pwms/v10nr_clust_public/snapshots/motifs-v10-nr.mgi-m0.00001-o0.0.tbl"))

## Position weight matrix in homer styles ----
# Get paths
inputfiles <- base::list.files(path = here::here("data/pwms/v10nr_clust_public/singletons"),
                         pattern = ".cb", full.names = TRUE)

# keep only files that are within the mouse annotation
file_names <- inputfiles %>%
  stringr::str_remove("/work/islets_multiome_hfd/data/pwms/v10nr_clust_public/singletons/") %>%
  stringr::str_remove(".cb")

keep <- intersect(file_names, motif_anno$motif)
keep_path <- paste0("/work/islets_multiome/data/pwms/v10nr_clust_public/singletons/", keep, ".cb")


# Import files that are not meta clusters ---------------------------------
# path to singlet files
singlet <- keep_path[-base::grep("metacluster", keep_path)]

# import
singlet_pwm <- purrr::map(singlet, read_homr2)

# clean up names
singlet_pwm <- singlet_pwm %>%
  purrr::set_names(~ .x %>%
                     stringr::str_remove("/work/islets_multiome_hfd/data/pwms/v10nr_clust_public/singletons/") %>%
                     stringr::str_remove(".cb"))
# Import metaclusters -----------------------------------------------------

# path to meta clusters
meta <- keep_path[base::grep("metacluster", keep_path)]

# import
meta_pwm <- purrr::map(singlet, read_homr2)

# clean up names
meta_pwm <- meta_pwm %>%
  purrr::set_names(~ .x %>%
                     stringr::str_remove("/work/islets_multiome_hfd/data/pwms/v10nr_clust_public/singletons/") %>%
                     stringr::str_remove(".cb"))

# meta clusters that contains more than 1 pwm, will be renamed to contain name of each motif
meta_pwm <- purrr::map(meta_pwm, function(x) {
  
  # If x is a list, extract names of the motifs
  if (is(x, "list")) {
    
    # Extract the name of each motif from the list elements (sometimes clusters of motif have the same name)
    names_clean <- purrr::map_chr(x, ~ .x@name) %>%
      stringr::str_remove(">") %>%
      make.unique()
    
    # Assign cleaned and unique names back to the elements of the list
    base::names(x) <- names_clean
  }
  
  # Return the modified object (whether list or not)
  return(x)
})

# unlist
meta_pwm <- base::unlist(meta_pwm)


# Combine pwms ------------------------------------------------------------
pwms <- c(singlet_pwm, meta_pwm)

# convert to PWMatrix
pwmlist <- universalmotif::convert_motifs(pwms, "TFBSTools-PWMatrix")

# convert to PWMatrixList
pwmlist <- do.call(TFBSTools::PWMatrixList, pwmlist)

# save list
base::saveRDS(pwmlist, here::here("data/scenicplus/scenic_PWMLIST_mouse.rds"))
