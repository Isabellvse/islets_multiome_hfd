# Description -------------------------------------------------------------
# Exporting excel results from various analysis

# Setup -------------------------------------------------------------------
base::source(here::here("R/set_up.R"))
set.seed(1000)

create_directories(c(here::here("data/export/")))
create_directories(c(here::here("data/export/supplementary_files")))

# Load --------------------------------------------------------------------
# Load old table
tables <- readxl::excel_sheets(here::here("data/export/supplementary_table.xlsx")) |> 
  purrr::set_names() |> 
  purrr::map(\(sheet){
    readxl::read_xlsx(here::here("data/export/supplementary_table.xlsx"), 
                      sheet = sheet,
                      .name_repair = snakecase::to_snake_case) 
  }) |> 
  purrr::compact()

# Load and rearrange the new data -----------------------------------------
table28 <- subcluster <- vroom::vroom(here::here("data/revisions/subcluster.csv"))
table29 <- vroom::vroom(here::here("data/revisions/cluster_composisiton.csv"))
table30 <- vroom::vroom(here::here("data/revisions/cluster3_markers_all.csv"))

table31 <- tables[["Table 29"]]
table32 <-  vroom::vroom(here::here("data/revisions/spatial_lfd_hfd/islet_data.csv"))
thresholds <- vroom::vroom(here::here("data/revisions/percentile/files/cell_data_status.csv")) |> 
  dplyr::select(unique_cell, tidyselect::starts_with("nucleus_p"), tidyselect::starts_with("peri_p"), tidyselect::ends_with("_or"))
table33 <- vroom::vroom(here::here("data/revisions/spatial_lfd_hfd/cell_data.csv")) |> 
  dplyr::left_join(thresholds)
table34 <- tables[["Table X"]]


# Edit overview -----------------------------------------------------------
overview <- tables[["Overview"]]

overview[29, 2] <- "Percentage of inflammatory response beta cells in subclusters"
overview[29, 3] <- "Percentage of inflammatory response beta cells in subclusters"
overview[30, 2] <- "DEGs between subcluster 3 and all other subcluters found in beta cells"
overview[30, 3] <- "Differentially expressed genes (DEGs) between various comparisons using pseudobulked counts, DESeq2 package and pairwise wald test. The comparisons include:  Cluster 3 vs all other subclusters"
overview[32, 1] <- "Table 32"
overview[32, 2] <- "Imaging analysis, islet level data"
overview[32, 3] <- "Imaging analysis, each row represents one segmented pancreatic islet. Intensity and area measurements for insulin, STAT1, and DAPI staining. Background refers to islet regions outside the nuclear and peri-nuclear compartments"
overview[33, 1] <- "Table 33"
overview[33, 2] <- "Imaging analysis, cell level data"
overview[33, 3] <- "Imaging analysis, each row represents one segmented cell inside an islet. Measurements the same as islets, as well as pixel distance from the islet edge. pX refers to percentile thresholds used for each compartment"
overview[34, 1] <- "Table 34"
overview[34, 2] <- "Overview of reagents and tools"
overview[34, 3] <- "Overview of reagents and tools used"


# Add to list -------------------------------------------------------------
tables[["Overview"]] <- overview
tables[["Table 28"]] <- table28
tables[["Table 29"]] <- table29
tables[["Table 30"]] <- table30
tables[["Table 31"]] <- table31
tables[["Table 32"]] <- table32
tables[["Table 33"]] <- table33
tables[["Table 34"]] <- table34
tables[["Table X"]] <- NULL


# Table 34 ----------------------------------------------------------------
table34 <- table34 %>%
  mutate(categories = if_else(is.na(source) & is.na(identifier), 
                              reagents_or_tools, 
                              NA_character_)) |> 
  fill(categories, .direction = "down") %>%
  filter(!(is.na(source) & is.na(identifier)))

## Add more info:
new_reagents <- tibble::tibble(
  reagents_or_tools = c(
    "Cover plate",
    "O.C.T",
    "Sequenza rack",
    "anti-STAT1",
    "anti-insulin",
    "DAPI",
    "Alexa Fluor 488-conjugated anti-guinea pig",
    "Alexa Fluor 546-conjugated anti-STAT1",
    "Vectashield antifade mounting medium",
    "Nikon Ti2 widefield microscope",
    "Teledyne Photometrics Kinetix sCMOS camera",
    "Andor Zyla 5.5 sCMOS camera"
  ),
  source = c(
    "Epredia",
    "Sakura Tissue-Tek",
    "Epredia",
    "Abcam",
    "Agilent",
    "Thermo Scientific",
    "Abcam",
    "Invitrogen",
    "Vector Laboratories",
    NA,
    NA,
    NA
  ),
  identifier = c(
    "72-110-017",
    "NC1862249",
    "73-310-017",
    "ab239360",
    "IR00251-2",
    "10116287",
    "ab150185",
    "A-11035",
    "H-1400-10",
    NA,
    NA,
    NA
  ),
  categories = c(
    "Other",
    "Chemicals and enzymes",
    "Other",
    "Dyes",
    "Dyes",
    "Dyes",
    "Dyes",
    "Dyes",
    "Chemicals and enzymes",
    "Other",
    "Other",
    "Other"
  )
)

table34 <- bind_rows(table34, new_reagents)
table34 <- table34 |> 
  dplyr::arrange(categories) |> 
  dplyr::relocate(categories)

# Add
tables[["Table 34"]] <- table34


# diets -------------------------------------------------------------------
table1 <- tables[["Table 1"]]
colnames(table1) <- c("item", "lfd_gm", "lfd_kcal", "hfd_gm", "hfd_kcal")
table1 <- table1[-1, ]  # drop the old header row that's now redundant
table1 <- table1 |> 
  mutate(categories = if_else(is.na(lfd_gm) & is.na(lfd_kcal) & is.na(hfd_gm) & is.na(hfd_kcal), 
                              item, 
                              NA_character_)) |> 
  fill(categories, .direction = "down") |> 
  dplyr::mutate(categories = if_else(is.na(categories), "Macronutrients", categories)) |> 
  dplyr::filter(!item == categories) |> 
  dplyr::relocate(categories) |> 
  dplyr::mutate(dplyr::across(c(-categories, -item), as.numeric))

tables[["Table 1"]] <- table1

# Save --------------------------------------------------------------------
# Save as tab seperated file 
tables |> 
  purrr::iwalk(\(table, name){
    vroom::vroom_write(table, paste0(here::here("data/export/supplementary_files/"), name, ".tsv"))
  })

# Compress
# tar -czvf supplementary_files.tar.gz ./supplementary_files

# save as xsls file 
openxlsx::write.xlsx(tables, here::here("data/export/supplementary_files.xlsx"))   


# Reorder tables again ----------------------------------------------------
tables <- readxl::excel_sheets(here::here("data/export/supplementary_files.xlsx")) |> 
  purrr::set_names() |> 
  purrr::map(\(sheet){
    readxl::read_xlsx(here::here("data/export/supplementary_files.xlsx"), 
                      sheet = sheet) 
  })

mapping <- c(
  "19" = 25, "20" = 26, "21" = 27, "22" = 28, "23" = 29, "24" = 30,
  "25" = 31, "26" = 32, "27" = 33, "28" = 19, "29" = 20, "30" = 21,
  "31" = 22, "32" = 23, "33" = 24, "34" = 34
)

overview <- tables[["Overview"]]
old_num <- as.numeric(gsub("Table ", "", overview$overview_of_tables))
new_num <- ifelse(as.character(old_num) %in% names(mapping),
                  mapping[as.character(old_num)],
                  old_num)
overview$overview_of_tables <- paste0("Table ", new_num)
overview <- overview[order(new_num), ]
tables[["Overview"]] <- overview

# --- rename and reorder the tables list itself ---
nm <- names(tables)
num <- as.numeric(gsub("Table ", "", nm))  # NA for "Overview"

new_num_list <- ifelse(!is.na(num) & as.character(num) %in% names(mapping),
                       mapping[as.character(num)],
                       num)

names(tables) <- ifelse(is.na(new_num_list), nm, paste0("Table ", new_num_list))

tables <- tables[order(ifelse(is.na(new_num_list), -Inf, new_num_list))]

names(tables)  # sanity check: "Overview" "Table 1" ... "Table 34"


# Save again --------------------------------------------------------------
# Save as tab seperated file 
tables |> 
  purrr::iwalk(\(table, name){
    vroom::vroom_write(table, paste0(here::here("data/export/supplementary_files/"), name, ".tsv"))
  })

# Compress
# tar -czvf supplementary_files.tar.gz ./supplementary_files

# save as xsls file 
openxlsx::write.xlsx(tables, here::here("data/export/supplementary_files.xlsx"))   
