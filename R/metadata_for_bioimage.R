# Description -------------------------------------------------------------
# Generating meta data to upload raw images to bioimage archive
mouse_info <- vroom::vroom(here::here("data/export/supplementary_files/Table 31.tsv")) |> 
  dplyr::select(sma_mouse_id = id, 
                age_weeks = age_w,
                cage,
                sex,
                line_strain)

df <- data.frame("Files" = list.files(path = here::here("data-raw/images/"), recursive = T, pattern = "*.nd2"))


# Add information ---------------------------------------------------------
df |> 
  tidyr::separate(Files, into = c("slide_name", "image_name"), remove = F, sep = "/") |> 
  tidyr::separate(slide_name, into = c("diet", "mouse_id", "slide_number"), remove = F, sep = "_") |> 
  dplyr::mutate(sma_mouse_id = paste0("SMA-", mouse_id),
                file_format = "nd2",
                FITC = "Insulin",
                TRITC = "STAT1",
                DAPI = "dapi") |> 
  dplyr::left_join(mouse_info) |> 
  dplyr::select(Files, file_format, slide_name, image_name, FITC, TRITC, DAPI, diet, mouse_id, age_weeks, sex, cage) |> 
  vroom::vroom_write(here::here("data/export/bioimage.tsv"))
