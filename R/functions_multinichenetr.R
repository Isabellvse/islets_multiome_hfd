make_circos_one_group_2 <- function (prioritized_tbl_oi, colors_sender, colors_receiver) 
{
  requireNamespace("dplyr")
  requireNamespace("ggplot2")
  requireNamespace("circlize")
  prioritized_tbl_oi = prioritized_tbl_oi %>% dplyr::ungroup()
  grid_col_tbl_ligand = tibble::tibble(sender = colors_sender %>% 
                                         names(), color_ligand_type = colors_sender)
  grid_col_tbl_receptor = tibble::tibble(receiver = colors_receiver %>% 
                                           names(), color_receptor_type = colors_receiver)
  groups_oi = prioritized_tbl_oi$group %>% unique()
  all_plots = groups_oi %>% lapply(function(group_oi) {
    title = group_oi
    circos_links_oi = prioritized_tbl_oi %>% dplyr::filter(group == 
                                                             group_oi)
    circos_links = circos_links_oi %>% dplyr::rename(weight = prioritization_score)
    df = circos_links
    ligand.uni = unique(df$ligand)
    for (i in 1:length(ligand.uni)) {
      df.i = df[df$ligand == ligand.uni[i], ]
      sender.uni = unique(df.i$sender)
      for (j in 1:length(sender.uni)) {
        df.i.j = df.i[df.i$sender == sender.uni[j], ]
        df.i.j$ligand = paste0(df.i.j$ligand, paste(rep(" ", 
                                                        j - 1), collapse = ""))
        df$ligand[df$id %in% df.i.j$id] = df.i.j$ligand
      }
    }
    receptor.uni = unique(df$receptor)
    for (i in 1:length(receptor.uni)) {
      df.i = df[df$receptor == receptor.uni[i], ]
      receiver.uni = unique(df.i$receiver)
      for (j in 1:length(receiver.uni)) {
        df.i.j = df.i[df.i$receiver == receiver.uni[j], 
        ]
        df.i.j$receptor = paste0(df.i.j$receptor, paste(rep(" ", 
                                                            j - 1), collapse = ""))
        df$receptor[df$id %in% df.i.j$id] = df.i.j$receptor
      }
    }
    intersecting_ligands_receptors = generics::intersect(unique(df$ligand), 
                                                         unique(df$receptor))
    while (length(intersecting_ligands_receptors) > 0) {
      df_unique = df %>% dplyr::filter(!receptor %in% intersecting_ligands_receptors)
      df_duplicated = df %>% dplyr::filter(receptor %in% 
                                             intersecting_ligands_receptors)
      df_duplicated = df_duplicated %>% dplyr::mutate(receptor = paste(" ", 
                                                                       receptor, sep = ""))
      df = dplyr::bind_rows(df_unique, df_duplicated)
      intersecting_ligands_receptors = generics::intersect(unique(df$ligand), 
                                                           unique(df$receptor))
    }
    circos_links = df
    circos_links = circos_links %>% dplyr::inner_join(grid_col_tbl_ligand) %>% 
      dplyr::inner_join(grid_col_tbl_receptor)
    links_circle = circos_links %>% dplyr::distinct(ligand, 
                                                    receptor, weight)
    ligand_color = circos_links %>% dplyr::distinct(ligand, 
                                                    color_ligand_type)
    grid_ligand_color = ligand_color$color_ligand_type %>% 
      magrittr::set_names(ligand_color$ligand)
    receptor_color = circos_links %>% dplyr::distinct(receptor, 
                                                      color_receptor_type)
    grid_receptor_color = receptor_color$color_receptor_type %>% 
      magrittr::set_names(receptor_color$receptor)
    grid_col = c(grid_ligand_color, grid_receptor_color)
    ligand_order = prioritized_tbl_oi$sender %>% unique() %>% 
      sort() %>% lapply(function(sender_oi) {
        ligands = circos_links %>% dplyr::filter(sender == 
                                                   sender_oi) %>% dplyr::arrange(ligand) %>% dplyr::distinct(ligand)
      }) %>% unlist()
    receptor_order = prioritized_tbl_oi$receiver %>% unique() %>% 
      sort() %>% lapply(function(receiver_oi) {
        receptors = circos_links %>% dplyr::filter(receiver == 
                                                     receiver_oi) %>% dplyr::arrange(receptor) %>% 
          dplyr::distinct(receptor)
      }) %>% unlist()
    order = c(ligand_order, receptor_order)
    width_same_cell_same_ligand_type = 0.275
    width_different_cell = 3
    width_ligand_receptor = 9
    width_same_cell_same_receptor_type = 0.275
    sender_gaps = prioritized_tbl_oi$sender %>% unique() %>% 
      sort() %>% lapply(function(sender_oi) {
        sector = rep(width_same_cell_same_ligand_type, times = (circos_links %>% 
                                                                  dplyr::filter(sender == sender_oi) %>% dplyr::distinct(ligand) %>% 
                                                                  nrow() - 1))
        gap = width_different_cell
        return(c(sector, gap))
      }) %>% unlist()
    sender_gaps = sender_gaps[-length(sender_gaps)]
    receiver_gaps = prioritized_tbl_oi$receiver %>% unique() %>% 
      sort() %>% lapply(function(receiver_oi) {
        sector = rep(width_same_cell_same_receptor_type, 
                     times = (circos_links %>% dplyr::filter(receiver == 
                                                               receiver_oi) %>% dplyr::distinct(receptor) %>% 
                                nrow() - 1))
        gap = width_different_cell
        return(c(sector, gap))
      }) %>% unlist()
    receiver_gaps = receiver_gaps[-length(receiver_gaps)]
    gaps = c(sender_gaps, width_ligand_receptor, receiver_gaps, 
             width_ligand_receptor)
    if (length(gaps) != length(union(circos_links$ligand, 
                                     circos_links$receptor) %>% unique())) {
      warning("Specified gaps have different length than combined total of ligands and receptors - This is probably due to duplicates in ligand-receptor names")
    }
    links_circle$weight[links_circle$weight == 0] = 0.01
    circos.clear()
    circos.par(gap.degree = gaps)
    chordDiagram(links_circle, directional = 1, order = order, 
                 link.sort = TRUE, link.decreasing = TRUE, grid.col = grid_col, 
                 diffHeight = 0.0075, direction.type = c("diffHeight", 
                                                         "arrows"), link.visible = links_circle$weight > 
                   0.01, annotationTrack = "grid", preAllocateTracks = list(track.height = 0.175), 
                 grid.border = "gray35", link.arr.length = 0.05, link.arr.type = "big.arrow", 
                 link.lwd = 1.25, link.lty = 1, link.border = "gray35", 
                 reduce = 0, scale = TRUE)
    circos.track(track.index = 1, panel.fun = function(x, 
                                                       y) {
      circos.text(CELL_META$xcenter, CELL_META$ylim[1], 
                  CELL_META$sector.index, facing = "clockwise", 
                  niceFacing = TRUE, adj = c(0, 0.5), cex = 0.5)
    }, bg.border = NA)
    title(title)
    p_circos = recordPlot()
    return(p_circos)
  })
  names(all_plots) = groups_oi
  plot(NULL, xaxt = "n", yaxt = "n", bty = "n", ylab = "", 
       xlab = "", xlim = 0:1, ylim = 0:1)
  legend = ComplexHeatmap::Legend(at = prioritized_tbl_oi$receiver %>% 
                                    unique() %>% sort(), type = "grid", legend_gp = grid::gpar(fill = colors_receiver[prioritized_tbl_oi$receiver %>% 
                                                                                                                        unique() %>% sort()]), title_position = "topleft", title = "Receiver")
  ComplexHeatmap::draw(legend, just = c("left", "bottom"))
  legend = ComplexHeatmap::Legend(at = prioritized_tbl_oi$sender %>% 
                                    unique() %>% sort(), type = "grid", legend_gp = grid::gpar(fill = colors_sender[prioritized_tbl_oi$sender %>% 
                                                                                                                      unique() %>% sort()]), title_position = "topleft", title = "Sender")
  ComplexHeatmap::draw(legend, just = c("left", "top"))
  p_legend = grDevices::recordPlot()
  all_plots$legend = p_legend
  return(all_plots)
}
<bytecode: 0x555bd326d138>
  <environment: namespace:multinichenetr>