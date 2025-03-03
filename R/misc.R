
# Levels ------------------------------------------------------------------

sample_levels <- c(
  "LFD_R1",
  "LFD_R2",
  "LFD_R3",
  "HFD_1_R1",
  "HFD_1_R2",
  "HFD_3_R1",
  "HFD_3_R2"
)


condition_levels <- c("LFD", "HFD_1", "HFD_3")

manual_anno_levels = c("Beta", "Alpha", "Delta", "Gamma", "Acinar",
                       "Endothelial", "Stellate", "Immune")
# Colors ------------------------------------------------------------------


sample_color <- c(  "LFD_R1" = "#004B7A",
                    "LFD_R2" = "#0484d4",
                    "LFD_R3" = "#40aff5",
                    "HFD_1_R1" = "#c71500",
                    "HFD_1_R2" = "#fa1f05",
                    "HFD_3_R1" = "#fcb21c",
                    "HFD_3_R2" = "#FFD700")
condition_color <- c(
  "LFD" = "#004B7A",
  "HFD_1" = "#c71500",
  "HFD_3" = "#fcb21c"
)

cluster_color <- c(
  "0" = "#FA8231",
  "1" = "#FABA3E",
  "2" = "#1B8235",
  "3" = "#1977FA",
  "4" = "#FA3C25",
  "5" = "#8F426D",
  "6" = "#8EE8D7",
  "7" = "#004B7A",
  "8" = "#89C75F",
  "9" = "#ff00ff",
  "10" = "#44B7C2",
  "11" ="#ba55d3",
  "12" = "#00ff00",
  "13" = "#d8a767",
  "14" = "#A83708",
  "15" = "#87cefa",
  "16" = "#30C28F",
  "17" = "grey",
  "18" = "#8A9FD1",
  "19" = "#9983BD",
  "20" = "#ffec00",
  "21" = "#e6c2dc",
  "22" = "#f37b7d",
  "23" = "#2f4f4f",
  "24" = "#8b4513",
  "25" = "#b8860b",
  "26" = "#dc143c",
  "27" = "#00bfff",
  "28" = "#ff7f50",
  "29" = "purple",
  "30" = "#7fffd4",
  "31" = "#ff4500",
  "32" = "#008b8b",
  "33" = "#c0c0c0",
  "34" = "#0000ff",
  "35" = "#00ced1",
  "36" = "#ff1493",
  "37" = "#9acd32",
  "38" = "#cd5c5c",
  "39" = "#ffe4c4",
  "40" = "#ffb6c1",
  "41" = "#1e90ff",
  "42" = "#8fbc8f",
  "43" = "#bdb76b",
  "44" = "#2e8b57",
  "45" = "#ffa500")

member_color <- c("1" = "#A83708",
                  "2" = "#004B7A",
                  "3" = "#2e8b57",
                  "4" = "#F8AD4B",
                  "5" = "#FADEB8")

cluster_anno <- c("Beta" = "#A83708",
                  "Alpha" = "#2e8b57",
                  "Delta" = "#004B7A",
                  "Gamma" = "#F8AD4B",
                  "Acinar" = "#1e90ff",
                  "Stellate" = "#ff1493",
                  "Endothelial" = "#ba55d3",
                  "Immune" = "#00ced1")
# Marker genes ------------------------------------------------------------
# list of markers short
markers_short <- list(
  "Beta" = c("Ins1", "Ins2", "Slc2a2"),
  "Alpha" = c("Gcg", "Slc7a2", "Ttr"),
  "Delta" = c("Sst", "Rbp4", "Ghsr"),
  "Gamma" = c("Ppy", "Etv1", "Arx"),
  "Acinar" = c("Prss2", "Cpa1", "Pnlip"),
  "Endothelial" = c("Plvap", "Flt1", "Esm1"),
  "P_stellate" = c("Rgs5", "Col1a2", "Col6a1"),
  "Immune" = c("Cd74", "Cd83", "Cd53"))


markers_short_tn <- list(
  "Beta" = c("Ins2"),
  "Alpha" = c("Gcg"),
  "Delta" = c("Sst"),
  "Gamma" = c("Ppy"),
  "Cycling" = c("Top2a"),
  "Acinar" = c("Prss2"),
  "Endothelial" = c("Plvap"),
  "Q_stellate" = c("Rgs5"),
  "A_stellate" = c("Col1a1"),
  "Immune" = c("Cd74"),
  "Ductal" = c("Krt19"))
