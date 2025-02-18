
# Sample levels -----------------------------------------------------------

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
