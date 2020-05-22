library(maxabpp)
library(tidyverse)
library(stringdist)
library(rscripting)
library(stringr)
library(ggplot2)
library(ggforce)
library(ggrepel)
library(grid)
library(RColorBrewer)
library(plyr)

setList <- make_proteinGroups_setList()
plot_Max_Venn(Max_Venn(setList, IndividualAnalysis = TRUE))
tiers <- make_tiers(setList)
plot_target(tiers, density = 500)

output1 <- pairwise_LFQ(
  raw = sh,
  metadata = metadata,
  name_probe_mod = c("sh(all)"),
  max_each_mod = 1,
  max_total_mods = 1,
  quantitation_level = "peptide",
  background_check = FALSE)
output2 <- append_ec_sites(output1, quantitation_level = "peptide")
plot_volcano(output2, "Orlistat10 _vs_ Orlistat1 _log2fold_change" , "Orlistat10 _vs_ Orlistat1 _-log10p-value", xlim = c(-5.5, 3), ylim = c(0, 5), "Gene.Names", 3, 0, "Example Inhibitor/Example Probe")