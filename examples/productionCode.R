library(maxabpp)
library(tidyverse)
library(stringdist)
library(stringr)
library(ggplot2)
library(ggforce)
library(ggrepel)
library(grid)
library(RColorBrewer)
library(plyr)


output1 <- pairwise_LFQ(
  raw = merge2,
  metadata = metadata,
  name_probe_mod = c("sh(all","sw(all)","sq(all)"),
  max_each_mod = 1,
  max_total_mods = 1,
  quantitation_level = "peptide",
  background_check = FALSE)
output2 <- append_ec_sites(output1, quantitation_level = "peptide")
plot_volcano(output2, "Orlistat10 _vs_ Orlistat1 _log2fold_change" , "Orlistat10 _vs_ Orlistat1 _-log10p-value", xlim = c(-5.5, 3), ylim = c(0, 5), "Gene Names", 3, 0, "Example Inhibitor/Example Probe", TRUE)

#Read MeLac-alkyne data
raw1 <-read_tsv("sw.txt")
raw2 <-read_tsv("sq.txt")
raw3 <-read_tsv("sh.txt")
common1 <- intersect(colnames(raw1), colnames(raw2))
common2 <- intersect(common1, colnames(raw3))
merge1 <- merge(raw1, raw2, by = common1, all = TRUE)
merge2 <- merge(merge1, raw3, by = common2, all = TRUE)
metadata <- read_tsv("metadata.txt")
peptide_lvl_table <- proc_mspTable(raw = merge2, metadata = metadata,
                                                          name_probe_mod = c("sh(all","sw(all)","sq(all)"),
                                                 max_each_mod = 1, max_total_mods = 1, quantitation_level = "peptide", background_check = FALSE, normalize_to = NULL)


all_peptides <- subset(filter(merge2, is.na(Reverse)) %>% select("Sequence", "Modifications", "Proteins", "Gene Names", "Protein Names", starts_with("Intensity ")))
#proteinGroups <- read_proteinGroups()
setList <- make_geneNames_setList()
setList[["peptide_lvl"]] <- as.list(peptide_lvl_table[["Gene Names"]])
plot_Max_Venn(Max_Venn(setList, IndividualAnalysis = FALSE))
tiers <- make_tiers2(setList)
tiers[["peptide_lvl"]] <- unique(all_peptides[["Gene Names"]])
plot_Max_Venn(Max_Venn(tiers, IndividualAnalysis = FALSE))


plot_target(tiers, density = 500)




