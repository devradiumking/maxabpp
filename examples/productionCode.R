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


LFQ_table_MeLac_noNorm <- pairwise_LFQ(
  raw = test,
  metadata = metadata,
  name_probe_mod = c("sh(all)","sw(all)","sq(all)"),
  max_each_mod = 1,
  max_total_mods = 1,
  quantitation_level = "peptide",
  background_check = FALSE)
LFQ_table_MeLac_noNorm <- append_ec_sites(LFQ_table_MeLac_noNorm, quantitation_level = "peptide")

LFQ_table_MeLac_normtoAll <- pairwise_LFQ(
  raw = test,
  metadata = metadata,
  name_probe_mod = c("sh(all)","sw(all)","sq(all)"),
  max_each_mod = 1,
  max_total_mods = 1,
  quantitation_level = "peptide",
  background_check = FALSE,
  normalize_to = "sum_all")
LFQ_table_MeLac_normtoAll <- append_ec_sites(LFQ_table_MeLac_normtoAll, quantitation_level = "peptide")

plot_volcano(output2, "Parthenolide10 _vs_ Parthenolide1 _log2fold_change" , "Parthenolide10 _vs_ Parthenolide1 _-log10p-value", xlim = c(-6, 1), ylim = c(0, 5), "Gene Names", 1.3, -1, "Parthenolide 10 uM vs 1 uM/MeLac-alkyne Probe", TRUE)

#Read MeLac-alkyne data
raw1 <-read_tsv("sw.txt")
raw2 <-read_tsv("sq.txt")
raw3 <-read_tsv("sh.txt")
common_colnames <- intersection(colnames(raw1), colnames(raw2), colnames(raw3))
raw1 <- subset(raw1) %>% select(common_colnames)
raw2 <- subset(raw2) %>% select(common_colnames)
raw3 <- subset(raw3) %>% select(common_colnames)
test <- rbind(raw1, raw2, raw3)


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




