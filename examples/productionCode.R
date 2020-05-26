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
plot_volcano(LFQ_table_MeLac_noNorm, "Orlistat10 _vs_ Orlistat1 _log2fold_change" , "Orlistat10 _vs_ Orlistat1 _-log10p-value", xlim = c(-5.5, 2), ylim = c(0, 4), "Gene Names", 1.3, 0, "Orlistat 10 uM vs 1 uM/MeLac-alkyne Probe", FALSE)
write_csv(LFQ_table_MeLac_noNorm, "LFQ_table_MeLac_noNorm.csv")
plot_volcano(LFQ_table_MeLac_noNorm, "Parthenolide10 _vs_ Parthenolide1 _log2fold_change" , "Parthenolide10 _vs_ Parthenolide1 _-log10p-value", xlim = c(-5.5, 1), ylim = c(0, 4), "Gene Names", 1.3, 0, "Parthenolide 10 uM vs 1 uM/MeLac-alkyne Probe", FALSE)

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
plot_volcano(LFQ_table_MeLac_normtoAll, "Orlistat10 _vs_ Orlistat1 _log2fold_change" , "Orlistat10 _vs_ Orlistat1 _-log10p-value", xlim = c(-5.5, 2), ylim = c(0, 4), "Gene Names", 1.3, 0, "Orlistat 10 uM vs 1 uM/MeLac-alkyne Probe (Normalized)", FALSE)
write_csv(LFQ_table_MeLac_normtoAll, "LFQ_table_MeLac_normtoAll.csv")
plot_volcano(LFQ_table_MeLac_normtoAll, "Parthenolide10 _vs_ Parthenolide1 _log2fold_change" , "Parthenolide10 _vs_ Parthenolide1 _-log10p-value", xlim = c(-5.5, 2), ylim = c(0, 4), "Gene Names", 1.3, 0, "Parthenolide 10 uM vs 1 uM/MeLac-alkyne Probe (Normalized)", FALSE)

#Read MeLac-alkyne data
raw1 <-read_tsv("sw.txt")
raw2 <-read_tsv("sq.txt")
raw3 <-read_tsv("sh.txt")
common_colnames <- intersection(colnames(raw1), colnames(raw2), colnames(raw3))
raw1 <- subset(raw1) %>% select(common_colnames)
raw2 <- subset(raw2) %>% select(common_colnames)
raw3 <- subset(raw3) %>% select(common_colnames)
test <- rbind(raw1, raw2, raw3)
raw4 <- read_tsv("sGSH.txt")


LFQ_table_GSH_MeLac_noNorm <- pairwise_LFQ(
  raw = raw4,
  metadata = metadata,
  name_probe_mod = c("sGSH"),
  max_each_mod = 1,
  max_total_mods = 1,
  quantitation_level = "peptide",
  background_check = FALSE)
LFQ_table_GSH_MeLac_noNorm <- append_ec_sites(LFQ_table_GSH_MeLac_noNorm, quantitation_level = "peptide")
plot_volcano(LFQ_table_GSH_MeLac_noNorm, "Orlistat10 _vs_ Orlistat1 _log2fold_change" , "Orlistat10 _vs_ Orlistat1 _-log10p-value", xlim = c(-5.5, 2), ylim = c(0, 4), "Gene Names", 1.3, 0, "Orlistat 10 uM vs 1 uM/GSH-MeLac Probe", FALSE)
write_csv(LFQ_table_GSH_MeLac_noNorm, "LFQ_table_GSH_MeLac_noNorm.csv")
plot_volcano(LFQ_table_GSH_MeLac_noNorm, "Parthenolide10 _vs_ Parthenolide1 _log2fold_change" , "Parthenolide10 _vs_ Parthenolide1 _-log10p-value", xlim = c(-5.5, 1), ylim = c(0, 4), "Gene Names", 1.3, 0, "Parthenolide 10 uM vs 1 uM/GSH-MeLac Probe", FALSE)

LFQ_table_GSH_MeLac_normtoAll <- pairwise_LFQ(
  raw = raw4,
  metadata = metadata,
  name_probe_mod = c("sGSH"),
  max_each_mod = 1,
  max_total_mods = 1,
  quantitation_level = "peptide",
  background_check = FALSE,
  normalize_to = "sum_all")
LFQ_table_GSH_MeLac_normtoAll <- append_ec_sites(LFQ_table_GSH_MeLac_normtoAll, quantitation_level = "peptide")
plot_volcano(LFQ_table_GSH_MeLac_normtoAll, "Orlistat10 _vs_ Orlistat1 _log2fold_change" , "Orlistat10 _vs_ Orlistat1 _-log10p-value", xlim = c(-5.5, 2), ylim = c(0, 4), "Gene Names", 1.3, 0, "Orlistat 10 uM vs 1 uM/GSH-MeLac Probe (Normalized)", FALSE)
write_csv(LFQ_table_GSH_MeLac_normtoAll, "LFQ_table_GSH_MeLac_normtoAll.csv")
plot_volcano(LFQ_table_GSH_MeLac_normtoAll, "Parthenolide10 _vs_ Parthenolide1 _log2fold_change" , "Parthenolide10 _vs_ Parthenolide1 _-log10p-value", xlim = c(-5.5, 1), ylim = c(0, 4), "Gene Names", 1.3, 0, "Parthenolide 10 uM vs 1 uM/GSH-MeLac Probe", FALSE)

metadata <- read_tsv("metadata.txt")
peptide_lvl_table <- proc_mspTable(raw = merge2, metadata = metadata,
                                   name_probe_mod = c("sh(all","sw(all)","sq(all)"),
                                   max_each_mod = 1, max_total_mods = 1, quantitation_level = "peptide", background_check = FALSE, normalize_to = NULL)


MeLac_peptides <- subset(test) %>% select("Sequence", "Modifications", "Proteins", "Gene Names", "Protein Names", starts_with("Intensity "))
MeLac_mod_peptides <- elemental_mod_subset(MeLac_peptides, c("sh(all)","sw(all)","sq(all)"), 1, 1)
GSH_MeLac_peptides <- subset(raw4) %>% select("Sequence", "Modifications", "Proteins", "Gene Names", "Protein Names", starts_with("Intensity "))
GSH_MeLac_mod_peptides <- elemental_mod_subset(GSH_MeLac_peptides, c("sGSH"), 1, 1)
setList_GSH_GSHMeLac <- list()
setList_GSH_GSHMeLac[["MeLac-alkyne"]] <- as.list(unique(MeLac_mod_peptides$Sequence))
setList_GSH_GSHMeLac[["GSH-MeLac"]] <- as.list(unique(GSH_MeLac_mod_peptides$Sequence))
setList_GSH_GSHMeLac_miscleavageGrp <- list()
setList_GSH_GSHMeLac_miscleavageGrp[["MeLac-alkyne"]] <- as.list(unique(LFQ_table_MeLac_normtoAll$parent_sequence))
setList_GSH_GSHMeLac_miscleavageGrp[["GSH-MeLac"]] <- as.list(unique(LFQ_table_GSH_MeLac_normtoAll$parent_sequence))
setList_GSH_GSHMeLac_Venn <- Max_Venn(setList_GSH_GSHMeLac, IndividualAnalysis = FALSE)
setList_GSH_GSHMeLac_Venn_miscleavageGrp <- Max_Venn(setList_GSH_GSHMeLac_miscleavageGrp, IndividualAnalysis = FALSE)
plot_Max_Venn(setList_GSH_GSHMeLac_Venn)
plot_Max_Venn(setList_GSH_GSHMeLac_Venn_miscleavageGrp)

#proteinGroups <- read_proteinGroups()
setList <- make_proteinGroups_setList()
#setList[["peptide_lvl"]] <- as.list(peptide_lvl_table[["Gene Names"]])
plot_Max_Venn(Max_Venn(setList, IndividualAnalysis = FALSE))
tiers <- make_tiers(setList)
tiers[["peptide_lvl"]] <- unique(all_peptides[["Gene Names"]])
plot_Max_Venn(Max_Venn(tiers, IndividualAnalysis = FALSE))
plot_target(tiers, density = 500)
dftier <- list2df(tiers)
df_GSH_GSHMeLac_Venn_miscleavageGrp <- list2df(setList_GSH_GSHMeLac_Venn_miscleavageGrp@IntersectionSets)
write_csv(dftier, "3tiers.csv")
write_csv(LFQ_table_MeLac_normtoAll, "LFQ_table_MeLac_normtoAll.csv")
write_csv(LFQ_table_GSH_MeLac_noNorm, "LFQ_table_GSH_MeLac_noNorm.csv")
write_csv(LFQ_table_GSH_MeLac_normtoAll, "LFQ_table_GSH_MeLac_normtoAll.csv")
write_csv(df_GSH_GSHMeLac_Venn_miscleavageGrp, "df_GSH_GSHMeLac_Venn_miscleavageGrp.csv")
