#' Pairwise label-free quantitation: convert MaxQuant modificationSpecificPeptides.txt file to a volcano plot-ready table
#' @param raw                 a dataframe by reading modificationSpecificPeptides.txt
#' @param metadata            a dataframe that maches the MaxQuant input. Column 1: Intensity (such as Intensity samplename, same as the column names in modificationSpecificPeptides.txt) name Column 2: Replicate group (use the same name for each group of replicates)
#' @param name_probe_mod      a string vector of chemical probe/modification names, such as c("Mod1", "Mod2"), must match MaxQuant input
#' @param max_each_mod        a integer as the maximal number of modifications on a single peptide, set for each chemical probe
#' @param max_total_mods      a integer as the maximal number of modifications on a single peptide, set for all chemical probes Note max_each_mod must not be less than max_total_mods
#' @param quantitation_level  a string, must be either "peptide" or "protein"
#' @param background_check    a boolean, FALSE = quantify probe-modified peptides, TRUE = quantify non-probe-modified peptides
#' @param normalize_to        a string, must be either "sum_all", "mean_all", (normalize to all peptides) "sum_background", or "mean_background" (normalize to background/non-probe-modified peptides).
#' @return a dataframe table of peptides, proteins and their corresponding -log10(p-value) and log2(FC)
#' @examples output <- pairwise_LFQ(raw = read.delim("modificationSpecificPeptides.txt", header=TRUE, sep="\t"), metadata = read.delim("metadata.txt", header=TRUE, sep="\t"), name_probe_mod = c("Mod1", "Mod2"), max_each_mod = 1, max_total_mods = 1, quantitation_level = "peptide" , background_check = FALSE, normalize_to = "mean_all")
#' @export
pairwise_LFQ <- function (raw = read_tsv("modificationSpecificPeptides.txt"), metadata = read_tsv("metadata.txt"),
name_probe_mod, max_each_mod = 1, max_total_mods = 1, quantitation_level = "peptide" , background_check = FALSE, normalize_to = NULL) {

#* Establish a new 'ArgCheck' object
Check <- ArgumentCheck::newArgCheck()
num_mods <- length(name_probe_mod)

#* Add an error if raw or metadata is missing
if (length(metadata[[1]]) != length(intersect(metadata[[1]], names(raw))))
  ArgumentCheck::addError(
    msg = "'Bad metadata! The metadata does not match the MaxQuant output. Doube check the MaxQuant output Intensity xxx column names and make sure the first column of metadata contains all of them",
    argcheck = Check
  )


#* Add an error if raw or metadata is missing
if (is.null(raw))
  ArgumentCheck::addError(
    msg = "'raw' is missing",
    argcheck = Check
  )
if (is.null(metadata))
  ArgumentCheck::addError(
    msg = "'metadata' is missing",
    argcheck = Check
  )
#* Add an error if number of input mods is less than 1
if (num_mods < 1 || num_mods > 10 || is.character(name_probe_mod) == FALSE)
  ArgumentCheck::addError(
    msg = "'name_probe_mod' must be a vector containing 1 to 10 string components",
    argcheck = Check
  )

#* Return errors and warnings (if any)
ArgumentCheck::finishArgCheck(Check)


#Internal function 1: appending a column named "parent sequence" of the longest miscleaved version of each peptide sequence
#If one cleaved peptide can map to multiple miscleaved peptides, the "parent sequence" will be a concatenated string of all
add_parent_sequence <- function (dataset) {
  seq <- dataset$Sequence
  dataset$parent_sequence <- ""
  for (i in 1:length(seq)) {
    related_seq <- unique(as.character(subset(seq, grepl(seq[i],seq))))
    length_seq <- nchar(related_seq)
    num_seq <- length(length_seq)
    if (num_seq == 1) {
      dataset$parent_sequence[i] <- as.character(related_seq)
    } else if (num_seq == 2) {
      dataset$parent_sequence[i] <- as.character(related_seq[2])
    } else {
      parents <- subset(related_seq, nchar(related_seq) == max(length_seq))
      if (length(parents) == 1) {
        dataset$parent_sequence[i] <- parents
      } else {
        dataset$parent_sequence[i] <- paste(parents, collapse = ";")
      }
    }
  }
  return(dataset)
}

#Internal function 2: replacing column names with metadata replicate groups
replace_column_names <- function (interim_data, metadata) {
  num_column_rename <- length(colnames(interim_data %>%  select(starts_with("Intensity "))))
  if (nrow(metadata) == num_column_rename) {
    for (k in 1:num_column_rename) {
      names(interim_data)[names(interim_data) == metadata[[1]][k]] <- metadata[[2]][k]
    }
    return(interim_data)
  } else {
    return(NULL)
  }
}

#Internal function 3: replace characters in a data frame
replace_val <- function (data, unwanted_value, replacement) {
  data <- apply(data,2,function(x)gsub(unwanted_value, replacement,x))
return(data)
}

#Internal function 4A: determine if two groups of measurements have the same variance using F-test
if_equal_variance <- function (group1, group2) {
  if (mean(group1) == 0 && mean(group2) == 0) {
    return(TRUE)
  } else {
    F_result <- var.test(group1, group2, ratio = 1,
                         alternative = c("two.sided"),
                         conf.level = 0.95)
    if (F_result$p.value <= 0.05) {
      return(FALSE)
    } else {
      return(TRUE)
    }
  }
}

#Internal function 4B: generate a table for volcano plot
generate_volcano_table <- function (interim_data) {
  result <- data.frame()
  result_p <- data.frame()
  result_fc <- data.frame()
  #interim_data_select <- subset(interim_data, select = -c(2,3,4,5))
  for (j in 1:nrow(interim_data)) {
    entry_data <- as.data.frame(t(interim_data[j,]))
    entry_data <- rownames_to_column(entry_data, "Group")
    entry_data <- separate(entry_data, Group, c("remain","trim"), sep = "\\.")
    for (k in 1:length(filtered_pairs)) {
      #Define group2 as the inhibitor-treated at a higher inhibitor concentration
      entry_data_group1 <- as.numeric(as.vector(subset(entry_data, remain == filtered_pairs$V1[k], stringasfactor = FALSE)[,3]))
      entry_data_group2 <- as.numeric(as.vector(subset(entry_data, remain == filtered_pairs$V2[k], stringasfactor = FALSE)[,3]))
      #Under the presumption that the higher-concentration inhibitor-treated group2 should have a lower mean intensity
      #Null hypothesis: group1 and group2 have the same intensity
      #Alternative hypothesis, group1 has greater intensity
      t_result <- NULL
       t_result <- t.test(x = entry_data_group1, y = entry_data_group2,
                        alternative = c("greater"),
                        #Call function 4A to determine variance equality
                        mu = 0, paired = FALSE, var.equal = if_equal_variance(entry_data_group1, entry_data_group2),
                        conf.level = 0.95)
      p_value <- t_result$p.value
      mean1 <- as.numeric(t_result$estimate[1])
      mean2 <- as.numeric(t_result$estimate[2])
      fold_change <- mean2/mean1
      result_p[j,k] <- -log10(p_value)
      result_fc[j,k] <- log2(fold_change)
      colnames(result_p)[k] <- paste(filtered_pairs$V2[k],"_vs_", filtered_pairs$V1[k],"_-log10p-value")
      colnames(result_fc)[k] <- paste(filtered_pairs$V2[k],"_vs_", filtered_pairs$V1[k],"_log2fold_change")
    }
    result <- cbind(result_p, result_fc)
  }
  if (colnames(interim_data)[1] == "parent_sequence") {
    result <- cbind(interim_data[,1:6], result)
  } else {
    result <- cbind(interim_data[,1:5], result)
  }
  return (result)
}

#Internal function 5A: Recursively calculate the total mods
criteria_total <- function (dataframe, probe_mods) {
  type_mods <- length(probe_mods)
  m <- dataframe[[probe_mods[type_mods]]]
  if (type_mods == 1) {
    return(m)
  } else {
    probe_mods <- head(probe_mods, type_mods - 1)
    return(m + criteria_total(dataframe, probe_mods))
  }
}

#Internal function 5B: Filter dataset according to the numbers of max mods
subset_mod <- function (dataframe, probe_mods, max_each_mod, max_total_mods) {
  type_mods <- length(probe_mods)
  criteria_total <- criteria_total(dataframe, probe_mods)
  subset1 <- dataframe[criteria_total <= max_total_mods , ]
  for (i in 1:type_mods) {
    criteria_single <- subset1[[probe_mods[i]]]
    subset1 <- subset1[criteria_single <= max_each_mod , ]
  }
  return(subset1)
}

#Interanl function 6: Normalize to background
normalize_to_background <- function (data, background) {
  intensity_colnames <- intersect(colnames(data), colnames(background))
  for (i in 1:length(intensity_colnames)) {
    data[intensity_colnames[i]] <- data[intensity_colnames[i]]/as.numeric(background[intensity_colnames[i]])
  }
  return(data)
}

#Interanl function 7: Filter dataset to meet mod requirements
#Buggy!
elemental_mod_subset <- function (dataset, probe_mods, max_each_mod, max_total_mods) {
  #Create an tibble with same column names and first row of data from the input dataset
  names <- names(dataset)
  filtered_dataset <- tibble() %>% tibble_add_column(dataset[1,])
  #This step removes first row of data and creates an empty tibble with same column names as the input dataset
  filtered_dataset <- filtered_dataset[-1,]
  for (row_num in 1:nrow(dataset)) {
    cell_must_contain <- str_split(dataset[["Modifications"]], "\\;")[[row_num]]
#Bad coding!
#"Carbamidomethyl (C)" "2 sw(all)"
#detected_mods <- intersect(cell_must_contain, name_probe_mod)
#detected_mods
#character(0)
    detected_mods <- subset(cell_must_contain, mgrepl(probe_mods, cell_must_contain, op = "|"))
    #Bad coding! str_extract failed to include number substring
    mod_counts <- str_extract(detected_mods, "[1-9]")
    if (length(mod_counts) == 0) {
      mod_counts <- 1
    }
    if (is.na(mod_counts)) {
      mod_counts <- as.integer()
      for (c in 1:length(detected_mods)) {
        mod_counts[c] <- 1
      }
    } else {mod_counts <- as.integer(mod_counts)}
    condition1 <- as.integer(sum(mod_counts) <=  max_total_mods)
    condition2 <- prod(as.integer((mod_counts <= max_each_mod)))
    if (length(detected_mods) > 0 & condition1*condition2 == 1) {
      filtered_dataset <- rbind(filtered_dataset, dataset[row_num,])
    }
  }
  return(filtered_dataset)
}

#Obsolete, use readr instead of read.delim (Replace space and dash characters with dots due to the limitation of column renaming in R)
#metadata <- replace_val(metadata,'\\s+', '.')
#metadata <- replace_val(metadata,'\\-', '.')
sample_groups <- as.character(unique(metadata[,2])[[1]])
#Replace special characters with dots within the probe modification string due to the limitation of column renaming in R
#probe_mods <- name_probe_mod
#probe_mods <- gsub("-", ".", gsub(" ", ".", gsub("[()]", ".", name_probe_mod)))
probe_mods <- str_replace_all(name_probe_mod,"\\(","\\\\(")
probe_mods <- str_replace_all(probe_mods,"\\)","\\\\)")
##Estimate background using non-probe-modified peptides
background_peptides <- subset(filter(raw, !mgrepl(probe_mods, Modifications, op ="|")), is.na(Reverse)) %>% select("Sequence", "Modifications", "Proteins", "Gene Names", "Protein Names", starts_with("Intensity "))
all_peptides <- subset(filter(raw, is.na(Reverse)) %>% select("Sequence", "Modifications", "Proteins", "Gene Names", "Protein Names", starts_with("Intensity ")))
sum_background <- summarise_all(background_peptides %>% select(starts_with("Intensity ")), sum)
mean_background <- summarise_all(background_peptides %>% select(starts_with("Intensity ")), mean)
#mad_background <- summarise_all(background_peptides %>% select(starts_with("Intensity ")), mad)
sum_all <- summarise_all(all_peptides %>% select(starts_with("Intensity ")), sum)
mean_all <- summarise_all(all_peptides %>% select(starts_with("Intensity ")), mean)
#mad_all <- summarise_all(all_peptides %>% select(starts_with("Intensity ")), mad)
###Conditional statement for the quantititaion of non-probe-modified peptides
if (background_check == TRUE) {
  filtered_step2 <- add_parent_sequence(background_peptides)
  mod_validation <- TRUE
} else {
###Conditional statement probe modification-specific quantitation
 #Extract a subset containing non-reverse peptides carrying probe modifications
 filtered_step1 <- subset(filter(subset(raw, mgrepl(probe_mods, raw[["Modifications"]], op = "|")), is.na(Reverse)))
 #Validate user defined probe modification
 column_names <- colnames(raw)
 #mod_col_index <- match(probe_mods,column_names,nomatch = 0)
 mod_validation <- TRUE
 #Filter out peptides carrying more or less than specified numbers of probe modifications on a single peptide
 #filtered_step2 <- subset_mod(filtered_step1, probe_mods, max_each_mod, max_total_mods)
 filtered_step2 <- elemental_mod_subset(filtered_step1, probe_mods, max_each_mod, max_total_mods)
 #Call function 2
 filtered_step2 <- add_parent_sequence(filtered_step2)
}
#Compare string similarity of sample groups and generate t.test pairs automatically
pairs <- t(combn(sample_groups,2))
V3 <- stringsim(pairs[,1],pairs[,2])
pairs <- cbind(pairs,V3)
filtered_pairs <- filter(as.data.frame(pairs,stringsAsFactors = FALSE), V3 > 0.8)
filtered_pairs <- filter(filtered_pairs, V3 < 1)
#Swap values so column V1 contains only shorter values
b <- str_length(filtered_pairs$V1) < str_length(filtered_pairs$V2)
for (i in 1:length(b)) {
  longer_value <- filtered_pairs$V1[i]
  shorter_value <- filtered_pairs$V2[i]
  if (b[i] == FALSE) {
    filtered_pairs$V1[i] <- shorter_value
    filtered_pairs$V2[i] <- longer_value
  }
}
if (quantitation_level == "peptide" && mod_validation == TRUE) {
  ##Quantitation at peptide level
  #Compress/summarize intensity values of each peptide for each sample, extract columns containing sequence information and intensity data
  by_sequence <- filtered_step2 %>%  select("parent_sequence", "Sequence", "Modifications", "Proteins", "Gene Names", "Protein Names", starts_with("Intensity ")) %>% group_by(parent_sequence)
  #Conditional calling normalization function 6
  if (!is.null(normalize_to) && normalize_to == "sum_background") {
    by_sequence <- normalize_to_background(by_sequence, sum_background)
  } else if (!is.null(normalize_to) && normalize_to == "mean_background") {
    by_sequence <- normalize_to_background(by_sequence, mean_background)
  } else if (!is.null(normalize_to) && normalize_to == "sum_all") {
    by_sequence <- normalize_to_background(by_sequence, sum_all)
  } else if (!is.null(normalize_to) && normalize_to == "mean_all") {
    by_sequence <- normalize_to_background(by_sequence, mean_all)
  } else {
    #do nothing
  }
  intensity_subset1 <- by_sequence %>% summarise_if(is.numeric, sum)
  #Collapse name strings with ; as the delimiter
  name_subset1 <- by_sequence %>% summarise_at(c("Sequence","Modifications", "Proteins", "Gene Names", "Protein Names"), ~paste(.x, collapse="; "))
  #Combine two parts as a new dataset
  interim_peptide_data <- cbind(name_subset1, subset(intensity_subset1, select = -parent_sequence))
  #Call function 2
  interim_peptide_data <- replace_column_names(interim_peptide_data, metadata)
  #Generate peptide level output
  peptide_volcano_table <- generate_volcano_table(interim_peptide_data)
  return (peptide_volcano_table)
 } else if (quantitation_level == "protein" && mod_validation == TRUE) {
   ##Quantitation at protein level
   #Compress/summarize intensity values of each protein/protein group for each sample, extract columns containing sequence information and intensity data
   by_proteins <- filtered_step2 %>%  select("Sequence", "Modifications", "Proteins", "Gene Names", "Protein Names", starts_with("Intensity ")) %>% group_by(Proteins)
   if (!is.null(normalize_to) && normalize_to == "sum_background") {
     by_proteins <- normalize_to_background(by_proteins, sum_background)
   } else if (!is.null(normalize_to) && normalize_to == "mean_background") {
     by_proteins <- normalize_to_background(by_proteins, mean_background)
   } else if (!is.null(normalize_to) && normalize_to == "sum_all") {
     by_proteins <- normalize_to_background(by_proteins, sum_all)
   } else if (!is.null(normalize_to) && normalize_to == "mean_all") {
     by_proteins <- normalize_to_background(by_proteins, mean_all)
   } else {
     #do nothing
   }
   intensity_subset2 <- by_proteins %>% summarise_if(is.numeric, sum)
   #Collapse name strings with ; as the delimiter
   name_subset2 <- by_proteins %>% summarise_at(c("Modifications", "Sequence", "Gene Names", "Protein Names"), ~paste(.x, collapse="; "))
   #Combine two parts as a new dataset
   interim_protein_data <- cbind(name_subset2, subset(intensity_subset2, select = -Proteins))
   #Call function 2
   interim_protein_data <- replace_column_names(interim_protein_data, metadata)
   #Generate protein level output
   protein_volcano_table <- generate_volcano_table(interim_protein_data)
  return (protein_volcano_table)
 } else {
   return(NULL)
 }
}

