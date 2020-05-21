#' @export
clean_list <- function (L) {
  L[] <- lapply(L, function(x) x[!x %in% ""])
  return (L)
}
#'@title Function read_proteinGroups
#'@description Read all .txt files (proteinGroups.txt renamed according to a specific experimental condition) a user-created folder.
#'@param folderName Default folder name is "proteinGroups", the fold must exist in working directory, use \code{\link{getwd}}() and \code{\link{dir}}() to check
#'@param selectColumnNames Default selected column names are
#'"Protein IDs",                          "Majority protein IDs",                 "Peptide counts (all)",
#'"Peptide counts (razor+unique)",        "Peptide counts (unique)",              "Protein names",
#'"Gene names",      "Number of proteins", "Peptides",                             "Razor + unique peptides",        "Unique peptides",
#'"Sequence coverage [%]",                "Unique + razor sequence coverage [%]", "Unique sequence coverage [%]",
#'"Mol. weight [kDa]",                    "Sequence length",                      "Sequence lengths",
#'"Q-value",                              "Score",                                    "Intensity",               "MS/MS count",
#'"Only identified by site",              "Potential contaminant"
#'@return a list of dataframes (each dataframe contains data from each file)
#'@export
read_proteinGroups <- function(folderName = "proteinGroups", selectColumnNames = defaultColumnNames) {
  paths <- dir(path = folderName, pattern = "\\.txt")
  names(paths) <- basename(paths)
  len <- length(paths)
  datasets <- sapply(1:len,function(x)NULL)
  names(datasets) <- names(paths)
  defaultColumnNames <- c(
    "Protein IDs",                          "Majority protein IDs",                 "Peptide counts (all)",
    "Peptide counts (razor+unique)",        "Peptide counts (unique)",              "Protein names",
    "Gene names",      "Number of proteins", "Peptides",                             "Razor + unique peptides",        "Unique peptides",
    "Sequence coverage [%]",                "Unique + razor sequence coverage [%]", "Unique sequence coverage [%]",
    "Mol. weight [kDa]",                    "Sequence length",                      "Sequence lengths",
    "Q-value",                              "Score",                                    "Intensity",               "MS/MS count",
    "Only identified by site",              "Potential contaminant")
  for (n_sets in 1:len)
  {
    dataset <- ldply(paste0("proteinGroups/",paths[[n_sets]]), read_tsv, col_types = cols(
      `Protein IDs` = col_character(),
      `Majority protein IDs` = col_character(),
      `Peptide counts (all)` = col_character(),
      `Peptide counts (razor+unique)` = col_character(),
      `Peptide counts (unique)` = col_character(),
      `Protein names` = col_character(),
      `Gene names` = col_character(),
      `Fasta headers` = col_character(),
      `Number of proteins` = col_integer(),
      Peptides = col_integer(),
      `Razor + unique peptides` = col_integer(),
      `Unique peptides` = col_integer(),
      `Sequence coverage [%]` = col_double(),
      `Unique + razor sequence coverage [%]` = col_double(),
      `Unique sequence coverage [%]` = col_double(),
      `Mol. weight [kDa]` = col_double(),
      `Sequence length` = col_integer(),
      `Sequence lengths` = col_character(),
      `Q-value` = col_double(),
      Score = col_double(),
      Intensity = col_double(),
      `MS/MS count` = col_integer(),
      `Only identified by site` = col_character(),
      Reverse = col_character(),
      `Potential contaminant` = col_character(),
      id = col_double(),
      `Peptide IDs` = col_character(),
      `Peptide is razor` = col_character(),
      `Mod. peptide IDs` = col_character(),
      `Evidence IDs` = col_character(),
      `MS/MS IDs` = col_character(),
      `Best MS/MS` = col_character(),
      `Oxidation (M) site IDs` = col_character(),
      `Oxidation (M) site positions` = col_character()))
    datasets[[n_sets]] <- subset(dataset, is.na(Reverse), select = selectColumnNames)
  }
  return(datasets)
}

#'@title Function make_proteinGroups_setList
#'@description Convert a list of datasets generated from \code{\link{read_proteinGroups}} to a setList of protein groups to be used for Venn Diagram \code{\link{Max_Venn}}()  or Tier \code{\link{make_tiers}}()  analysis
#'@export
make_proteinGroups_setList <- function (datasets = read_proteinGroups()) {
  setList <- sapply(1:length(datasets),function(x)NULL)
  names(setList) <- str_replace_all(names(datasets), "\\.txt", "")
  for (i in 1:length(datasets)) {
    setList[[i]] <- clean_list(datasets[[i]]$`Protein IDs`)
  }
  return(setList)
}

#'@title Function make_geneNames_setList
#'@description Convert a list of datasets generated from \code{\link{read_proteinGroups}} to a setList of protein groups to be used for Venn Diagram \code{\link{Max_Venn}}()  or Tier \code{\link{make_tiers}}()  analysis
#'@export
make_geneNames_setList <- function (datasets = read_proteinGroups()) {
  setList <- sapply(1:length(datasets),function(x)NULL)
  names(setList) <- str_replace_all(names(datasets), "\\.txt", "")
  for (i in 1:length(datasets)) {
    setList[[i]] <- clean_list(datasets[[i]]$`Gene names`)
  }
  return(setList)
}

#'@title Function make_proteinNames_setList
#'@description Convert a list of datasets generated from \code{\link{read_proteinGroups}} to a setList of protein groups to be used for Venn Diagram \code{\link{Max_Venn}}()  or Tier \code{\link{make_tiers}}()  analysis
#'@export
make_proteinNames_setList <- function (datasets = read_proteinGroups()) {
  setList <- sapply(1:length(datasets),function(x)NULL)
  names(setList) <- str_replace_all(names(datasets), "\\.txt", "")
  for (i in 1:length(datasets)) {
    setList[[i]] <- clean_list(datasets[[i]]$`Protein names`)
  }
  return(setList)
}

#'@title Function make_majorityproteinIDs_setList
#'@description Convert a list of datasets generated from \code{\link{read_proteinGroups}} to a setList of protein groups to be used for Venn Diagram \code{\link{Max_Venn}}()  or Tier \code{\link{make_tiers}}()  analysis
#'@export
make_majorityproteinIDs_setList <- function (datasets = read_proteinGroups()) {
  setList <- sapply(1:length(datasets),function(x)NULL)
  names(setList) <- str_replace_all(names(datasets), "\\.txt", "")
  for (i in 1:length(datasets)) {
    setList[[i]] <- clean_list(datasets[[i]]$`Majority protein IDs`)
  }
  return(setList)
}

#'@title Function make_custom_setList
#'@param columnName a string that specifies the name of the column, which user wants to keep for making the setList.
#'@description Convert a list of datasets generated from \code{\link{read_proteinGroups}} to a setList of protein groups to be used for Venn Diagram \code{\link{Max_Venn}}()  or Tier \code{\link{make_tiers}}()  analysis
#'@export
make_custom_setList <- function (datasets = read_proteinGroups(), columnName) {
  setList <- sapply(1:length(datasets),function(x)NULL)
  names(setList) <- str_replace_all(names(datasets), "\\.txt", "")
  for (i in 1:length(datasets)) {
    setList[[i]] <- clean_list(datasets[[i]][[columnName]])
  }
  return(setList)
}
