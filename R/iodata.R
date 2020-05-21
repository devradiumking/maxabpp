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
    dataset <- ldply(paste0("proteinGroups/",paths[[n_sets]]), read_tsv)
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
