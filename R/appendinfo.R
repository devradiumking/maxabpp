#' Append columns of enzyme class and active site information to the pairwise_LFQ output
#' @seealso \code{\link{pairwise_LFQ}}
#' @param LFQ_table     a dataframe saved as the pairwise_LFQ() output
#' @param quantitation_level  a string, must be either "peptide" or "protein", must be match the pairwise_LFQ argument
#' @return             a volcano plot-ready dataframe annotated with two extra columns of enzyme class and active sites
#' @examples output2 <- append_ec_sites(output1, quantitation_level = "peptide")
#' @export
append_ec_sites <- function(LFQ_table, quantitation_level) {
#proteome_database <-read.delim("homo_sapiens_ec_sites_database_05132020.txt", header=TRUE, sep="\t")
#use_data(proteome_database)
#Internal function 1: extract an ec number for the all proteins in a protein group from the database
database_search_ec <- function (protein_group) {
  protein_hits <- vector()
  ec <- vector()
  proteins <- str_split(protein_group, ";")[[1]]
  #Remove isoform-related suffix of each protein entry
  deisoformed_proteins <- str_replace(proteins, "\\-.*", "")
  ec <- vector()
  for (i in 1:length(deisoformed_proteins)) {
  index <- grep(deisoformed_proteins[i], proteome_database$Entry)
  each_ec <- as.character(proteome_database$"EC number"[index])
  each_ec <- str_trunc(each_ec, 1, "right", ellipsis = "")
  ec <- c(ec, each_ec)
  }
  ec <- subset(ec, ec != "")
  ec <- unique(ec)
  if (length(ec) > 0) {
  return(ec)
  } else {return ("na")}
}

#Internal function 2: translate ec numbers to the class names and attach
attach_enzyme_class <- function (LFQ_table) {
  len <- nrow(LFQ_table)
  ec_list <- vector()
  ec <- NA
  for (i in 1:len) {
    ec_group <- database_search_ec(LFQ_table$Proteins[i])
    ec_sublist <- vector()
    for (j in 1:length(ec_group)) {
    switch(ec_group[j],
           "1" = {ec_sublist[j] <- "Oxidoreductase"},
           "2" = {ec_sublist[j] <- "Transferase"},
           "3" = {ec_sublist[j] <- "Hydrolase"},
           "4" = {ec_sublist[j] <- "Lyase"},
           "5" = {ec_sublist[j] <- "Isomerase"},
           "6" = {ec_sublist[j] <- "Ligase"},
           "7" = {ec_sublist[j] <- "Translocase"},
           ec_sublist[j] <- ""
          )
    }
    if (length(unique(ec_sublist)) > 1) {
      ec_list[i] <- paste0(ec_sublist, collapse = "/")
    }
    else if (length(ec_sublist) == 1 && ec_sublist != "") {
      ec_list[i] <- ec_sublist
    } else {
      ec_list[i] <- "Others"
    }
  }
  return(cbind(LFQ_table, ec_list))
}

#Internal function 3: extract numeric site infomation
extract_sites <- function (column_of_sites, prefix) {
  test <- str_trim(str_split(column_of_sites, ";")[[1]])
  criteria_test <- test %in% grep(prefix, test, value = TRUE)
  return (parse_number(subset(test, criteria_test)))
}

#Internal function 4: extract a protein sequence for each protein in a protein group from the database
database_search_sites <- function (protein_group) {
  protein_hits <- vector()
  binding_sites <- vector()
  active_sites <- vector()
  other_sites <- vector()
  for (n in 1:length(protein_group)) {
    index <- grep(protein_group[n], proteome_database$Entry)
    protein_hits[n] <- as.character(proteome_database$Sequence[index])
    s1 <- extract_sites(proteome_database$"Binding site"[index], "BINDING")
    if (length(s1) > 0) {
      binding_sites[n] <- s1
    } else {
      binding_sites[n] <- "no"
    }
    s2 <- extract_sites(proteome_database$"Active site"[index], "ACT_SITE")
    if (length(s2) > 0) {
      active_sites[n] <- s2
    } else {
      active_sites[n] <- "no"
    }
    s3 <- extract_sites(proteome_database$"Site"[index], "SITE")
    if (length(s3) > 0) {
      other_sites[n] <- s3
    } else {
      other_sites[n] <- "no"
    }
  }
  return(cbind(protein_hits, binding_sites, active_sites, other_sites))
}

#Internal function 5: compare reported sites to tryptic peptides, get site coverage of each tryptic peptide
compare_sites <- function (site_positions, n_term, c_term) {
  num_sites <- length(site_positions)
  site <- NULL
  output <- NULL
  for (i in 1:num_sites) {
    site <- as.integer(site_positions[i])
    if (is.na(site) | is.na(n_term) | is.na(c_term) | n_term == c_term){
      output <- paste0(output, ";", "na")
    } else if (n_term <= site && c_term >= site) {
      output <- paste0(output, ";", site)
    } else {
      output <- paste0(output, ";", "na")
    }
  }
  return(str_sub(output, 2))
}

#Internal function 6: append site information
match_mod_peptides <- function (LFQ_table) {
  output <- data.frame(binding_sites = NULL, active_sites = NULL, other_sites = NULL)
  for (i in 1:nrow(LFQ_table)) {
    proteins_list <- LFQ_table$Proteins
    peptides_list <- LFQ_table$parent_sequence
    binding_sites <- data.frame()
    active_sites <- data.frame()
    other_sites <- data.frame()
    all_sites <- data.frame()
    empty_site_table <- data.frame(protein_hits = "", binding_sites = "", active_sites = "", other_sites = "")
    sites_table <- tryCatch(as.data.frame(database_search_sites(str_split(proteins_list[i], ";")[[1]])), error = function(c) empty_site_table)
    protein_sequence <- as.vector(sites_table[1,1])
    peptide_sequence <- peptides_list[i]
    mod_peptide_locations <- str_locate(protein_sequence, peptide_sequence)
    start <- as.integer(mod_peptide_locations[1])
    end <- as.integer(mod_peptide_locations[2])
    binding_site_positions <- as.integer(paste(sites_table[,2]))
    active_site_positions <- as.integer(paste(sites_table[,3]))
    other_site_positions <- as.integer(paste(sites_table[,4]))
    binding_sites <- compare_sites(binding_site_positions, start, end)
    active_sites <- compare_sites(active_site_positions, start, end)
    other_sites <- compare_sites(other_site_positions, start, end)
    all_sites <- cbind(binding_sites, active_sites, other_sites)
    output <- rbind(output, all_sites)
  }
  return(cbind(LFQ_table, output))
}

if (quantitation_level == "peptide") {
return(cbind(attach_enzyme_class(match_mod_peptides(LFQ_table))))
} else {return(attach_enzyme_class(LFQ_table))}
}
