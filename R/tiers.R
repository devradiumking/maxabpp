

#'@title Function make_tiers
#'@description Generate a tier list from a set of lists of protein groups identified by MaxQuant for cross-list comparisons
#'@export
make_tiers <- function(setLists) {
 #Call Max_Venn on setList
 maxVennOutput <- Max_Venn(setLists, IndividualAnalysis = TRUE)
 #Remove the null debugging field
 IntersectionSets <- maxVennOutput@IntersectionSets[-1]
 #Get len1 as the number of venn fields
 len1 <- length(IntersectionSets)
 n_levels <- sapply(1:len1,function(x)NULL)
 #Get len2 as the number of levels/tiers
 for (n_field in 1:len1)
 {
   field_name <- names(IntersectionSets)[n_field]
   n_levels[[n_field]] <- sum(as.integer(unlist(strsplit(field_name,""))))
 }
 len2 <- max(unlist(n_levels))
 levels <- sapply(1:len2,function(x)NULL)
 for (n_field in 1:len1)
 {
   field_name <- names(IntersectionSets)[n_field]
   field_level <- sum(as.integer(unlist(strsplit(field_name,""))))
   levels[[field_level]] <- union(levels[[field_level]], V3tags_new@IntersectionSets[[field_name]])
 }
 #The smaller the tier number is the higher the level, best tier has highest level
 n_tiers <- length(levels)
 tiers <- sapply(1:n_tiers,function(x)NULL)
 j <- n_tiers
 for (i in 1 : n_tiers) {
   tier_name <- paste("Tier ", as.roman(j), sep = "")
   names(tiers)[j] <- tier_name
   tiers[[j]] <- levels[[i]]
   j = j - 1
 }
 return(tiers)
}
