#' multi_volcano_plots: directly converts MaxQuant output to multiple volcano plots
#' @seealso \code{\link{pairwise_LFQ}}
#'          \code{\link{append_ec_sites}}
#'          \code{\link{plot_volcano}}
#' @param raw                 a dataframe by reading modificationSpecificPeptides.txt
#' @param metadata            a dataframe that maches the MaxQuant input. Column 1: Intensity (such as Intensity samplename, same as the column names in modificationSpecificPeptides.txt) name Column 2: Replicate group (use the same name for each group of replicates)
#' @param name_probe_mod      a string vector of chemical probe/modification names, such as c("Mod1", "Mod2"), must match MaxQuant input
#' @param max_each_mod        a integer as the maximal number of modifications on a single peptide, set for each chemical probe
#' @param max_total_mods      a integer as the maximal number of modifications on a single peptide, set for all chemical probes Note max_each_mod must not be less than max_total_mods
#' @param quantitation_level  a string, must be either "peptide" or "protein"
#' @param background_check    a boolean, FALSE = quantify probe-modified peptides, TRUE = quantify non-probe-modified peptides
#' @param normalize_to        a string, must be either "sum_all", "mean_all", (normalize to all peptides) "sum_background", or "mean_background" (normalize to background/non-probe-modified peptides).
#' @param xlim                a integer vector, such as c(-5, 5) for an x axis range of -5 to 5
#' @param ylim                a integer vector, such as c(0, 5) for an y axis range of 0 to 5
#' @param label_col_name      the input column name for labeling volcano plot data points such as "Gene.Names"
#' @param pCutoff             the p-Value cutoff, for instance, default p-value = 0.05
#' @param FCcutoff            the fold change cutoff, Note for ABPP, we are only interested in negative fold change (Lower intensity at higher inhibitor concentration)
#' @return  volcano plots
#' @examples  multi_volcano_plots(raw = raw, meta = meta, name_probe_mod = c("Mod"),
#'                    max_each_mod = 1, max_total_mods = 1, quantitation_level = "peptide" , background_check = FALSE, normalize_to = "mean_all",
#'                    xlim = c(-10, 3), ylim = c(0, 5), label_col_name = "Gene.Names", pCutoff = 0.05, FCcutoff = -2)
#' @export
multi_volcano_plots <- function(raw = read.delim("modificationSpecificPeptides.txt", header=TRUE, sep="\t"), metadata = read.delim("metadata.txt", header=TRUE, sep="\t"),
                                name_probe_mod = c("Mod"), max_each_mod = 1, max_total_mods = 1, quantitation_level = "peptide" , background_check = FALSE, normalize_to = NULL,
                                xlim = c(-10, 3), ylim = c(0, 5), label_col_name = "Gene.Names", pCutoff = 0.05, FCcutoff = -2) {

# Internal function generating multiple plots
multiplot <- function(plots, file, cols=2, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
 # plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


 log10p <- NULL
 log2fc <- NULL
 log10p <- -log10(pCutoff)
 log2fc <- -log2(abs(FCcutoff))
 LFQ_table <- pairwise_LFQ(raw = raw, meta = meta, name_probe_mod, max_each_mod = max_each_mod, max_total_mods = max_total_mods, quantitation_level = quantitation_level , background_check = background_check, normalize_to = normalize_to)
 LFQ_table_ec <- append_ec_sites(LFQ_table, quantitation_level = quantitation_level)
 inhibitors <- as.character(unique(gsub("[0-9]", "", meta$Replicate.group)))
 pvalue_index <- grep("p-value", colnames(LFQ_table_ec))
 fc_index <- grep("fold_change", colnames(LFQ_table_ec))
 plots <- NULL
 x <- NULL
 y <- NULL
 for (i in 1:length(inhibitors)) {
   pair_index <- grep(inhibitors[i], colnames(LFQ_table_ec))
   if (length(pair_index) > 0) {
     y <- colnames(LFQ_table_ec)[intersect(pair_index, pvalue_index)]
     x <- colnames(LFQ_table_ec)[intersect(pair_index, fc_index)]
     plots[[i]] <- plot_volcano(LFQ_table_ec, x, y, xlim = xlim, ylim = ylim, label_col_name = label_col_name, pCutoff = log10p, FCcutoff = log2fc, title = paste0(inhibitors[i], "/", name_probe_mod))
   }
 }
 plots[sapply(plots, is.null)] <- NULL
 return (multiplot(plots, cols = 2))
}


