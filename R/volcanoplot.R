#' Generates a volcano plot from the output of append_ec_sites()
#' @seealso \code{\link{pairwise_LFQ}}
#'          \code{\link{append_ec_sites}}
#' @param LFQ_table_ec        a dataframe saved as the append_ec_sites() output
#' @param x                   log2FC column name
#' @param y                   -log10pvalue column name
#' @param xlim                a integer vector, such as c(-5, 5) for an x axis range of -5 to 5
#' @param ylim                a integer vector, such as c(0, 5) for an y axis range of 0 to 5
#' @param label_col_name      the input column name for labeling volcano plot data points such as "Gene.Names"
#' @param pCutoff             the p-Value cutoff on -log10 scale, for instance, use 1.3 for p-value = 0.05
#' @param FCcutoff            the fold change cutoff on log2 scale, for instance, use -1 for -2 fold change, Note for ABPP, we are only interested in negative fold change
#' @param title               a string as the title of the volcano plot
#' @return                    a volcano plot
#' @examples plot_volcano(output2, x, y, xlim, ylim, label_col_name, pCutoff, FCcutoff, title)
#' @export
plot_volcano <- function(LFQ_table_ec, x, y, xlim, ylim, label_col_name, pCutoff, FCcutoff, title) {
  #rownames(LFQ_table_ec) <- LFQ_table_ec$protein$protein
  #rownames(LFQ_table_ec) #<- paste(LFQ_table_ec$protein$protein,":",LFQ_table_ec$`protein names`$`protein names`)
  # create custom key-value pairs for "Oxidoreductases" "Hydrolases" "Ligases" "Transferases" "Isomerases" "Lyases" "Translocases" "Other proteins"
for (nrow_LFQ_table in 1:nrow(LFQ_table_ec)) {
  if (LFQ_table_ec[[label_col_name]][nrow_LFQ_table] == "")
    LFQ_table_ec[[label_col_name]][nrow_LFQ_table] <- LFQ_table_ec[["Proteins"]][nrow_LFQ_table]
}
  major_label <- function (dataframe, label_col_name) {
    labs <- dataframe[[label_col_name]]
    active_sites <- dataframe[["active_sites"]]
    binding_sites <- dataframe[["binding_sites"]]
    other_sites <- dataframe[["other_sites"]]
    lab <- NULL
    if (length(labs) < 1) {
      return(NULL)
    } else {
    for (n_labs in 1:length(labs)) {
      lab[n_labs] <- str_split(labs, ";")[[n_labs]][1]
      active_site_hit <- str_extract(active_sites[[n_labs]], "[0-9]")
      binding_site_hit <- str_extract(binding_sites[[n_labs]], "[0-9]")
      other_site_hit <- str_extract(other_sites[[n_labs]], "[0-9]")
      if (!is.na(active_site_hit)) {lab[n_labs] <- paste0(lab[n_labs], "(Active Site)")}
      if (!is.na(binding_site_hit)) {lab[n_labs] <- paste0(lab[n_labs], "(Binding Site)")}
      if (!is.na(other_site_hit)) {lab[n_labs] <- paste0(lab[n_labs], "(Other Site)")}
    }
    return(lab)
    }
  }

  #Calculate ID summary
  identified <- nrow(LFQ_table_ec)
  isfinite_x <- LFQ_table_ec[[x]]
  finite_x <- LFQ_table_ec[is.finite(isfinite_x) , ]
  overinhibited_x <- isfinite_x == -Inf & !is.na(isfinite_x)
  infinite_x <- LFQ_table_ec[overinhibited_x , ]
  na_y <- finite_x[[y]]
  nna_y <- subset(finite_x, is.na(na_y) == FALSE)
  criteria_x <- nna_y[[x]]
  sub1 <- finite_x[criteria_x <= FCcutoff , ]
  criteria_y <- sub1[[y]]
  sub2 <- sub1[criteria_y >= pCutoff , ]
  criteria_y2 <- infinite_x[[y]]
  overinhibited <- infinite_x[criteria_y2 >= pCutoff , ]
  num_overinhibited <- length(na.exclude(major_label(overinhibited, label_col_name)))
  quantified <- nrow(nna_y) + num_overinhibited
  significant <- length(na.exclude(major_label(sub2, label_col_name))) + num_overinhibited
  n_active <- as.character(summary(str_detect(LFQ_table_ec$active_sites, "[0-9]"))[3])
  n_binding <- as.character(summary(str_detect(LFQ_table_ec$binding_sites, "[0-9]"))[3])
  n_others <- as.character(summary(str_detect(LFQ_table_ec$other_sites, "[0-9]"))[3])
  # set the base colour as 'black'
  keyvals <- rep('black', nrow(finite_x))

  names(keyvals) <- rep('Other proteins', nrow(finite_x))

  keyvals[which(finite_x$ec_list == "Oxidoreductase")] <- 'firebrick'
  names(keyvals)[(finite_x$ec_list == "Oxidoreductase")] <- 'Oxidoreductase'

  keyvals[which(finite_x$ec_list == "Transferase")] <- 'chocolate'
  names(keyvals)[(finite_x$ec_list == "Transferase")] <- 'Transferase'

  keyvals[which(finite_x$ec_list == "Hydrolase")] <- 'darkturquoise'
  names(keyvals)[(finite_x$ec_list == "Hydrolase")] <- 'Hydrolase'

  keyvals[which(finite_x$ec_list == "Lyase")] <- 'yellow3'
  names(keyvals)[(finite_x$ec_list == "Lyase")] <- 'Lyase'

  keyvals[which(finite_x$ec_list == "Isomerase")] <- 'deeppink2'
  names(keyvals)[(finite_x$ec_list == "Isomerase")] <- 'Isomerase'

  keyvals[which(finite_x$ec_list == "Ligase")] <- 'limegreen'
  names(keyvals)[(finite_x$ec_list == "Ligase")] <- 'Ligase'

  keyvals[which(finite_x$ec_list == "Translocase")] <- 'purple3'
  names(keyvals)[(finite_x$ec_list == "Translocase")] <- 'Translocase'

  EnhancedVolcano <- function (toptable, lab, x, y, selectLab = NULL, xlim = c(min(toptable[,
                                                                                            x], na.rm = TRUE), max(toptable[, x], na.rm = TRUE)), ylim = c(0,
                                                                                                                                                           max(toptable[, y], na.rm = TRUE) + 5), xlab = bquote(~Log[2] ~
                                                                                                                                                                                                                  "FC"), ylab = bquote(~-Log[10] ~ italic(P) ~ "-value"), axisLabSize = 18,
                               title = "Volcano plot", subtitle = NULL,
                               caption = paste0("Total = ", nrow(toptable), " Identified"),
                               titleLabSize = 18, subtitleLabSize = 14, captionLabSize = 14,
                               pCutoff = 1.3, pLabellingCutoff = pCutoff, FCcutoff = -1,
                               cutoffLineType = "longdash", cutoffLineCol = "black", cutoffLineWidth = 0.4,
                               transcriptPointSize = 0.8, transcriptLabSize = 3, transcriptLabCol = "black",
                               transcriptLabFace = "plain", transcriptLabhjust = 0, transcriptLabvjust = 1.5,
                               boxedlabels = FALSE, shape = 19, shapeCustom = NULL, col = c("grey30",
                                                                                            "forestgreen", "royalblue", "red2"), colCustom = NULL,
                               colAlpha = 1/2, legend = c("NS", "Log2 FC", "P", "P & Log2 FC"),
                               legendPosition = "top", legendLabSize = 14, legendIconSize = 4,
                               legendVisible = TRUE, shade = NULL, shadeLabel = NULL, shadeAlpha = 1/2,
                               shadeFill = "grey", shadeSize = 0.01, shadeBins = 2, drawConnectors = FALSE,
                               widthConnectors = 0.5, typeConnectors = "closed", endsConnectors = "first",
                               lengthConnectors = unit(0.01, "npc"), colConnectors = "grey10",
                               hline = NULL, hlineType = "longdash", hlineCol = "black",
                               hlineWidth = 0.4, vline = NULL, vlineType = "longdash", vlineCol = "black",
                               vlineWidth = 0.4, gridlines.major = TRUE, gridlines.minor = TRUE,
                               border = "partial", borderWidth = 0.8, borderColour = "black")
  {
    if (!requireNamespace("ggplot2")) {
      stop("Please install ggplot2 first.", call. = FALSE)
    }
    if (!requireNamespace("ggrepel")) {
      stop("Please install ggrepel first.", call. = FALSE)
    }
    if (!is.numeric(toptable[, x])) {
      stop(paste(x, " is not numeric!", sep = ""))
    }
    if (!is.numeric(toptable[, y])) {
      stop(paste(y, " is not numeric!", sep = ""))
    }

    finite_toptable <- toptable[, x] > -1e100 & toptable[, x] < 1e100
    toptable <- subset(toptable, finite_toptable)

    i <- xvals <- yvals <- Sig <- NULL
    toptable <- as.data.frame(toptable)
    toptable$Sig[(toptable[, x] > -1e100)] <- "NS"
    toptable$Sig[(toptable[, x] < FCcutoff & toptable[, x] > -1e100)] <- "FC"
    toptable$Sig[(toptable[, y] < pCutoff & toptable[, x] > -1e100)] <- "P"
    toptable$Sig[(toptable[, y] < pCutoff) & (toptable[, x] > FCcutoff & toptable[, x] > -1e100)] <- "FC_P"
    toptable$Sig <- factor(toptable$Sig, levels = c("NS", "FC",
                                                    "P", "FC_P"))
    if (min(toptable[, y], na.rm = TRUE) == 0) {
      warning(paste("One or more P values is 0.", "Converting to minimum possible value..."),
              call. = FALSE)
      toptable[which(toptable[, y] == 0), y] <- .Machine$double.xmin
    }

    toptable$lab <- major_label(toptable, lab)
    toptable$xvals <- toptable[, x]
    toptable$yvals <- toptable[, y]
    if (!is.null(selectLab)) {
      names.new <- rep(NA, length(toptable$lab))
      indices <- which(toptable$lab %in% selectLab)
      names.new[indices] <- toptable$lab[indices]
      toptable$lab <- names.new
    }
    th <- theme_bw(base_size = 24) + theme(legend.background = element_rect(),
                                           plot.title = element_text(angle = 0, size = titleLabSize,
                                                                     face = "bold", vjust = 1), plot.subtitle = element_text(angle = 0,
                                                                                                                             size = subtitleLabSize, face = "plain", vjust = 1),
                                           plot.caption = element_text(angle = 0, size = captionLabSize,
                                                                       face = "plain", vjust = 1), axis.text.x = element_text(angle = 0,
                                                                                                                              size = axisLabSize, vjust = 1), axis.text.y = element_text(angle = 0,
                                                                                                                                                                                         size = axisLabSize, vjust = 1), axis.title = element_text(size = axisLabSize),
                                           legend.position = legendPosition, legend.key = element_blank(),
                                           legend.key.size = unit(0.5, "cm"), legend.text = element_text(size = legendLabSize),
                                           title = element_text(size = legendLabSize), legend.title = element_blank())
    if (!is.null(colCustom) & !is.null(shapeCustom)) {
      plot <- ggplot(toptable, aes(x = xvals, D)) +
        th + guides(colour = guide_legend(order = 1, override.aes = list(size = legendIconSize)),
                    shape = guide_legend(order = 2, override.aes = list(size = legendIconSize))) +
        geom_point(aes(color = factor(names(colCustom)),
                       shape = factor(names(shapeCustom))), alpha = colAlpha,
                   size = transcriptPointSize, na.rm = TRUE) + scale_color_manual(values = colCustom) +
        scale_shape_manual(values = shapeCustom)
    }
    else if (!is.null(colCustom) & is.null(shapeCustom) & length(shape) ==
             1) {
      plot <- ggplot(toptable, aes(x = xvals, y = yvals)) +
        th + guides(colour = guide_legend(order = 1, override.aes = list(size = legendIconSize)),
                    shape = guide_legend(order = 2, override.aes = list(size = legendIconSize))) +
        geom_point(aes(color = factor(names(colCustom))),
                   alpha = colAlpha, shape = shape, size = transcriptPointSize,
                   na.rm = TRUE) + scale_color_manual(values = colCustom) +
        scale_shape_manual(guide = TRUE)
    }
    else if (!is.null(colCustom) & is.null(shapeCustom) & length(shape) ==
             4) {
      plot <- ggplot(toptable, aes(x = xvals, y = yvals)) +
        th + guides(colour = guide_legend(order = 1, override.aes = list(size = legendIconSize)),
                    shape = guide_legend(order = 2, override.aes = list(size = legendIconSize))) +
        geom_point(aes(color = factor(names(colCustom)),
                       shape = factor(Sig)), alpha = colAlpha, size = transcriptPointSize,
                   na.rm = TRUE) + scale_color_manual(values = colCustom) +
        scale_shape_manual(values = c(NS = shape[1], FC = shape[2],
                                      P = shape[3], FC_P = shape[4]), labels = c(NS = legend[1],
                                                                                 FC = paste(legend[2], sep = ""), P = paste(legend[3],
                                                                                                                            sep = ""), FC_P = paste(legend[4], sep = "")),
                           guide = TRUE)
    }
    else if (is.null(colCustom) & !is.null(shapeCustom)) {
      plot <- ggplot(toptable, aes(x = xvals, y = yvals)) +
        th + guides(colour = guide_legend(order = 1, override.aes = list(size = legendIconSize)),
                    shape = guide_legend(order = 2, override.aes = list(size = legendIconSize))) +
        geom_point(aes(color = factor(Sig), shape = factor(names(shapeCustom))),
                   alpha = colAlpha, size = transcriptPointSize,
                   na.rm = TRUE) + scale_color_manual(values = c(NS = col[1],
                                                                 FC = col[2], P = col[3], FC_P = col[4]), labels = c(NS = legend[1],
                                                                                                                     FC = paste(legend[2], sep = ""), P = paste(legend[3],
                                                                                                                                                                sep = ""), FC_P = paste(legend[4], sep = ""))) +
        scale_shape_manual(values = shapeCustom)
    }
    else if (is.null(colCustom) & is.null(shapeCustom) & length(shape) ==
             1) {
      plot <- ggplot(toptable, aes(x = xvals, y = yvals)) +
        th + guides(colour = guide_legend(order = 1, override.aes = list(shape = shape,
                                                                         size = legendIconSize))) + geom_point(aes(color = factor(Sig)),
                                                                                                               alpha = colAlpha, shape = shape, size = transcriptPointSize,
                                                                                                               na.rm = TRUE, show.legend = legendVisible) + scale_color_manual(values = c(NS = col[1],
                                                                                                                                                                                          FC = col[2], P = col[3], FC_P = col[4]), labels = c(NS = legend[1],
                                                                                                                                                                                                                                              FC = paste(legend[2], sep = ""), P = paste(legend[3],
                                                                                                                                                                                                                                                                                         sep = ""), FC_P = paste(legend[4], sep = "")))
    }
    else if (is.null(colCustom) & is.null(shapeCustom) & length(shape) ==
             4) {
      plot <- ggplot(toptable, aes(x = xvals, y = yvals)) +
        th + guides(colour = guide_legend(order = 1, override.aes = list(shape = c(NS = shape[1],
                                                                                   FC = shape[2], P = shape[3], FC_P = shape[4]), size = legendIconSize))) +
        geom_point(aes(color = factor(Sig), shape = factor(Sig)),
                   alpha = colAlpha, size = transcriptPointSize,
                   na.rm = TRUE, show.legend = legendVisible) +
        scale_color_manual(values = c(NS = col[1], FC = col[2],
                                      P = col[3], FC_P = col[4]), labels = c(NS = legend[1],
                                                                             FC = paste(legend[2], sep = ""), P = paste(legend[3],
                                                                                                                        sep = ""), FC_P = paste(legend[4], sep = ""))) +
        scale_shape_manual(values = c(NS = shape[1], FC = shape[2],
                                      P = shape[3], FC_P = shape[4]), guide = FALSE)
    }
    plot <- plot + xlab(xlab) + ylab(ylab) + xlim(xlim[1], xlim[2]) +
      ylim(ylim[1], ylim[2]) + geom_vline(xintercept = c(FCcutoff
      ), linetype = cutoffLineType, colour = cutoffLineCol,
      size = cutoffLineWidth) + geom_hline(yintercept = pCutoff,
                                           linetype = cutoffLineType, colour = cutoffLineCol, size = cutoffLineWidth)
    plot <- plot + labs(title = title, subtitle = subtitle, caption = caption)
    if (!is.null(vline)) {
      plot <- plot + geom_vline(xintercept = vline, linetype = vlineType,
                                colour = vlineCol, size = vlineWidth)
    }
    if (!is.null(hline)) {
      plot <- plot + geom_hline(yintercept = -log10(hline),
                                linetype = hlineType, colour = hlineCol, size = hlineWidth)
    }
    if (border == "full") {
      plot <- plot + theme(panel.border = element_rect(colour = borderColour,
                                                       fill = NA, size = borderWidth))
    }
    else if (border == "partial") {
      plot <- plot + theme(axis.line = element_line(size = borderWidth,
                                                    colour = borderColour), panel.border = element_blank(),
                           panel.background = element_blank())
    }
    else {
      stop("Unrecognised value passed to 'border'. Must be 'full' or 'partial'")
    }
    if (gridlines.major == TRUE) {
      plot <- plot + theme(panel.grid.major = element_line())
    }
    else {
      plot <- plot + theme(panel.grid.major = element_blank())
    }
    if (gridlines.minor == TRUE) {
      plot <- plot + theme(panel.grid.minor = element_line())
    }
    else {
      plot <- plot + theme(panel.grid.minor = element_blank())
    }
    if (boxedlabels == FALSE) {
      if (drawConnectors == TRUE && is.null(selectLab)) {
        plot <- plot + geom_text_repel(data = subset(toptable, toptable[, y] > pLabellingCutoff & toptable[, x] < FCcutoff & toptable[, x] > -1e100), aes(label = subset(toptable, toptable[, y] > pLabellingCutoff & toptable[, x] < FCcutoff & toptable[, x] > -1e100)[, "lab"]), size = transcriptLabSize,
                                       segment.color = colConnectors, segment.size = widthConnectors,
                                       arrow = arrow(length = lengthConnectors, type = typeConnectors,
                                                     ends = endsConnectors), hjust = transcriptLabhjust,
                                       vjust = transcriptLabvjust, colour = transcriptLabCol,
                                       fontface = transcriptLabFace, na.rm = TRUE)
      }
      else if (drawConnectors == TRUE && !is.null(selectLab)) {
        plot <- plot + geom_text_repel(data = subset(toptable, !is.na(toptable[, "lab"])), aes(label = subset(toptable, !is.na(toptable[, "lab"]))[, "lab"]), size = transcriptLabSize,
                                       segment.color = colConnectors, segment.size = widthConnectors,
                                       arrow = arrow(length = lengthConnectors, type = typeConnectors,
                                                     ends = endsConnectors), hjust = transcriptLabhjust,
                                       vjust = transcriptLabvjust, colour = transcriptLabCol,
                                       fontface = transcriptLabFace, na.rm = TRUE)
      }
      else if (drawConnectors == FALSE && !is.null(selectLab)) {
        plot <- plot + geom_text(data = subset(toptable,  !is.na(toptable[, "lab"])), aes(label = subset(toptable, !is.na(toptable[, "lab"]))[, "lab"]), size = transcriptLabSize,
                                 check_overlap = TRUE, hjust = transcriptLabhjust,
                                 vjust = transcriptLabvjust, colour = transcriptLabCol,
                                 fontface = transcriptLabFace, na.rm = TRUE)
      }
      else if (drawConnectors == FALSE && is.null(selectLab)) {
        plot <- plot + geom_text(data = subset(toptable, toptable[, y] > pLabellingCutoff & toptable[, x] < FCcutoff & toptable[, x] > -1e100), aes(label = subset(toptable, toptable[, y] > pLabellingCutoff & toptable[, x] < FCcutoff & toptable[, x] > -1e100)[, "lab"]), size = transcriptLabSize,
                                 check_overlap = TRUE, hjust = transcriptLabhjust,
                                 vjust = transcriptLabvjust, colour = transcriptLabCol,
                                 fontface = transcriptLabFace, na.rm = TRUE)
      }
    }
    else {
      if (drawConnectors == TRUE && is.null(selectLab)) {
        plot <- plot + geom_label_repel(data = subset(toptable, toptable[, y] > pLabellingCutoff & toptable[, x] < FCcutoff & toptable[, x] > -1e100), aes(label = subset(toptable, toptable[, y] > pLabellingCutoff & toptable[, x] < FCcutoff)[, "lab"] & toptable[, x] > -1e100), size = transcriptLabSize,
                                        arrow = arrow(length = lengthConnectors, type = typeConnectors,
                                                      ends = endsConnectors), hjust = transcriptLabhjust,
                                        vjust = transcriptLabvjust, colour = transcriptLabCol,
                                        fontface = transcriptLabFace, na.rm = TRUE)
      }
      else if (drawConnectors == TRUE && !is.null(selectLab)) {
        plot <- plot + geom_label_repel(data = subset(toptable, !is.na(toptable[, "lab"])), aes(label = subset(toptable,  !is.na(toptable[, "lab"]))[, "lab"]), size = transcriptLabSize,
                                        segment.color = colConnectors, segment.size = widthConnectors,
                                        arrow = arrow(length = lengthConnectors, type = typeConnectors,
                                                      ends = endsConnectors), hjust = transcriptLabhjust,
                                        vjust = transcriptLabvjust, colour = transcriptLabCol,
                                        fontface = transcriptLabFace, na.rm = TRUE)
      }
      else if (drawConnectors == FALSE && !is.null(selectLab)) {
        plot <- plot + geom_label(data = subset(toptable, !is.na(toptable[, "lab"])), aes(label = subset(toptable, !is.na(toptable[, "lab"]))[, "lab"]), size = transcriptLabSize,
                                  hjust = transcriptLabhjust, vjust = transcriptLabvjust,
                                  colour = transcriptLabCol, fontface = transcriptLabFace,
                                  na.rm = TRUE)
      }
      else if (drawConnectors == FALSE && is.null(selectLab)) {
        plot <- plot + geom_label(data = subset(toptable, toptable[, y] > pLabellingCutoff & toptable[, x] < FCcutoff & toptable[, x] > -1e100), aes(label = subset(toptable, toptable[, y] > pLabellingCutoff & toptable[, x] < FCcutoff)[, "lab"] & toptable[, x] > -1e100), size = transcriptLabSize,
                                  hjust = transcriptLabhjust, vjust = transcriptLabvjust,
                                  colour = transcriptLabCol, fontface = transcriptLabFace,
                                  na.rm = TRUE)
      }
    }
    if (!is.null(shade)) {
      plot <- plot + stat_density2d(data = subset(toptable, rownames(toptable) %in% shade), fill = shadeFill,
                                    alpha = shadeAlpha, geom = "polygon", contour = TRUE,
                                    size = shadeSize, bins = shadeBins, show.legend = FALSE,
                                    na.rm = TRUE) + scale_fill_identity(name = shadeLabel,
                                                                        labels = shadeLabel, guide = "legend")
    }
    return(plot)
  }


plot <- EnhancedVolcano(LFQ_table_ec,
                  lab = label_col_name,
                  x = x,
                  y = y,
                  xlim = xlim,
                  ylim = ylim,
                  selectLab = NULL,
                  caption = paste0("Summary: ", identified, " Identified, ", quantified, " Quantified, ", num_overinhibited, " Overinhibited, ", significant, " Significant", "\n" , "Probe-modified peptides cover ", n_active, " Active Sites, ", n_binding, " Binding Sites, ", n_others, " Other Sites."),
                  colCustom = keyvals,
                  title = title,
                  pCutoff = pCutoff,
                  FCcutoff = FCcutoff,
                  transcriptPointSize = 1.5,
                  transcriptLabSize = 3,
                  col=c('black',"black","black" ,'red3'),
                  colAlpha = 1,
                  legendVisible = FALSE,
                  drawConnectors = TRUE,
                  widthConnectors = 0.5,
                  colConnectors = 'darkgray')



return(plot)
}

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
#' @examples  multi_volcano_plots(raw = raw, metadata = metadata, name_probe_mod = c("Mod"),
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
  LFQ_table <- pairwise_LFQ(raw = raw, metadata = metadata, name_probe_mod, max_each_mod = max_each_mod, max_total_mods = max_total_mods, quantitation_level = quantitation_level , background_check = background_check, normalize_to = normalize_to)
  LFQ_table_ec <- append_ec_sites(LFQ_table, quantitation_level = quantitation_level)
  inhibitors <- as.character(unique(gsub("[0-9]", "", metadata$Replicate.group)))
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
