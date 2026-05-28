#' Volcano Plot of Maxquant Data
#'
#' This function generates a volcano plot from a subset of MaxQuant data. The
#' data frame should contain the following columns: Protein.IDs, Gene.names,
#' Protein.names, meas.ratio, p.value, and neg.log10.p.value. The function will
#' add a colour code column to the data frame based on the p-value and fold
#' change thresholds, and then use ggplot2 to create the volcano plot.
#'
#' @param df_subset a data frame containing the subset of the MaxQuant data that
#'   will be used to generate the plot
#' @param meas string indicating the measurement type (e.g. "LFQ.intensity",
#'   "Intensity", etc.)
#' @param threshold_p numeric indicating the p-value threshold
#' @param threshold_fc numeric indicating the fold change threshold (in log2
#'   space, i.e. 1 is a 2-fold change, 2 is a 4-fold change, etc.)
#' @param vp_colours a named vector of colours for the volcano plot. The names
#'   of the vector should be the integers 0 through 5, which correspond to the
#'   different combinations of p-value and fold change thresholds. If NULL,
#'   default colours will be used.
#' @param fsize numeric indicating the font size
#'
#' @returns ggplot object containing the volcano plot
#'
#' @import ggplot2
#'
#' @export
volcano_plot_maxquant <- function(df_subset = NULL, meas = "LFQ.intensity",
                                  threshold_p = 0.05, threshold_fc = 1,
                                  vp_colours = NULL, fsize = 8) {
  # satisfy R CMD check
  meas.ratio <- neg.log10.p.value <- vp_colorcode <- vp_label <- NULL

  # check that meas ends in . and if not, add it
  if (!grepl("\\.$", meas)) {
    meas <- paste0(meas, ".")
  }

  if(!validate_df_subset(df_subset)) {
    stop("Data frame does not contain the required columns to generate the volcano plot.")
  }

  groups <- determine_groups_from_df_subset(df_subset, meas)

  group1 <- groups[1]
  group2 <- groups[2]

  # add a colour code column for the volcano plot we'll do it using bit setting:
  # bit 0 is p-value < threshold_p, bit 1 is meas.ratio <= -1, bit 2 is
  # meas.ratio >= 1 then change to character for the purposes of applying
  # colours to points
  df_subset$vp_colorcode <- ifelse(df_subset$p.value < threshold_p, 1, 0) +
    ifelse(df_subset$meas.ratio <= -threshold_fc, 2, 0) +
    ifelse(df_subset$meas.ratio >= threshold_fc, 4, 0)
  df_subset$vp_colorcode <- as.character(df_subset$vp_colorcode)

  if(is.null(vp_label)) {
    # integers 0 through 5 are possible
    vp_colours <- c("0" = "#a0a0a0", "1" = "#808080",
                    "2" = "#606060", "3" = "#8080ff",
                    "4" = "#606060", "5" = "#ff80ff")
  }

  p <- ggplot(df_subset, aes(x = meas.ratio, y = neg.log10.p.value)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey") +
    geom_hline(yintercept = -log10(threshold_p), linetype = "dashed", colour = "grey") +
    geom_point(aes(colour = vp_colorcode),
               size = 1, shape = 16, alpha = 0.5) +
    scale_colour_manual(values = vp_colours) +
    xlab(bquote(.(group1) ~ - ~ .(group2) ~ (Log[2]))) +
    ylab(bquote(P~value~(-Log[10]))) +
    theme_classic(fsize) +
    theme(legend.position = "none")

  return (p)
}




