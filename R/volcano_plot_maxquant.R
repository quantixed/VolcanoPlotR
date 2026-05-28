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
#' @param p_line boolean indicating whether to add a horizontal line at the
#'   p-value threshold
#' @param zero_line boolean indicating whether to add a vertical line at x=0
#' @param x_line boolean indicating whether to add vertical lines at the fold
#'   change thresholds
#' @param vp_colours a named vector of colours for the volcano plot. The names
#'   of the vector should be the integers 0 through 5, which correspond to the
#'   different combinations of p-value and fold change thresholds. If NULL,
#'   default colours will be used.
#' @param groups a vector of two strings containing the names of the
#'   experimental (these are used for the x-axis label). The purpose is to allow
#'   the user to substitute a better label than that used in the MaxQuant data.
#'   If NULL the values from the data are deduced and used.
#' @param x_label string specifying the label for the x-axis. If NULL, a default
#'   label will be generated based on the group names.
#' @param y_label string specifying the label for the y-axis. If NULL, a default
#'   label will be generated.
#' @param fsize numeric indicating the font size to use for the plot (default is
#'   8)
#'
#' @returns ggplot object containing the volcano plot
#'
#' @import ggplot2
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # load the MaxQuant data
#' df <- load_maxquant()
#' # process the MaxQuant data to get the subset for the volcano plot
#' df_subset <- process_maxquant(df, group1 = "Treatment", group2  = "Control")
#' # generate the volcano plot
#' volcano_plot_maxquant(df_subset)
#' }
volcano_plot_maxquant <- function(df_subset = NULL,
                                  meas = "LFQ.intensity",
                                  threshold_p = 0.05,
                                  threshold_fc = 1,
                                  p_line = TRUE,
                                  zero_line = TRUE,
                                  x_line = FALSE,
                                  vp_colours = NULL,
                                  groups = NULL,
                                  x_label = NULL,
                                  y_label = NULL,
                                  fsize = 8) {
  # satisfy R CMD check
  meas.ratio <- neg.log10.p.value <- vp_colorcode <- vp_label <- NULL

  # if user has not supplied a data frame, stop and print an error message
  if (is.null(df_subset)) {
    stop("No data frame supplied. Please provide a data frame containing the subset of the MaxQuant data to be used for the volcano plot.")
  }

  # check that meas ends in . and if not, add it
  if (!grepl("\\.$", meas)) {
    meas <- paste0(meas, ".")
  }

  if(!validate_df_subset(df_subset)) {
    stop("Data frame does not contain the required columns to generate the volcano plot.")
  }

  if(is.null(groups)) {
    groups <- determine_groups_from_df_subset(df_subset, meas)

    group1 <- groups[1]
    group2 <- groups[2]
  } else {
    # user defined groups for x axis label
    group1 <- groups[1]
    group2 <- groups[2]
  }

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

  p <- ggplot()
  if(zero_line) {
    p <- p + geom_vline(xintercept = 0, linetype = "dashed", colour = "grey")
  }
  if(p_line) {
    p <- p + geom_hline(yintercept = -log10(threshold_p), linetype = "dashed", colour = "grey")
  }
  if(x_line) {
    p <- p + geom_vline(xintercept = c(-threshold_fc, threshold_fc), linetype = "dashed", colour = "grey")
  }
  # adding the points to the plot
  p <- p + geom_point(df_subset, aes(x = meas.ratio, y = neg.log10.p.value,
                              colour = vp_colorcode),
               size = 1, shape = 16, alpha = 0.5) +
    scale_colour_manual(values = vp_colours)
  # label the axes
  if(is.null(x_label)) {
    p <- p + xlab(bquote(.(group1) ~ - ~ .(group2) ~ (Log[2])))
  } else {
    p <- p + xlab(x_label)
  }
  if(is.null(y_label)) {
    p <- p + ylab(bquote(P~value~(-Log[10])))
  } else {
    p <- p + ylab(y_label)
  }
  # theme for the plot
  p <- p +
    theme_classic(fsize) +
    theme(legend.position = "none")

  return (p)
}




