#' Mean Plot of Maxquant Data
#'
#' This function generates a plot of the mean value of each protein for Group 1
#' (y) vs Group 2 (x) from a processed set of MaxQuant data. The data frame
#' should be the same as one that would be passed to `volcano_plot_maxquant()`
#' i.e. should contain the following columns: Protein.IDs, Gene.names,
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
#' @param xy_line boolean indicating whether to add a line to indicate x = y
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
#' @param label_points string indicating which points to label on the plot.
#'   Options are "none" (default) for no labels, "all" (attempt to label all
#'   points), a character vector of selected proteins to label (e.g. c("P12345",
#'   "Q67890", "IPO5") can be Gene.names, Protein.names, Protein.ID or a mix but
#'   the values will be used to label), or a code for auto-labelling (e.g.
#'   "3_all" or "5_10" which would label all significantly de-enriched points
#' (colour code 3) or the top 10 significantly enriched proteins (colour code
#' 5), a code of "top_20" will label the top 20 proteins by manhattan distance
#' regardless of colour code).
#' @param point_args a list of arguments to be passed to `geom_point()` for the
#'   points in the volcano plot (default is `list(size = 1, shape = 16, alpha =
#'   0.5)`). If the user specifies additional arguments or changes a default
#'   argument, the defaults will be merged with the user-specified entries.
#' @param label_args a list of arguments to be passed to
#'   `ggrepel::geom_text_repel()` for the labels in the volcano plot (default is
#'   `list(size = 1.5, max.overlaps = 25, segment.alpha = 0.5, segment.size =
#'   0.2, colour = "black")`). If the user specifies additional arguments or
#'   changes a default argument, the defaults will be merged with the
#'   user-specified entries.
#'
#' @returns ggplot object containing the volcano plot
#'
#' @import ggplot2
#' @import ggrepel
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
#' mean_plot_maxquant(df_subset)
#' }
mean_plot_maxquant <- function(df_subset = NULL,
                               meas = "LFQ.intensity",
                               threshold_p = 0.05,
                               threshold_fc = 1,
                               xy_line = FALSE,
                               vp_colours = NULL,
                               groups = NULL,
                               x_label = NULL,
                               y_label = NULL,
                               fsize = 8,
                               label_points = "none",
                               point_args = list(size = 1,
                                                 shape = 16,
                                                 alpha = 0.5),
                               label_args = list(size = 1.5,
                                                 max.overlaps = 25,
                                                 segment.alpha = 0.5,
                                                 segment.size = 0.2,
                                                 colour = "black")) {
  # satisfy R CMD check
  meas.ratio <- neg.log10.p.value <- vp_colorcode <- point_labels <- NULL
  meas.group1.mean <- meas.group2.mean <- NULL

  # if user has not supplied a data frame, stop and print an error message
  if (is.null(df_subset)) {
    stop("No data frame supplied. Please provide a data frame containing the subset of the MaxQuant data to be used for the mean plot.")
  }

  # check that meas ends in . and if not, add it
  if (!grepl("\\.$", meas)) {
    meas <- paste0(meas, ".")
  }

  if(!validate_df_subset(df_subset)) {
    stop("Data frame does not contain the required columns to generate the mean plot.")
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

  if(is.null(vp_colours)) {
    # integers 0 through 5 are possible
    vp_colours <- c("0" = "#a0a0a0", "1" = "#808080",
                    "2" = "#606060", "3" = "#8080ff",
                    "4" = "#606060", "5" = "#ff80ff")
  }

  # add a column for the point labels based on the label_points argument
  df_subset <- add_label_column(df_subset, label_points)

  # Merge user-supplied geom args with defaults so partial overrides keep
  # unspecified defaults.
  point_args <- utils::modifyList(
    list(size = 1, shape = 16, alpha = 0.5),
    point_args
  )
  label_args <- utils::modifyList(
    list(size = 1.5,
         max.overlaps = 25,
         segment.alpha = 0.5,
         segment.size = 0.2,
         colour = "black"),
    label_args
  )

  ## Generate the volcano plot

  p <- ggplot(df_subset, aes(x = meas.group2.mean, y = meas.group1.mean,
                             colour = vp_colorcode, label = point_labels))
  if(xy_line) {
    p <- p + geom_abline(slope = 1, intercept = 0,
                         linetype = "dashed", colour = "grey")
  }
  # adding the points to the plot
  p <- p + do.call(geom_point, point_args) +
    scale_colour_manual(values = vp_colours)
  # add labels if requested
  if(any(nzchar(df_subset$point_labels))) {
    p <- p + do.call(ggrepel::geom_text_repel, label_args)
  }
  # label the axes
  if(is.null(x_label)) {
    p <- p + xlab(group2)
  } else {
    p <- p + xlab(x_label)
  }
  if(is.null(y_label)) {
    p <- p + ylab(group1)
  } else {
    p <- p + ylab(y_label)
  }
  # theme for the plot
  p <- p +
    theme_classic(fsize) +
    theme(legend.position = "none")

  return (p)
}
