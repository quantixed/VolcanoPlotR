#' PCA Plot of MaxQuant Sample Columns
#'
#' This function performs principal component analysis on the measurement
#' columns present in a processed MaxQuant subset and plots PC1 against PC2,
#' with each sample column shown as a point. The input data frame should be the
#' same processed subset that would be passed to `volcano_plot_maxquant()`.
#'
#' @param df_subset a data frame containing the subset of the MaxQuant data that
#'   will be used to generate the plot
#' @param meas string indicating the measurement type (e.g. "LFQ.intensity",
#'   "Intensity", etc.)
#' @param groups optional character vector used to relabel the experimental
#'   groups in the legend. If `NULL`, group names are deduced from the
#'   measurement column names.
#' @param by_protein boolean indicating whether to perform PCA on the protein
#'   rows (default is FALSE). If TRUE, PCA will be performed on the protein rows
#'   rather than the sample columns.
#' @param threshold_p numeric indicating the p-value threshold. Used if
#'   by_protein is TRUE to colour the points in the PCA plot based on
#'   significance.
#' @param threshold_fc numeric indicating the fold change threshold (in log2
#'   space, i.e. 1 is a 2-fold change, 2 is a 4-fold change, etc.). Used if
#'   by_protein is TRUE to colour the points in the PCA plot based on
#'   significance.
#' @param vp_colours a named vector of colours for the PCA plot.
#' @param x_label string specifying the label for the x-axis. If `NULL`, a
#'   default label showing the variance explained by PC1 will be generated.
#' @param y_label string specifying the label for the y-axis. If `NULL`, a
#'   default label showing the variance explained by PC2 will be generated.
#' @param xy_line boolean indicating whether to add a line to indicate x = y = 0
#' @param fsize numeric indicating the font size to use for the plot (default is
#'   8)
#' @param label_points boolean indicating whether to label sample points
#'   (default is NULL and triggers a setting of FALSE if `by_protein` is TRUE,
#'   and TRUE if `by_protein` is is FALSE). If TRUE, sample points will be
#'   labelled with their sample names. If FALSE, no labels will be added.
#' @param point_args a list of arguments to be passed to `geom_point()` for the
#'   points in the PCA plot (default is `list(size = 3, alpha = 0.8)`). If the
#'   user specifies additional arguments or changes a default argument, the
#'   defaults will be merged with the user-specified entries.
#' @param label_args a list of arguments to be passed to
#'   `ggrepel::geom_text_repel()` for the labels in the PCA plot (default is
#'   `list(size = 3, max.overlaps = Inf, segment.alpha = 0.5, segment.size =
#'   0.2, colour = "black")`). If the user specifies additional arguments or
#'   changes a default argument, the defaults will be merged with the
#'   user-specified entries.
#'
#' @returns ggplot object containing the PCA plot
#'
#' @import ggplot2
#' @import ggrepel
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # load the MaxQuant data
#' df <- load_maxquant(file = "proteinGroups.txt", datadir = "inst/extdata")
#' # process the MaxQuant data to get the subset for plotting
#' df_subset <- process_maxquant(df, group1 = "WT", group2 = "Control")
#' # generate the PCA plot
#' pca_plot_maxquant(df_subset)
#' }
pca_plot_maxquant <- function(df_subset = NULL,
                              meas = "LFQ.intensity",
                              groups = NULL,
                              by_protein = FALSE,
                              threshold_p = 0.05,
                              threshold_fc = 1,
                              vp_colours = NULL,
                              x_label = NULL,
                              y_label = NULL,
                              xy_line = TRUE,
                              fsize = 8,
                              label_points = NULL,
                              point_args = list(size = 3,
                                                alpha = 0.8),
                              label_args = list(size = 3,
                                                max.overlaps = Inf,
                                                segment.alpha = 0.5,
                                                segment.size = 0.2,
                                                colour = "black")) {
  PC1 <- PC2 <- group <- sample <- NULL

  if (is.null(df_subset)) {
    stop("No data frame supplied. Please provide a data frame containing the subset of the MaxQuant data to be used for the PCA plot.")
  }

  if (!grepl("\\.$", meas)) {
    meas <- paste0(meas, ".")
  }

  if (!validate_df_subset(df_subset)) {
    stop("Data frame does not contain the required columns to generate the PCA plot.")
  }

  meas_columns <- colnames(df_subset)[startsWith(colnames(df_subset), meas)]

  if (length(meas_columns) < 2) {
    stop("At least two measurement columns are required to generate the PCA plot.")
  }

  sample_matrix <- as.matrix(df_subset[, meas_columns, drop = FALSE])
  storage.mode(sample_matrix) <- "numeric"

  if (anyNA(sample_matrix) || any(!is.finite(sample_matrix))) {
    stop("Measurement columns contain missing or non-finite values; please process or clean the data before running PCA.")
  }

  if (nrow(sample_matrix) < 2) {
    stop("At least two proteins are required to generate the PCA plot.")
  }

  if(by_protein) {
    sample_matrix <- t(sample_matrix)
  }
  pca_fit <- stats::prcomp(t(sample_matrix), center = TRUE, scale. = TRUE)

  if (ncol(pca_fit$x) < 2) {
    stop("PCA did not produce two principal components to plot.")
  }

  observed_groups <- determine_groups_from_df_subset(df_subset, meas)
  sample_names <- substring(meas_columns, nchar(meas) + 1)
  sample_groups <- gsub("[0-9]+$", "", sample_names)
  sample_groups <- gsub("[._]$", "", sample_groups)

  if (!is.null(groups)) {
    if (length(groups) != length(observed_groups)) {
      stop("'groups' must contain the same number of entries as the groups present in the data.")
    }
    rename_map <- stats::setNames(groups, observed_groups)
    sample_groups <- unname(rename_map[sample_groups])
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


  if(by_protein) {
    sample_names <- df_subset$Gene.names
    sample_groups <- df_subset$vp_colorcode
  }


  # create data frame for ggplot
  plot_df <- data.frame(
    sample = sample_names,
    group = sample_groups,
    PC1 = pca_fit$x[, 1],
    PC2 = pca_fit$x[, 2],
    stringsAsFactors = FALSE
  )

  # process graphing arguments
  if(is.null(label_points)) {
    label_points <- !by_protein
  }

  point_args <- utils::modifyList(
    list(size = 3, alpha = 0.8),
    point_args
  )
  label_args <- utils::modifyList(
    list(size = 3,
         max.overlaps = Inf,
         segment.alpha = 0.5,
         segment.size = 0.2,
         colour = "black"),
    label_args
  )

  explained <- summary(pca_fit)$importance[2, 1:2] * 100

  p <- ggplot(plot_df, aes(x = PC1, y = PC2, colour = group, label = sample))
  if (xy_line) {
    p <- p + geom_hline(yintercept = 0, linetype = "dashed", colour = "grey") +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "grey")
  }
  p <- p + do.call(geom_point, point_args)
  if (by_protein) {
    p <- p + scale_colour_manual(values = vp_colours)
  } else {
    p <- p + scale_colour_discrete()
  }

  if (isTRUE(label_points)) {
    p <- p + do.call(ggrepel::geom_text_repel, label_args)
  }

  if (is.null(x_label)) {
    p <- p + xlab(sprintf("PC1 (%.1f%%)", explained[1]))
  } else {
    p <- p + xlab(x_label)
  }
  if (is.null(y_label)) {
    p <- p + ylab(sprintf("PC2 (%.1f%%)", explained[2]))
  } else {
    p <- p + ylab(y_label)
  }

  p <- p +
    theme_classic(fsize) +
    labs(colour = "Group") +
    theme(legend.position = "none")

  return(p)
}
