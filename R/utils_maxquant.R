#' Validate a data frame to generate a volcano plot
#'
#' Checks if we have all the necessary components to generate a volcano plot.
#'
#' @param df_subset data frame containing the subset of the MaxQuant data that
#'   will be used to generate the volcano plot. This data frame should contain
#'   the columns "Protein.IDs", "Gene.names", "Protein.names", "meas.ratio",
#'   "p.value", and "neg.log10.p.value". The "Gene.names" and "Protein.names"
#'   columns are optional, but if they are not present, a warning will be
#'   printed.
#'
#' @returns A boolean indicating whether the data frame contains the required
#'   columns to generate a volcano plot, and a warning if the optional columns
#'   are not present.
#' @keywords internal
validate_df_subset <- function(df_subset) {
  required_columns <- c("Protein.IDs", "Gene.names", "Protein.names",
                        "meas.ratio", "p.value", "neg.log10.p.value")
  optional_columns <- c("Gene.names", "Protein.names")

  # check if df_subset has the required columns to make the volcano plot
  if (!all(required_columns %in% colnames(df_subset))) {
    return(FALSE)
  }

  # print a warning if the optional columns are not present, but don't stop the
  # function
  if (!all(optional_columns %in% colnames(df_subset))) {
    warning("Data frame is missing optional columns: Gene.names and/or Protein.names.")
  }

  return(TRUE)
}


#' Determine the experimental groups in a data frame subset
#'
#' This function takes a data frame containing a subset of the MaxQuant data and
#' determines the experimental groups that are present in the data frame by looking
#' at the column names.
#'
#' @param df_subset a data frame containing the subset of the MaxQuant data that will be used to
#' generate the volcano plot.
#' @param meas string indicating the measurement type
#'
#' @returns vector of two strings containing the names of the experimental groups
#' @keywords internal
determine_groups_from_df_subset <- function(df_subset, meas) {
  # find the experimental groups in the data frame using the names of the columns
  groups <- unique(gsub(meas, "", grep(meas, colnames(df_subset), value = TRUE)))
  # if no groups are found, stop and print an error message
  if (length(groups) == 0) {
    return(FALSE)
  } else {
    # these groups should be foo.1, foo.2, and bar.baz.1, bar.baz.2 so get the
    # group names by removing the .1 and .2 (or any integer) from the group
    # names
    groups <- unique(gsub("[0-9]+$", "", groups))
    # if there are trailing "." or "_" in the group names, remove them
    groups <- gsub("[._]$", "", groups)
  }

  return(groups)
}

#' Add a label column to the data frame subset for the volcano plot
#'
#' @param df_subset the data frame to be used for the volcano plot, which should
#'   contain the columns "Gene.names", "Protein.names", and "Protein.IDs" and a
#'   column "vp_colorcode" which contains the colour code for each point in the
#'   volcano plot (0, 1, 2, 3, 4, or 5)
#' @param label_points character string specifying which points to label in the
#'   volcano plot. This can be one of the following: "none", "all", a character
#'   vector of selected proteins to label, a code for auto-labelling
#'
#' @returns the input data frame with an additional column "point_labels" containing the
#'  labels for the points to be labelled in the volcano plot. The "point_labels"
#'
#' @keywords internal

add_label_column <- function(df_subset, label_points) {
  label_method <- label_count <- NULL
  n_label <- length(label_points)
  if(n_label > 1) {
    label_method <- "multiple"
  } else if(n_label == 1 && grepl("_", label_points)) {
    # label_points is a code, so split it into the method and the number of points to label
    label_method <- unlist(strsplit(label_points, "_"))[1]
    label_count <- suppressWarnings(as.numeric(unlist(strsplit(label_points, "_"))[2]))
  } else if(n_label == 1 && label_points %in% c("none", "all")) {
    label_method <- "allornone"
  } else {
    label_method <- "single"
  }

  if(is.null(label_method)) {
    warning("Invalid label method specified in label_points. No points will be labelled.")
    label_points <- "none"
    label_method <- "allornone"
  }

  # sort the data frame by manhattan distance in decreasing order
  output_df <- df_subset[order(df_subset$manhattan.distance,
                               decreasing = TRUE), ]

  # here we use the coded label method to determine which points to label.
  if(label_method == "top") {
    # select the first n rows of output_df where n is the label count
    last_row <- min(label_count, NROW(output_df))
    top_points <- output_df[1:last_row, ]
    label_points <- top_points$Gene.names
  } else if(label_method %in% c("0", "1", "2", "3", "4", "5")) {
    # subset the points with the specified colour code
    output_df <- output_df[output_df$vp_colorcode == label_method, ]
    # if there are no points with the specified colour code, label none of the points
    if(NROW(output_df) == 0) {
      label_points <- "none"
      label_method <- "allornone"
    }
    # user may have specified all points with that colour code, so if
    # label_count is NA, label all of them, otherwise label the top n points
    # with that colour code
    if(is.na(label_count)) {
      label_points <- output_df$Gene.names
    } else {
      last_row <- min(label_count, NROW(output_df))
      top_points <- output_df[1:last_row, ]
      label_points <- top_points$Gene.names
    }
  }

  # at this stage label_points will be either "none", "all", or a character
  # vector of gene names to label, so label the points accordingly

  if(label_method == "allornone") {
    if(label_points == "none") {
      df_subset$point_labels <- ""
    } else if(label_points == "all") {
      df_subset$point_labels <- df_subset$Gene.names
    }
  } else {
    # we should have a vector of gene names to label in label_points, so label
    # those points and leave the rest blank label_points may match Gene.names,
    # Protein.names, or Protein.IDs, so check all three columns for matches and
    # label the points accordingly
    df_subset$point_labels <- ifelse(df_subset$Gene.names %in% label_points, df_subset$Gene.names,
                                     ifelse(df_subset$Protein.names %in% label_points, df_subset$Protein.names,
                                            ifelse(df_subset$Protein.IDs %in% label_points, df_subset$Protein.IDs, "")))
  }
  return(df_subset)
}
