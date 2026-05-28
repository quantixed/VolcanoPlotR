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
