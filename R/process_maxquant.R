#' Process MaxQuant data for volcano plot
#'
#' This function takes a data frame containing MaxQuant output and processes it
#' to prepare for a volcano plot. It identifies the experimental groups, allows
#' the user to select which groups to compare, imputes missing values,
#' calculates log2 fold changes and p-values for each protein, and returns a
#' processed data frame ready for plotting. It uses equivalent processing as
#' Perseus, including imputation of missing values from a Gaussian distribution
#' with specified downshift and width parameters. Note, group1 is compared to
#' group2, so the log2 fold change is calculated as group1 - group2, meaning
#' that a positive value indicates enrichment in group1 and a negative value
#' indicates enrichment in group2.
#'
#' @param df data frame containing MaxQuant output, typically loaded using
#'   `load_maxquant()`.
#' @param group1 character string specifying the name of the first experimental
#'   group to compare (e.g., "Treatment"). If NULL, the user will be prompted to
#'   select from the groups found in the data frame.
#' @param group2 character string specifying the name of the second experimental
#'   group to compare (e.g., "Control"). If NULL, the user will be prompted to
#'   select from the groups found in the data frame.
#' @param meas character string specifying the prefix of the measurement columns
#'   to use for the comparison (default is "LFQ.intensity").
#' @param baseval numeric value specifying the base value to impute for missing
#'   values (default is 0). This is typically 0 for MaxQuant LFQ intensity data,
#'   but may be different for other types of data.
#' @param width numeric value specifying the width parameter for imputation of
#'   missing values from a Gaussian distribution (default is 0.3). This is
#'   typically 0.3 for MaxQuant LFQ intensity data, but may be different for
#'   other types of data.
#' @param downshift numeric value specifying the downshift parameter for
#'   imputation of missing values from a Gaussian distribution (default is 1.8).
#'   This is typically 1.8 for MaxQuant LFQ intensity data, but may be different
#'   for other types of data.
#'
#' @returns A data frame containing the processed MaxQuant data suitable for
#'   volcano plot visualization, including columns for log2 fold change and
#'   p-values.
#' @export
#'
#' @examples
process_maxquant <- function(df = NULL,
                             group1 = NULL, group2 = NULL,
                             meas = "LFQ.intensity", baseval = 0,
                             width = 0.3, downshift = 1.8) {
  # check that meas ends in . and if not, add it
  if (!grepl("\\.$", meas)) {
    meas <- paste0(meas, ".")
  }

  # find the experimental groups in the data frame using the names of the columns
  groups <- unique(gsub(meas, "", grep(meas, colnames(df), value = TRUE)))
  # if no groups are found, stop and print an error message
  if (length(groups) == 0) {
    stop(paste("No groups found in data frame with meas =", meas))
  } else {
    # these groups should be foo.1, foo.2, and bar.baz.1, bar.baz.2 so get the
    # group names by removing the .1 and .2 (or any integer) from the group
    # names
    groups <- unique(gsub("[0-9]+$", "", groups))
    # if there are trailing "." or "_" in the group names, remove them
    groups <- gsub("[._]$", "", groups)
  }

  # if group1 and group2 are already specified check that they are in groups
  if (is.null(group1) || is.null(group2)) {
    cat("No groups specified for comparison. Please select from the following:\n")
    for (i in seq_along(groups)) {
      cat(i, ":", groups[i], "\n")
    }
    group1_index <- as.integer(readline(prompt = "Enter the number corresponding to the first group to compare: "))
    group2_index <- as.integer(readline(prompt = "Enter the number corresponding to the second group to compare: "))
    group1 <- groups[group1_index]
    group2 <- groups[group2_index]
    cat("You have selected:", group1, "versus", group2, "\n")
  } else {
    if (group1 %in% groups && group2 %in% groups) {
      cat("Using specified groups:", group1, "versus", group2, "\n")
    } else {
      cat("Specified groups not found in data. Please select from the following:\n")
      for (i in seq_along(groups)) {
        cat(i, ":", groups[i], "\n")
      }
      group1_index <- as.integer(readline(prompt = "Enter the number corresponding to the first group to compare: "))
      group2_index <- as.integer(readline(prompt = "Enter the number corresponding to the second group to compare: "))
      group1 <- groups[group1_index]
      group2 <- groups[group2_index]
      cat("You have selected:", group1, "versus", group2, "\n")
    }
  }

  # create a subset of the data frame that contains only the columns
  # "Gene.names", "Protein.names", and the columns that contain "LFQ.intensity."
  # (or the meas name) with group1 or 2 in their names
  df_subset <- df[, c(
    "Protein.IDs", "Gene.names", "Protein.names",
    grep(paste0(meas, group1), colnames(df), value = TRUE),
    grep(paste0(meas, group2), colnames(df), value = TRUE)
  )]
  # because the file may have more conditions, while we only compare 2, it is
  # possible that a protein is only identified in the other conditions, so we
  # need to remove rows that have all NA values in the columns that contain
  # "LFQ.intensity." with group1 or 2 in their names remove rows that have all 0
  # values in the columns that contain "LFQ.intensity." with group1 or 2 in
  # their names
  df_subset <- df_subset[!apply(df_subset[, grep(meas, colnames(df_subset))], 1, function(x) all(x == baseval)), ]
  # now we will impute the base value working on each LFQ.intensity column in
  # turn, mask off the 0 values using conversion to NA then take the sd and mean
  # of the non-NA values in the column and then use this to impute values from a
  # gaussian back to the NA values
  df_subset[, grep(meas, colnames(df_subset))] <- lapply(df_subset[, grep(meas, colnames(df_subset))], function(x) {
    x <- as.numeric(x)
    x[x == baseval] <- NA
    # log2 transform the data
    x <- log2(x)
    sd_x <- sd(x, na.rm = TRUE)
    mean_x <- mean(x, na.rm = TRUE) - (sd_x * downshift)
    sd_x <- sd_x * width
    x[is.na(x)] <- rnorm(sum(is.na(x)), mean = mean_x, sd = sd_x)
    return(x)
  })

  ## ratio calculation (log2FC)

  # calculate the mean for each row of the LFQ.intensity.group1 as a new column and the mean of LFQ.intensity.group2 as a new column
  df_subset$meas.group1.mean <- rowMeans(df_subset[, grep(paste0(meas, group1), colnames(df_subset))], na.rm = TRUE)
  df_subset$meas.group2.mean <- rowMeans(df_subset[, grep(paste0(meas, group2), colnames(df_subset))], na.rm = TRUE)
  # the ratio is group1 - group2 because we are in log2 space
  df_subset$meas.ratio <- df_subset$meas.group1.mean - df_subset$meas.group2.mean

  ## p-value calculation

  # for each row, do a t-test between the values for LFQ.intensity.group1 and LFQ.intensity.group2
  for (i in 1:nrow(df_subset)) {
    x <- df_subset[i, grep(paste0(meas, group1), colnames(df_subset))]
    y <- df_subset[i, grep(paste0(meas, group2), colnames(df_subset))]
    df_subset$p.value[i] <- t.test(x, y, var.equal = TRUE)$p.value
  }
  # do a -log10 transformation of the p-values and add it as a new column
  df_subset$neg.log10.p.value <- -log10(df_subset$p.value)

  return(df_subset)
}
