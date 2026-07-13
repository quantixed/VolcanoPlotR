#' Process multiple MaxQuant datasets for volcano plot
#'
#' This function takes a list of data frames containing MaxQuant output and
#' processes it to prepare for a volcano plot. It identifies the experimental
#' groups, allows the user to select which groups to compare, imputes missing
#' values, calculates log2 fold changes and p-values for each protein, and
#' returns a processed data frame ready for plotting. It uses equivalent
#' processing as Perseus, including imputation of missing values from a Gaussian
#' distribution with specified downshift and width parameters. Note, group1 is
#' compared to group2, so the log2 fold change is calculated as group1 - group2,
#' meaning that a positive value indicates enrichment in group1 and a negative
#' value indicates enrichment in group2.
#'
#' @param data list of data frame containing MaxQuant outputs, typically loaded
#'   using `load_multiple_maxquant()`.
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
#' @param seed numeric value specifying the random seed to use for imputation of
#'   missing values from a Gaussian distribution (default is 123). Setting a
#'   seed ensures that the imputation is reproducible.
#' @param var.equal boolean indicating whether to assume equal variances for the
#'   t-test (default is TRUE). This is typically TRUE for MaxQuant LFQ intensity
#'   data, but may be different for other types of data.
#' @param paired boolean indicating whether to perform a paired t-test (default
#'   is FALSE ). Note, paired assumes the columns for group1 and group2 are in
#'   the same order and correspond to each other. No check is made to match the
#'   columns by name.
#' @param ratio boolean indicating whether the enrichment in group 1 vs group 2
#'   is calculated for each protein for each run and then averaged (default is
#'   FALSE). If FALSE, the mean of the measure of a group 1 protein is compared
#'   with the mean of the measure for group 2. The ratio method is useful when
#'   files have been processed in different batches and the means of the groups
#'   are not comparable. The ratio method is more robust to batch effects, but
#'   it requires that the number of replicates in group 1 and group 2 are the
#'   same.
#'
#' @returns A data frame containing the processed MaxQuant data suitable for
#'   volcano plot visualization, including columns for log2 fold change and
#'   p-values.
#' @importFrom stats rnorm t.test sd
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # load the MaxQuant data
#' df_list <- load_multiple_maxquant()
#' # process the MaxQuant data to get the subset for the volcano plot
#' df_subset <- process_multiple_maxquant(data = df_list, group1 = "Treatment", group2  = "Control")
#' }
process_multiple_maxquant <- function(data = NULL,
                                      group1 = NULL, group2 = NULL,
                                      meas = "LFQ.intensity", baseval = 0,
                                      width = 0.3, downshift = 1.8,
                                      seed = 123,
                                      var.equal = TRUE,
                                      paired = FALSE,
                                      ratio = FALSE) {

  # check for list and if not present, stop and print an error message
  if (is.null(data)) {
    stop("No list of data frames supplied.
    Please provide a list of data frames containing the MaxQuant data to be
         processed.")
  } else {
    # do we have a list of data frames or a single data frame?
    if (is.data.frame(data)) {
      stop("A single data frame has been supplied.
         Please provide a list of data frames containing the MaxQuant data to be
         processed.")
    } else if (!is.list(data)) {
      stop("The supplied data is not a list of data frames.
         Please provide a list of data frames containing the MaxQuant data to be
         processed.")
    }
  }

  # skip the check that a processed data frame has passed accidentally, will add
  # back if required

  # check that meas ends in . and if not, add it
  if (!grepl("\\.$", meas)) {
    meas <- paste0(meas, ".")
  }

  # print the groups found for each data frame in the list
  all_groups <- lapply(data, function(x){
    unique(gsub(meas, "", grep(meas, colnames(x), value = TRUE)))
  })
  cat("The group-measures in each data frame are:\n")
  print(all_groups)

  groups <- lapply(all_groups, function(x){
    # these groups should be foo.1, foo.2, and bar.baz.1, bar.baz.2 so get the
    # group names by removing the .1 and .2 (or any integer) from the group
    # names
    grps <- unique(gsub("[0-9]+$", "", x))
    # if there are trailing "." or "_" in the group names, remove them
    grps <- gsub("[._]$", "", grps)
  })

  cat("The groups are:\n")
  print(groups)

  groups <- unique(unlist(groups))

  # if group1 and group2 are already specified check that they are in groups
  if (is.null(group1) || is.null(group2)) {
    cat("No groups specified for comparison. Please select from the following:\n")
    for (i in seq_along(groups)) {
      cat(i, ":", groups[i], "\n")
    }
    group1_index <- as.integer(scan(text = readline(
      prompt = "Enter the number(s) corresponding to the first group (separate with a space): ")))
    group2_index <- as.integer(scan(text = readline(
      prompt = "Enter the number(s) corresponding to the second group (separate with a space): ")))
    group1 <- groups[group1_index]
    group2 <- groups[group2_index]
    cat("You have selected:", group1, "versus", group2, "\n")
  } else {
    # test all items in group1 and group2 are present in groups
    if (all(c(group1,group2) %in% groups)) {
      cat("Using specified groups:", group1, "versus", group2, "\n")
    } else {
      cat("Specified groups not found in data. Please select from the following:\n")
      for (i in seq_along(groups)) {
        cat(i, ":", groups[i], "\n")
      }
      group1_index <- as.integer(scan(text = readline(
        prompt = "Enter the number(s) corresponding to the first group (separate with a space): ")))
      group2_index <- as.integer(scan(text = readline(
        prompt = "Enter the number(s) corresponding to the second group (separate with a space): ")))
      group1 <- groups[group1_index]
      group2 <- groups[group2_index]
      cat("You have selected:", group1, "versus", group2, "\n")
    }
  }

  # It is unlikely that the user will select one group1 name and one group2 name AND that the repeats are uniquely labelled
  # We probably have foo.1, foo.2  in one data frame and foo.1, foo.2 in another. So we need to rename and relabel them
  g1_names <- grep(paste(group1, collapse = "|"), unlist(all_groups), value = TRUE)
  g2_names <- grep(paste(group2, collapse = "|"), unlist(all_groups), value = TRUE)
  g1_aliases <- c(paste("test", 1:length(g1_names), sep = "."))
  g2_aliases <- c(paste("control", 1:length(g2_names), sep = "."))

  # walk through the list of data frames and rename the columns for group1 and group2 to the aliases
  # it would be more efficient to generate 1 name and 1 alias list and do this once, but I'm struggling to see how to do this
  for (i in seq_along(g1_names)) {
    for(j in 1:length(data)) {
      if (!is.na(grep(g1_names[i], colnames(data[[j]]))[1])) {
        colnames(data[[j]]) <- gsub(g1_names[i], g1_aliases[i], colnames(data[[j]]))
        break
      }
    }
  }

  for (i in seq_along(g2_names)) {
    for(j in 1:length(data)) {
      if (!is.na(grep(g2_names[i], colnames(data[[j]]))[1])) {
        colnames(data[[j]]) <- gsub(g2_names[i], g2_aliases[i], colnames(data[[j]]))
        break
      }
    }
  }

  # rename group1 and group2 now
  cat("Renaming group1 and group2 to 'test' and 'control' for processing.\n")
  group1 <- "test"
  group2 <-  "control"

  # create a subset of the data frames that contains only the columns
  # "Gene.names", "Protein.names", and the columns that contain "LFQ.intensity."
  # (or the meas name) with test or control in their names
  list_subset <- lapply(data, function(df) {
    df_subset <- df[, c(
      "Protein.IDs", "Gene.names", "Protein.names",
      grep(paste0(meas, group1), colnames(df), value = TRUE),
      grep(paste0(meas, group2), colnames(df), value = TRUE)
    )]
  })
  # merge all data frames in list
  df_subset <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = c(
    "Protein.IDs", "Gene.names", "Protein.names"), all.x = TRUE),
    list_subset)
  # replace all NA values with baseval
  df_subset[is.na(df_subset)] <- baseval
  # because the file may have more conditions, while we only compare 2, it is
  # possible that a protein is only identified in the other conditions, so we
  # need to remove rows that have all NA values in the columns that contain
  # "LFQ.intensity." with group1 or 2 in their names remove rows that have all 0
  # values in the columns that contain "LFQ.intensity." with group1 or 2 in
  # their names
  df_subset <- df_subset[!apply(df_subset[, grep(meas, colnames(df_subset))],
                                1, function(x) all(x == baseval)), ]
  # now we will impute the base value working on each LFQ.intensity column in
  # turn, mask off the 0 values using conversion to NA then take the sd and mean
  # of the non-NA values in the column and then use this to impute values from a
  # Gaussian back to the NA values
  set.seed(seed)
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
  group1_data <- df_subset[, grep(paste0(meas, group1), colnames(df_subset))]
  group2_data <- df_subset[, grep(paste0(meas, group2), colnames(df_subset))]
  # check that the number of columns in group1_data and group2_data are the same
  if(ncol(group1_data) != ncol(group2_data)) {
    unequal <- TRUE
  }
  if(ratio & unequal) {
    warning("The number of columns in group1_data and group2_data are not the same. The ratio will be calculated as the mean of group1_data minus the mean of group2_data.")
    ratio <- FALSE
  }
  if(paired & unequal) {
    warning("The number of columns in group1_data and group2_data are not the same. Paired t-test will not be used.")
    paired <- FALSE
  }
  # calculate the mean for each row of the LFQ.intensity.group1 as a new column and the mean of LFQ.intensity.group2 as a new column
  df_subset$meas.group1.mean <- rowMeans(group1_data, na.rm = TRUE)
  df_subset$meas.group2.mean <- rowMeans(group2_data, na.rm = TRUE)
  if(ratio) {
    # the ratio is group1 - group2 because we are in log2 space
    df_subset$meas.ratio <- rowMeans(group1_data - group2_data, na.rm = TRUE)
  } else {
    df_subset$meas.ratio <- df_subset$meas.group1.mean - df_subset$meas.group2.mean
  }

  ## p-value calculation

  # for each row, do a t-test between the values for LFQ.intensity.group1 and LFQ.intensity.group2
  for (i in 1:nrow(df_subset)) {
    x <- group1_data[i, ]
    y <- group2_data[i, ]
    df_subset$p.value[i] <- t.test(x, y,
                                   var.equal = var.equal, paired = paired)$p.value
  }
  # do a -log10 transformation of the p-values and add it as a new column
  df_subset$neg.log10.p.value <- -log10(df_subset$p.value)
  # calculate the Manhattan distance for each row from the origin (0,0) in the volcano plot space and add it as a new column
  df_subset$manhattan.distance <- abs(df_subset$meas.ratio) + abs(df_subset$neg.log10.p.value)

  return(df_subset)
}
