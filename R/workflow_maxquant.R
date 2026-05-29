#' Workflow for MaxQuant data analysis
#'
#' The purpose of this function is to provide a single function that will take a
#' `proteinGroups.txt` produced by MaxQuant, perform the necessary data
#' processing, and generate a volcano plot. This function will call the
#' `load_maxquant()`, `process_maxquant()`, and `volcano_plot_maxquant()`
#' functions in turn. The user will be asked to specify which groups are group 1
#' or group 2.
#'
#' @param ... ellipsis to allow the user to specify arguments for each function
#'   in the workflow. The accepted arguments for each function are as follows:
#' \itemize{
#'  \item `load_maxquant()`: `file`, `datadir`, `clean`
#'  \item `process_maxquant()`: `group1`, `group2`, `meas`, `baseval`, `width`,
#'  `downshift`, `seed`, `var.equal`, `paired`
#'  \item `volcano_plot_maxquant()`: `meas`, `threshold_p`, `threshold_fc`,
#'  `p_line`, `zero_line`, `x_line`, `vp_colours`, `groups`, `x_label`,
#'  `y_label`, `fsize`
#'  }
#'
#' @returns A ggplot object containing the volcano plot generated from the
#'   MaxQuant data.
#' @export
#'
#' @examples
#' \dontrun{
#' workflow_maxquant()
#' }
workflow_maxquant <- function(...) {
  # get ellipsis arguments
  args <- list(...)
  # R CMD CHECK satisfaction
  raw_data <- processed_data <- p <- NULL

  ## Loading
  accepted_args <- c("file", "datadir", "clean")
  filtered_args <- args[names(args) %in% accepted_args]
  if(length(filtered_args) > 0) {
    raw_data <- do.call(load_maxquant, filtered_args)
  } else {
    raw_data <- load_maxquant()
  }
  # check raw_data is not null before continuing
  if(is.null(raw_data)) {
    stop("Error: raw_data is null")
  }

  ## Processing
  accepted_args <- c("group1", "group2", "meas", "baseval", "width", "downshift",
                     "seed", "var.equal", "paired")
  filtered_args <- args[names(args) %in% accepted_args]
  if(length(filtered_args) > 0) {
    processed_data <- do.call(process_maxquant, c(list(df = raw_data), filtered_args))
  } else {
    processed_data <- process_maxquant(df = raw_data)
  }
  # check processed_data is not null before continuing
  if(is.null(processed_data)) {
    stop("Error: processed_data is null")
  }

  ## Volcano plot
  accepted_args <- c("meas", "threshold_p", "threshold_fc",
                     "p_line", "zero_line", "x_line", "vp_colours", "groups",
                     "x_label", "y_label", "fsize")
  filtered_args <- args[names(args) %in% accepted_args]
  if(length(filtered_args) > 0) {
    p <- do.call(volcano_plot_maxquant, c(list(df = processed_data), filtered_args))
  } else {
    p <- volcano_plot_maxquant(df_subset = processed_data)
  }
  # check p is not null before returning
  if(is.null(p)) {
    stop("Error: volcano plot is null")
  }

  return(p)
}
