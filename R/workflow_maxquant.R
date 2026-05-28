#' Workflow for MaxQuant data analysis
#'
#' The purpose of this function is to provide a single function that will take a
#' `proteinGroups.txt` produced by MaxQuant, perform the necessary data
#' processing, and generate a volcano plot. This function will call the
#' `load_maxquant()`, `process_maxquant()`, and `volcano_plot_maxquant()`
#' functions in turn. The user will be asked to specify which groups are group 1 or group 2.
#'
#' @returns A ggplot object containing the volcano plot generated from the MaxQuant data.
#' @export
#'
#' @examples
#' \dontrun{
#' workflow_maxquant()
#' }
workflow_maxquant <- function() {
  # R CMD CHECK satisfaction
  raw_data <- processed_data <- p <- NULL
  # perform check that there is a file called "proteinGroups.txt" in the "Data" directory and if not, stop and print an error message
  if (!file.exists("Data/proteinGroups.txt")) {
    stop("File not found: Data/proteinGroups.txt. Please make sure the file is
 in the correct location and try again.")
  }
  raw_data <- load_maxquant()
  processed_data <- process_maxquant(raw_data)
  p <- volcano_plot_maxquant(processed_data)
  return(p)
}
