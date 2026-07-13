#' Load multiple MaxQuant output files
#'
#' Load several `proteinGroups.txt` files from MaxQuant and perform some
#' cleaning steps to prepare the data for combined analysis. If only one file is
#' to be analyzed, use `load_maxquant()` instead. For automatic loading, place
#' `proteinGroups.txt` files in subfolders within `datadir`. Proteins that are
#' only identified by site, reverse hits, and potential contaminants are
#' removed, and missing gene and protein names are filled in based on the
#' "Fasta.headers" column. Note: in the case of multiple identified proteins,
#' the first one in the list is kept.
#'
#' @param file character string specifying the name of the MaxQuant output file
#'   to load (default is "proteinGroups.txt").
#' @param datadir character string specifying the directory where the MaxQuant
#'   output files are located (default is "Data"). The function performs a
#'   recursive search in this directory for all files matching the specified
#'   `file` name.
#' @param clean boolean indicating whether to perform cleaning steps on the data
#'   (default is TRUE). If TRUE, rows with a "+" in the
#'   "Only.identified.by.site", "Reverse", or "Potential.contaminant" columns
#'   will be removed, and missing gene and protein names will be filled in based
#'   on the "Fasta.headers" column.
#'
#' @returns A list of data frames containing the MaxQuant data.
#' @importFrom utils read.delim
#' @export
#'
#' @examples
#' \dontrun{
#' load_multiple_maxquant()
#' }
load_multiple_maxquant <- function(file = "proteinGroups.txt",
                                   datadir = "Data",
                                   clean = TRUE) {
  # check datadir exists
  if (!dir.exists(datadir)) {
    stop(paste("Directory not found:", datadir))
  }
  # recursive search for all files matching the specified file name
  files <- list.files(path = datadir, pattern = file,
                      recursive = TRUE, full.names = TRUE)
  # check if any files were found
  if (length(files) == 0) {
    stop(paste("No files found matching", file, "in directory", datadir))
  }

  # for each file found, load the data and store in a list
  data_list <- lapply(files, function(f) {
    # Read the MaxQuant output file
    df <- load_maxquant(filepath = f, clean = clean)
    df
  })
  # Return the loaded data list
  return(data_list)
}
