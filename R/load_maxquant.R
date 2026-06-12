#' Load MaxQuant output file
#'
#' Load a `proteinGroups.txt` file from MaxQuant and perform some cleaning steps
#' to prepare the data for analysis. Proteins that are only identified by site,
#' reverse hits, and potential contaminants are removed, and missing gene and
#' protein names are filled in based on the "Fasta.headers" column. Note: in the
#' case of multiple identified proteins, the first one in the list is kept.
#'
#' @param file character string specifying the name of the MaxQuant output file
#'   to load (default is "proteinGroups.txt").
#' @param datadir character string specifying the directory where the MaxQuant
#'   output file is located (default is "Data").
#' @param clean boolean indicating whether to perform cleaning steps on the data
#'   (default is TRUE). If TRUE, rows with a "+" in the
#'   "Only.identified.by.site", "Reverse", or "Potential.contaminant" columns
#'   will be removed, and missing gene and protein names will be filled in based
#'   on the "Fasta.headers" column.
#'
#' @returns A data frame containing the MaxQuant data.
#' @importFrom utils read.delim
#' @export
#'
#' @examples
#' \dontrun{
#' load_maxquant()
#' }
load_maxquant <- function(file = "proteinGroups.txt", datadir = "Data", clean = TRUE) {
  # construct the full file path
  filepath <- file.path(datadir, file)
  # chack file exists
  if (!file.exists(filepath)) {
    stop(paste("File not found:", filepath))
  }
  # Read the MaxQuant output file
  df <- read.delim(filepath, stringsAsFactors = FALSE)

  # if requested, remove rows that have a "+" in either of these columns:
  # "Only.identified.by.site", "Reverse", "Potential.contaminant"
  # ensure columns are character before checking
  if(clean) {
    cols_of_interest <- c("Only.identified.by.site",
                          "Reverse",
                          "Potential.contaminant")
    df[, cols_of_interest] <- lapply(df[, cols_of_interest], as.character)
    # if there is a "+" in any of the columns, remove the row
    df <- df[!apply(df[, cols_of_interest], 1,
                    function(x) any(grepl("\\+", x))), ]
  }

  # if there are any blanks in Gene.names, replace them with the text that is
  # after "GN=" in the "Fasta.headers" column
  df$Gene.names <- ifelse(
    df$Gene.names == "",
    sapply(df$Fasta.headers, function(x) {
      match <- regexpr("GN=([^ ]+)", x)
      if (match != -1) {
        return(substr(x, match + 3, match + attr(match, "match.length") - 1))
      } else {
        return("")
      }
    }), df$Gene.names)
  # if there are any rows in Gene.names that have a ";" in them, keep only the
  # first gene name (before the ";")
  df$Gene.names <- sapply(df$Gene.names,
                          function(x) strsplit(x,";")[[1]][1])
  # if there are any blanks in Protein.names, replace them with the text that is
  # after the 1st space and up to the occurrence of " OS=" in the "Fasta.headers"
  # column
  df$Protein.names <- ifelse(df$Protein.names == "",
                             sapply(df$Fasta.headers, function(x) {
                               start <- regexpr(" ", x)
                               stop <- regexpr(" OS=", x)
                               if (stop != -1) {
                                 return(substr(x, start + 1, stop - 1))
                               } else {
                                 return("")
                               }
                             }), df$Protein.names)
  # if there are any rows in Protein.names that have a ";" in them, keep only
  # the first gene name (before the ";")
  df$Protein.names <- sapply(df$Protein.names,
                             function(x) strsplit(x,";")[[1]][1])
  # if there are any rows in Protein.IDs that have a ";" in them, keep only the
  # first gene name (before the ";")
  df$Protein.IDs <- sapply(df$Protein.IDs,
                           function(x) strsplit(x,";")[[1]][1])

  # Return the loaded data
  return(df)
}
