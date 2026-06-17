# Load MaxQuant output file

Load a `proteinGroups.txt` file from MaxQuant and perform some cleaning
steps to prepare the data for analysis. Proteins that are only
identified by site, reverse hits, and potential contaminants are
removed, and missing gene and protein names are filled in based on the
"Fasta.headers" column. Note: in the case of multiple identified
proteins, the first one in the list is kept.

## Usage

``` r
load_maxquant(file = "proteinGroups.txt", datadir = "Data", clean = TRUE)
```

## Arguments

- file:

  character string specifying the name of the MaxQuant output file to
  load (default is "proteinGroups.txt").

- datadir:

  character string specifying the directory where the MaxQuant output
  file is located (default is "Data").

- clean:

  boolean indicating whether to perform cleaning steps on the data
  (default is TRUE). If TRUE, rows with a "+" in the
  "Only.identified.by.site", "Reverse", or "Potential.contaminant"
  columns will be removed, and missing gene and protein names will be
  filled in based on the "Fasta.headers" column.

## Value

A data frame containing the MaxQuant data.

## Examples

``` r
if (FALSE) { # \dontrun{
load_maxquant()
} # }
```
