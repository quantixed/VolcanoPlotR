# Load multiple MaxQuant output files

Load several `proteinGroups.txt` files from MaxQuant and perform some
cleaning steps to prepare the data for combined analysis. If only one
file is to be analyzed, use
[`load_maxquant()`](https://quantixed.github.io/VolcanoPlotR/reference/load_maxquant.md)
instead. For automatic loading, place `proteinGroups.txt` files in
subfolders within `datadir`. Proteins that are only identified by site,
reverse hits, and potential contaminants are removed, and missing gene
and protein names are filled in based on the "Fasta.headers" column.
Note: in the case of multiple identified proteins, the first one in the
list is kept.

## Usage

``` r
load_multiple_maxquant(
  file = "proteinGroups.txt",
  datadir = "Data",
  clean = TRUE
)
```

## Arguments

- file:

  character string specifying the name of the MaxQuant output file to
  load (default is "proteinGroups.txt").

- datadir:

  character string specifying the directory where the MaxQuant output
  files are located (default is "Data"). The function performs a
  recursive search in this directory for all files matching the
  specified `file` name.

- clean:

  boolean indicating whether to perform cleaning steps on the data
  (default is TRUE). If TRUE, rows with a "+" in the
  "Only.identified.by.site", "Reverse", or "Potential.contaminant"
  columns will be removed, and missing gene and protein names will be
  filled in based on the "Fasta.headers" column.

## Value

A list of data frames containing the MaxQuant data.

## Examples

``` r
if (FALSE) { # \dontrun{
load_multiple_maxquant()
} # }
```
