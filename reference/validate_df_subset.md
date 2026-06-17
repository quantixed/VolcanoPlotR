# Validate a data frame to generate a volcano plot

Checks if we have all the necessary components to generate a volcano
plot.

## Usage

``` r
validate_df_subset(df_subset)
```

## Arguments

- df_subset:

  data frame containing the subset of the MaxQuant data that will be
  used to generate the volcano plot. This data frame should contain the
  columns "Protein.IDs", "Gene.names", "Protein.names", "meas.ratio",
  "p.value", and "neg.log10.p.value". The "Gene.names" and
  "Protein.names" columns are optional, but if they are not present, a
  warning will be printed.

## Value

A boolean indicating whether the data frame contains the required
columns to generate a volcano plot, and a warning if the optional
columns are not present.
