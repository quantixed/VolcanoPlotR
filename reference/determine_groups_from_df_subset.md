# Determine the experimental groups in a data frame subset

This function takes a data frame containing a subset of the MaxQuant
data and determines the experimental groups that are present in the data
frame by looking at the column names.

## Usage

``` r
determine_groups_from_df_subset(df_subset, meas)
```

## Arguments

- df_subset:

  a data frame containing the subset of the MaxQuant data that will be
  used to generate the volcano plot.

- meas:

  string indicating the measurement type

## Value

vector of two strings containing the names of the experimental groups
