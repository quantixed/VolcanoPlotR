# Process MaxQuant data for volcano plot

This function takes a data frame containing MaxQuant output and
processes it to prepare for a volcano plot. It identifies the
experimental groups, allows the user to select which groups to compare,
imputes missing values, calculates log2 fold changes and p-values for
each protein, and returns a processed data frame ready for plotting. It
uses equivalent processing as Perseus, including imputation of missing
values from a Gaussian distribution with specified downshift and width
parameters. Note, group1 is compared to group2, so the log2 fold change
is calculated as group1 - group2, meaning that a positive value
indicates enrichment in group1 and a negative value indicates enrichment
in group2.

## Usage

``` r
process_maxquant(
  df = NULL,
  group1 = NULL,
  group2 = NULL,
  meas = "LFQ.intensity",
  baseval = 0,
  width = 0.3,
  downshift = 1.8,
  seed = 123,
  var.equal = TRUE,
  paired = FALSE
)
```

## Arguments

- df:

  data frame containing MaxQuant output, typically loaded using
  [`load_maxquant()`](https://quantixed.github.io/VolcanoPlotR/reference/load_maxquant.md).

- group1:

  character string specifying the name of the first experimental group
  to compare (e.g., "Treatment"). If NULL, the user will be prompted to
  select from the groups found in the data frame.

- group2:

  character string specifying the name of the second experimental group
  to compare (e.g., "Control"). If NULL, the user will be prompted to
  select from the groups found in the data frame.

- meas:

  character string specifying the prefix of the measurement columns to
  use for the comparison (default is "LFQ.intensity").

- baseval:

  numeric value specifying the base value to impute for missing values
  (default is 0). This is typically 0 for MaxQuant LFQ intensity data,
  but may be different for other types of data.

- width:

  numeric value specifying the width parameter for imputation of missing
  values from a Gaussian distribution (default is 0.3). This is
  typically 0.3 for MaxQuant LFQ intensity data, but may be different
  for other types of data.

- downshift:

  numeric value specifying the downshift parameter for imputation of
  missing values from a Gaussian distribution (default is 1.8). This is
  typically 1.8 for MaxQuant LFQ intensity data, but may be different
  for other types of data.

- seed:

  numeric value specifying the random seed to use for imputation of
  missing values from a Gaussian distribution (default is 123). Setting
  a seed ensures that the imputation is reproducible.

- var.equal:

  boolean indicating whether to assume equal variances for the t-test
  (default is TRUE). This is typically TRUE for MaxQuant LFQ intensity
  data, but may be different for other types of data.

- paired:

  boolean indicating whether to perform a paired t-test (default is
  FALSE ). Note, paired assumes the columns for group1 and group2 are in
  the same order and correspond to each other. No check is made to match
  the columns by name.

## Value

A data frame containing the processed MaxQuant data suitable for volcano
plot visualization, including columns for log2 fold change and p-values.

## Examples

``` r
if (FALSE) { # \dontrun{
# load the MaxQuant data
df <- load_maxquant()
# process the MaxQuant data to get the subset for the volcano plot
df_subset <- process_maxquant(df, group1 = "Treatment", group2  = "Control")
} # }
```
