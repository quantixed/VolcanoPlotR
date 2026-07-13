# PCA Plot of MaxQuant Sample Columns

This function performs principal component analysis on the measurement
columns present in a processed MaxQuant subset and plots PC1 against
PC2, with each sample column shown as a point. The input data frame
should be the same processed subset that would be passed to
[`volcano_plot_maxquant()`](https://quantixed.github.io/VolcanoPlotR/reference/volcano_plot_maxquant.md).

## Usage

``` r
pca_plot_maxquant(
  df_subset = NULL,
  meas = "LFQ.intensity",
  groups = NULL,
  by_protein = FALSE,
  threshold_p = 0.05,
  threshold_fc = 1,
  vp_colours = NULL,
  x_label = NULL,
  y_label = NULL,
  xy_line = TRUE,
  fsize = 8,
  label_points = NULL,
  point_args = list(size = 3, alpha = 0.8),
  label_args = list(size = 3, max.overlaps = Inf, segment.alpha = 0.5, segment.size =
    0.2, colour = "black")
)
```

## Arguments

- df_subset:

  a data frame containing the subset of the MaxQuant data that will be
  used to generate the plot

- meas:

  string indicating the measurement type (e.g. "LFQ.intensity",
  "Intensity", etc.)

- groups:

  optional character vector used to relabel the experimental groups in
  the legend. If `NULL`, group names are deduced from the measurement
  column names.

- by_protein:

  boolean indicating whether to perform PCA on the protein rows (default
  is FALSE). If TRUE, PCA will be performed on the protein rows rather
  than the sample columns.

- threshold_p:

  numeric indicating the p-value threshold. Used if by_protein is TRUE
  to colour the points in the PCA plot based on significance.

- threshold_fc:

  numeric indicating the fold change threshold (in log2 space, i.e. 1 is
  a 2-fold change, 2 is a 4-fold change, etc.). Used if by_protein is
  TRUE to colour the points in the PCA plot based on significance.

- vp_colours:

  a named vector of colours for the PCA plot.

- x_label:

  string specifying the label for the x-axis. If `NULL`, a default label
  showing the variance explained by PC1 will be generated.

- y_label:

  string specifying the label for the y-axis. If `NULL`, a default label
  showing the variance explained by PC2 will be generated.

- xy_line:

  boolean indicating whether to add a line to indicate x = y = 0

- fsize:

  numeric indicating the font size to use for the plot (default is 8)

- label_points:

  boolean indicating whether to label sample points (default is NULL and
  triggers a setting of FALSE if `by_protein` is TRUE, and TRUE if
  `by_protein` is is FALSE). If TRUE, sample points will be labelled
  with their sample names. If FALSE, no labels will be added.

- point_args:

  a list of arguments to be passed to `geom_point()` for the points in
  the PCA plot (default is `list(size = 3, alpha = 0.8)`). If the user
  specifies additional arguments or changes a default argument, the
  defaults will be merged with the user-specified entries.

- label_args:

  a list of arguments to be passed to
  [`ggrepel::geom_text_repel()`](https://ggrepel.slowkow.com/reference/geom_text_repel.html)
  for the labels in the PCA plot (default is
  `list(size = 3, max.overlaps = Inf, segment.alpha = 0.5, segment.size = 0.2, colour = "black")`).
  If the user specifies additional arguments or changes a default
  argument, the defaults will be merged with the user-specified entries.

## Value

ggplot object containing the PCA plot

## Examples

``` r
if (FALSE) { # \dontrun{
# load the MaxQuant data
df <- load_maxquant(file = "proteinGroups.txt", datadir = "inst/extdata")
# process the MaxQuant data to get the subset for plotting
df_subset <- process_maxquant(df, group1 = "WT", group2 = "Control")
# generate the PCA plot
pca_plot_maxquant(df_subset)
} # }
```
