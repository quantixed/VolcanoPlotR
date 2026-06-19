# Mean Plot of Maxquant Data

This function generates a plot of the mean value of each protein for
Group 1 (y) vs Group 2 (x) from a processed set of MaxQuant data. The
data frame should be the same as one that would be passed to
[`volcano_plot_maxquant()`](https://quantixed.github.io/VolcanoPlotR/reference/volcano_plot_maxquant.md)
i.e. should contain the following columns: Protein.IDs, Gene.names,
Protein.names, meas.ratio, p.value, and neg.log10.p.value. The function
will add a colour code column to the data frame based on the p-value and
fold change thresholds, and then use ggplot2 to create the volcano plot.

## Usage

``` r
mean_plot_maxquant(
  df_subset = NULL,
  meas = "LFQ.intensity",
  threshold_p = 0.05,
  threshold_fc = 1,
  xy_line = FALSE,
  vp_colours = NULL,
  groups = NULL,
  x_label = NULL,
  y_label = NULL,
  fsize = 8,
  label_points = "none",
  point_args = list(size = 1, shape = 16, alpha = 0.5),
  label_args = list(size = 1.5, max.overlaps = 25, segment.alpha = 0.5, segment.size =
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

- threshold_p:

  numeric indicating the p-value threshold

- threshold_fc:

  numeric indicating the fold change threshold (in log2 space, i.e. 1 is
  a 2-fold change, 2 is a 4-fold change, etc.)

- xy_line:

  boolean indicating whether to add a line to indicate x = y

- vp_colours:

  a named vector of colours for the volcano plot. The names of the
  vector should be the integers 0 through 5, which correspond to the
  different combinations of p-value and fold change thresholds. If NULL,
  default colours will be used.

- groups:

  a vector of two strings containing the names of the experimental
  (these are used for the x-axis label). The purpose is to allow the
  user to substitute a better label than that used in the MaxQuant data.
  If NULL the values from the data are deduced and used.

- x_label:

  string specifying the label for the x-axis. If NULL, a default label
  will be generated based on the group names.

- y_label:

  string specifying the label for the y-axis. If NULL, a default label
  will be generated.

- fsize:

  numeric indicating the font size to use for the plot (default is 8)

- label_points:

  string indicating which points to label on the plot. Options are
  "none" (default) for no labels, "all" (attempt to label all points), a
  character vector of selected proteins to label (e.g. c("P12345",
  "Q67890", "IPO5") can be Gene.names, Protein.names, Protein.ID or a
  mix but the values will be used to label), or a code for
  auto-labelling (e.g. "3_all" or "5_10" which would label all
  significantly de-enriched points (colour code 3) or the top 10
  significantly enriched proteins (colour code 5), a code of "top_20"
  will label the top 20 proteins by manhattan distance regardless of
  colour code).

- point_args:

  a list of arguments to be passed to `geom_point()` for the points in
  the volcano plot (default is
  `list(size = 1, shape = 16, alpha = 0.5)`). If the user specifies
  additional arguments or changes a default argument, the defaults will
  be merged with the user-specified entries.

- label_args:

  a list of arguments to be passed to
  [`ggrepel::geom_text_repel()`](https://ggrepel.slowkow.com/reference/geom_text_repel.html)
  for the labels in the volcano plot (default is
  `list(size = 1.5, max.overlaps = 25, segment.alpha = 0.5, segment.size = 0.2, colour = "black")`).
  If the user specifies additional arguments or changes a default
  argument, the defaults will be merged with the user-specified entries.

## Value

ggplot object containing the volcano plot

## Examples

``` r
if (FALSE) { # \dontrun{
# load the MaxQuant data
df <- load_maxquant()
# process the MaxQuant data to get the subset for the volcano plot
df_subset <- process_maxquant(df, group1 = "Treatment", group2  = "Control")
# generate the volcano plot
mean_plot_maxquant(df_subset)
} # }
```
