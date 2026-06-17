# Add a label column to the data frame subset for the volcano plot

Add a label column to the data frame subset for the volcano plot

## Usage

``` r
add_label_column(df_subset, label_points)
```

## Arguments

- df_subset:

  the data frame to be used for the volcano plot, which should contain
  the columns "Gene.names", "Protein.names", and "Protein.IDs" and a
  column "vp_colorcode" which contains the colour code for each point in
  the volcano plot (0, 1, 2, 3, 4, or 5)

- label_points:

  character string specifying which points to label in the volcano plot.
  This can be one of the following: "none", "all", a character vector of
  selected proteins to label, a code for auto-labelling

## Value

the input data frame with an additional column "point_labels" containing
the labels for the points to be labelled in the volcano plot. The
"point_labels"
