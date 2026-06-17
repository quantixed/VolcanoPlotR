# Workflow for MaxQuant data analysis

The purpose of this function is to provide a single function that will
take a `proteinGroups.txt` produced by MaxQuant, perform the necessary
data processing, and generate a volcano plot. This function will call
the
[`load_maxquant()`](https://quantixed.github.io/VolcanoPlotR/reference/load_maxquant.md),
[`process_maxquant()`](https://quantixed.github.io/VolcanoPlotR/reference/process_maxquant.md),
and
[`volcano_plot_maxquant()`](https://quantixed.github.io/VolcanoPlotR/reference/volcano_plot_maxquant.md)
functions in turn. The user will be asked to specify which groups are
group 1 or group 2.

## Usage

``` r
workflow_maxquant(...)
```

## Arguments

- ...:

  ellipsis to allow the user to specify arguments for each function in
  the workflow. The accepted arguments for each function are as follows:

  - [`load_maxquant()`](https://quantixed.github.io/VolcanoPlotR/reference/load_maxquant.md):
    `file`, `datadir`, `clean`

  - [`process_maxquant()`](https://quantixed.github.io/VolcanoPlotR/reference/process_maxquant.md):
    `group1`, `group2`, `meas`, `baseval`, `width`, `downshift`, `seed`,
    `var.equal`, `paired`

  - [`volcano_plot_maxquant()`](https://quantixed.github.io/VolcanoPlotR/reference/volcano_plot_maxquant.md):
    `meas`, `threshold_p`, `threshold_fc`, `p_line`, `zero_line`,
    `x_line`, `vp_colours`, `groups`, `x_label`, `y_label`, `fsize`

## Value

A ggplot object containing the volcano plot generated from the MaxQuant
data.

## Examples

``` r
if (FALSE) { # \dontrun{
workflow_maxquant()
} # }
```
