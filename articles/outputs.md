# Outputs from VolcanoPlotR

The main output is obviously, the volcano plot itself, which is returned
as a ggplot object from the
[`volcano_plot_maxquant()`](https://quantixed.github.io/VolcanoPlotR/reference/volcano_plot_maxquant.md)
function. This can be customised further using ggplot2 functions.

## Text output

To generate a text output of the enriched proteins, you can set the
`text_output` parameter to `TRUE` in the
[`volcano_plot_maxquant()`](https://quantixed.github.io/VolcanoPlotR/reference/volcano_plot_maxquant.md)
function. This will generate a text file containing the enriched
proteins and their associated statistics. The output is saved as a
tab-separated text file in the `text_output_dir` directory (default is
“Output/Data/”).

It contains a ranked list of proteins where the ranking is calculated by
the manhattan distance from the origin (0,0) in the volcano plot. This
list can be useful for other applications or for generating a table for
publication.

## Alternative processing and statistical tests

### Changing the processing of MaxQuant data

By default,
[`process_maxquant()`](https://quantixed.github.io/VolcanoPlotR/reference/process_maxquant.md)
will analysis LFQ intensity values. To use a different measurement, the
`meas` parameter can be set to one of the other measurements in the
MaxQuant output (e.g. “iBAQ”, “Intensity”, “MS/MS.counts”).

The `baseval` parameter can be set to a value other than 0 to use a
different base value for log2 transformation. The `width` and
`downshift` parameters can be set to different values to change the
imputation of missing values. The `seed` parameter can be set to a
different value to change the random seed used for imputation.

### Statistical tests

The standard method for calculating p-values is via an unpaired two
sample t-test assuming equal variance. To change to paired t-test, you
can set the `paired` parameter to `TRUE` in the
[`process_maxquant()`](https://quantixed.github.io/VolcanoPlotR/reference/process_maxquant.md)
function, similarly the `var.equal` parameter can be set to `FALSE` to
use Welch’s t-test instead of Student’s t-test. The default settings
mirror Perseus processing.
