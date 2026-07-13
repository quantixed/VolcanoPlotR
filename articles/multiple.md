# Analyzing multiple MaxQuant files with VolcanoPlotR

To combine more than one MaxQuant file, you can use the
[`load_multiple_maxquant()`](https://quantixed.github.io/VolcanoPlotR/reference/load_multiple_maxquant.md)
function to load them all into R.

The best approach is to place each `proteinGroups.txt` into its own
subdirectory in the `Data` directory; the command will find and load any
`proteinGroups.txt` in the `Data` directory.

You will be given the option of specifying which group(s) corresponds to
group1 (right side of the volcano plot) and which group(s) corresponds
to group2 (left side of the volcano plot). The function will then
combine the data from all the files and return a single data frame that
can be used for plotting.

In a simple example, two MaxQuant files have two groups each, “WT” and
“Control”. The runs in the both files are called WT.1, WT.2, WT.3,
Control.1, Control.2, Control.3. The function presents the user with a
list of the groups found in the files and asks which group(s) should be
used for group1 and which group(s) should be used for group2. If we
specify WT as group1 and Control as group2, the function will combine
them. In doing so it will rename the runs to test.1, test.2, test.3,
test.4, test5, test.6, control.1, control.2, control.3, control.4,
control.5, control.6 to avoid any duplicate names.

In a more complicated example, in another set, the groups are called
“WT” and “Control” in one file and “wildtype” and “control” in the
other. The user can select from a list to specify the groups in group1
and those in group2. They would see a list of 1 WT, 2 Control, 3
wildtype, 4 control. If they select “1 3” for group1 and “2 4” for
group2, the function will combine them and rename the runs to test.1,
test.2, test.3, test.4, test5, test.6, control.1, control.2, control.3,
control.4, control.5, control.6.

Note that the group names need to be unique and logical, for example if
WT and Control are mislabelled in one file, or if names are shared
between those in group1 and group2 then they cannot be combined. The
solution here is to rerun MaxQuant with appropriate labelling.

In both of these examples, the user can bypass the selection step by
specifying the groups in the function call, that is using
`load_multiple_maxquant(group1 = "WT", group2 = "Control")` for the
first example and
`load_multiple_maxquant(group1 = c("WT", "wildtype"), group2 = c("Control", "control"))`
for the second example.

## Processing and plotting

The next steps are the same as if a single MaxQuant file had been used.
