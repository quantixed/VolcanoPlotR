volcano_plot_maxquant <- function(df_subset = NULL, meas = "LFQ.intensity",
                                  threshold_p = 0.05, threshold_fc = 1,
                                  vp_colours = NULL, fsize = 8) {
  # check that meas ends in . and if not, add it
  if (!grepl("\\.$", meas)) {
    meas <- paste0(meas, ".")
  }

  if(!validate_df_subset(df_subset)) {
    stop("Data frame does not contain the required columns to generate the volcano plot.")
  }

  groups <- determine_groups_from_df_subset(df_subset, meas)

  group1 <- groups[1]
  group2 <- groups[2]

  # add a colour code column for the volcano plot we'll do it using bit setting:
  # bit 0 is p-value < threshold_p, bit 1 is meas.ratio <= -1, bit 2 is
  # meas.ratio >= 1 then change to character for the purposes of applying
  # colours to points
  df_subset$vp_colorcode <- ifelse(df_subset$p.value < threshold_p, 1, 0) +
    ifelse(df_subset$meas.ratio <= -threshold_fc, 2, 0) +
    ifelse(df_subset$meas.ratio >= threshold_fc, 4, 0)
  df_subset$vp_colorcode <- as.character(df_subset$vp_colorcode)

  if(is.null(vp_label)) {
    # integers 0 through 5 are possible
    vp_colours <- c("0" = "#a0a0a0", "1" = "#808080",
                    "2" = "#606060", "3" = "#8080ff",
                    "4" = "#606060", "5" = "#ff80ff")
  }

  p <- ggplot(df_subset, aes(x = meas.ratio, y = neg.log10.p.value)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey") +
    geom_hline(yintercept = -log10(threshold_p), linetype = "dashed", colour = "grey") +
    geom_point(aes(colour = vp_colorcode),
               size = 1, shape = 16, alpha = 0.5) +
    scale_colour_manual(values = vp_colours) +
    # scale_x_continuous(limits = c(-11, 11), breaks = seq(-10, 10, 2)) +
    # scale_y_continuous(limits = c(0, 20), breaks = seq(0, 20, 2),
    #                    labels = vp_label) +
    xlab(bquote(.(group1) ~ - ~ .(group2) ~ (Log[2]))) +
    ylab(bquote(P~value~(-Log[10]))) +
    theme_classic(fsize) +
    theme(legend.position = "none")

  return (p)
}




