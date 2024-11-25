# Clear
rm(list = ls(all.names = TRUE)) # will clear all (including hidden) objects.
invisible(gc()) #free up memory

# Packages for plotting
libplt <- c("devtools", "ggh4x", "ggplot2", "latex2exp", "reshape", "scales")
lapply(libplt, require, character.only = TRUE)

# devtools::install_github("stefano-meschiari/latex2exp") # to get ell in TeX

# Load workspace
load("markup_dependence_workspace_80_MC_200_to_1000_n_16_to_128_p_diagonal_only.RData")

## == HELPER FUNCTIONS ##
# Function providing locations of methods to plot
which_methods <- function(met_plt) {
  nummet_plt <- length(met_plt)
  where_met_plt <- numeric(nummet_plt)
  for (thismet in 1:nummet_plt) {
    where_met_plt[thismet] <- which(methods == met_plt[thismet])
  }
  return(where_met_plt)
}

# Function providing locations of designs to plot
which_designs <- function(des_plt) {
  numdes_plt <- length(des_plt)
  where_des_plt <- numeric(numdes_plt)
  for (thisdesign in 1:numdes_plt) {
    where_des_plt[thisdesign] <- which(designs == des_plt[thisdesign])
  }
  return(where_des_plt)
}

# Create labels (via functions) for plots

clab <- function(string) {
  sprintf("c = %s", string) # note: labeller expects a string
}

plab <- function(string) {
  sprintf("p = %s", string) # note: labeller expects a string
}

des_lab <- function(string) {
  des_list <-
    c("Diagonal" = TeX("A$: Diagonal, IID Normal"),
      "NearBand" = TeX("$B$: Near Band, IID Normal"),
      "BlockDiag" = TeX("$C$: Block Diag, IID Normal"),
      "Correlated" = TeX("$D$: Diag, Corr Normal"),
      "HeavyTailed" = TeX("$E$: Diag, Corr Student"),
      "Heteroskedastic" = TeX("$F$: Diag, Hetero Normal"),
      "NearUnity" = TeX("$G$: Near Unity, IID Normal")
    )
  return(des_list[string])
}

met_lab <- function(string) {
  met_list <- c("Lasso" = TeX("Weighted Lasso"),
                "SqrtLasso" = TeX("Sqrt-Lasso"),
                "PostLasso" = TeX("Post-Lasso"),
                "PostSqrtLasso" = TeX("Post-Sqrt-Lasso")
              )
  return(met_list[string])
}

# Function providing string for statistic (included in y-label)
stat_string <- function(stat) {
  switch(stat,
         "max" = "Maximum",
         "mean" = "Average",
         "median" = "Median",
         "q90" = "90th Quantile",
         "q95" = "95th Quantile",
         "var" = "Variance",
         "sd" = "Standard Deviation",
         stop("Invalid statistic"))
}

# Function providing (function for) statistic
stat_fctn <- function(stat) {
  switch(stat,
         "max" = max,
         "mean" = mean,
         "median" = median,
         "q90" = function(x) quantile(x, probs = 0.9),
         "q95" = function(x) quantile(x, probs = 0.95),
         "var" = var,
         "sd" = sd,
         stop("Invalid statistic"))
}

# Color-blind friendly color palette
cb_palette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Shape order
shape_order <- c(0, 1, 2, 5, 4, 3)
# 0=square, 1=circle, 2=triangle, 5=diamond, 4=cross, 3=plus

# Linetype order
linetype_order <- c("solid", "dashed", "dotdash",
                    "dotted", "longdash", "twodash")

# Function to plot error statistics

t_str_stat <- " of $\\max_{i\\in[p]}\\|\\widehat{\\beta}_i-\\beta_{{0}{i}}\\|_{\\ell_2}$"

plot_stat_p_by_method <- function(stat, met_plt, design) {
  where_met_plt <- which_methods(met_plt)
  where_des <- which_designs(design)
  # dimnames(max_ell2_errors)[["p"]] <- plab(dimnames(max_ell2_errors)[["p"]])
  # dimnames(max_ell2_errors)[["design"]] <-
    # des_lab(dimnames(max_ell2_errors)[["design"]])
  # dimnames(max_ell2_errors)[["method"]] <-
    # met_lab(dimnames(max_ell2_errors)[["method"]])
  dimnames(max_ell2_errors)[["c"]] <- clab(dimnames(max_ell2_errors)[["c"]])
  stat_max_ell2 <- apply(max_ell2_errors, c("n", "p", "design", "method", "c"),
                         stat_fctn(stat))
  stat_max_ell2_plt <- stat_max_ell2[, , where_des, where_met_plt, ]
  df <- reshape::melt(stat_max_ell2_plt)
  g <- guide_legend(nrow = 1)
  plt <- ggplot(df, aes(x = n, y = value)) +
    geom_line(aes(linetype = factor(c), color = factor(c)), linewidth = 0.5) +
    geom_point(aes(shape = factor(c), color = factor(c)), size = 2) +
    scale_shape(solid = FALSE) +
    ggh4x::facet_grid2(method ~ p,
                       labeller = labeller(method = met_lab, p = plab),
                       scales = "free_y") +
    labs(
      x = TeX(r"($n$)"),
      y = TeX(paste0("Monte Carlo ", stat_string(stat), t_str_stat))
    ) +
    theme_bw() +
    theme(legend.position = "bottom",
          legend.box.background = element_rect(colour = "black"),
          legend.title = element_blank(),
          legend.text = element_text(size = rel(.75)),
          legend.margin = margin(c(0, 0, 0, 0), unit = "mm"),
          legend.spacing = unit(0, "mm"),
          legend.spacing.x = unit(0, "mm"),
          legend.spacing.y = unit(0, "mm"),
          panel.spacing.y = unit(3, "mm"),
          plot.title = element_text(size = rel(0.5)),
          plot.caption = element_text(hjust = 0, size = rel(0.5))
    ) +
    scale_x_continuous(breaks = seq(200, max(nvec), 200),
                       limits = c(200, max(nvec))) +
    scale_colour_manual(values = cb_palette) +
    scale_shape_manual(values = shape_order) +
    scale_linetype_manual(values = linetype_order) +
    guides(color = g, shape = g, linetype = g)
}

# Plotting statistics of maximum ell_2 errors as function of the score markup
met_plt <- c("Lasso", "PostLasso", "SqrtLasso")
stats_plt <- c("mean")
des_plt <- c("Diagonal")
for (thisdes in seq_along(des_plt)) {
  design <- des_plt[thisdes]
  for (thisstat in seq_along(stats_plt)) {
    plt <- plot_stat_p_by_method(stats_plt[thisstat],
                                 met_plt = met_plt,
                                 design = design)
    file_name <- paste0("markup_dependence_", stats_plt[thisstat], "_error_",
                        design, "_design_", nummc, "_MC_",
                        min(nvec), "_to_", max(nvec), "_n_",
                        min(pvec), "_to_", max(pvec), "_p")
    ggsave(plt, filename = paste0("img/", file_name, ".png"),
           width = 8, height = 6)
  }
}