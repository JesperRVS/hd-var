# Clear
rm(list = ls(all.names = TRUE)) # will clear all (including hidden) objects.
invisible(gc()) #free up memory

# Packages for plotting
libplt <- c("devtools", "ggh4x", "ggplot2", "latex2exp", "reshape", "scales")
lapply(libplt, require, character.only = TRUE)
# devtools::install_github("stefano-meschiari/latex2exp") # to get ell in TeX

# Load workspace
# load("simulations_workspace_10_MC_100_to_500_n_4_to_16_p_3_h.RData")
# load("simulations_workspace_80_MC_100_to_500_n_16_to_128_p_3_h_LINUX_with_num_upd.RData")
load("simulations_workspace_250_MC_200_to_1000_n_16_to_128_p_3_h_LINUX_with_num_upd.RData")

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
plab <- function(string) {
  sprintf("p = %s", string) # note: labeller expects a string - not integer
}

des_lab <- function(string) {
  des_list <- list("Diagonal" = "A: Diagonal, IID Gaussian",
                   "BlockDiag" = "B: Block diag., IID Gaussian",
                   "NearBand" = "C: Near band, IID Gaussian",
                   "Correlated" = "D: Diag., Correlated Gaussian",
                   "Heteroskedastic" = "E: Diag., Heterosked. Gaussian",
                   "HeavyTailed" = "F: Diag., IID Student's t")
  return(des_list[string])
}

met_lab <- function(string) {
  met_list <- list("Lasso" = "Lasso",
                   "AICLasso" = "AIC Lasso",
                   "BICLasso" = "BIC Lasso",
                   "HQICLasso" = "HQIC Lasso",
                   "SqrtLasso" = "Square-Root Lasso",
                   "PostLasso" = "Post Lasso",
                   "PostAICLasso" = "Post AIC Lasso",
                   "PostBICLasso" = "Post BIC Lasso",
                   "PostHQICLasso" = "Post HQIC Lasso",
                   "PostSqrtLasso" = "Post Square-Root Lasso")
  return(met_list[string])
}

# Specify statistics to plot and methods to include
stats_plt <- c("mean", "median", "q90") # statistics to plot
# stats_plt <- c("max", "mean", "median", "q90", "sd") # statistics to plot

# Function providing string for statistic (included in y-label)
stat_string <- function(stat) {
  switch(stat,
         "max" = "Maximum",
         "mean" = "Average",
         "median" = "Median",
         "q90" = "90th Quantile",
         "var" = "Variance",
         "sd" = "Standard Deviation",
         stop("Invalid statistic"))
}

# Function providing function for statistic
stat_fctn <- function(stat) {
  switch(stat,
         "max" = max,
         "mean" = mean,
         "median" = median,
         "q90" = function(x) quantile(x, probs = 0.9),
         "var" = var,
         "sd" = sd,
         stop("Invalid statistic"))
}

# Color-blind friendly color palette 
cb_palette <- c("#E69F00", "#56B4E9", "#000000", "#009E73",
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Shape order
shape_order <- c(0, 1, 2, 5, 4) # square, circle, triangle, diamond, cross

# Linetype order
linetype_order <- c("solid", "dashed", "dotted", "dotdash", "longdash")

# Function to plot error statistics as method (cols) by design (rows)
t_str <- " of $\\max_{i\\in[p]}\\|\\widehat{\\beta}_i-\\beta_{{0}{i}}\\|_{\\ell_2}$"
plot_stat_method_by_design <- function(stat, met_plt = methods, des_plt = designs) {
  where_des_plt <- which_designs(des_plt)
  where_met_plt <- which_methods(met_plt)
  stat_max_ell2_errors <- apply(max_ell2_errors, c(2, 3, 4, 5), stat_fctn(stat))
  df <- reshape::melt(stat_max_ell2_errors[, , design = where_des_plt, method = where_met_plt])
  g <- guide_legend(title = TeX("$p =$"))
  plt <- ggplot(df, aes(x = n, y = value)) +
    geom_line(aes(linetype = as.factor(p), color = as.factor(p)), linewidth = 0.5) +
    geom_point(aes(shape = as.factor(p), color = as.factor(p)), size = 2) +
    scale_shape(solid = FALSE) +
    ggh4x::facet_grid2(design ~ method,
                       labeller = labeller(p = plab,
                                           method = met_lab,
                                           design = des_lab),
                       scales = "free_y") +
    facetted_pos_scales(
      y = list(
        design == "BlockDiag" ~
          scale_y_continuous(limits = c(NA, 0.55),
                             breaks = scales::breaks_extended(n = 5)),
        design != "BlockDiag" ~
          scale_y_continuous(breaks = scales::breaks_extended(n = 5))
      )
    ) +
    labs(
      x = TeX(r"($n$)"),
      y = TeX(paste0("Monte Carlo ", stat_string(stat), t_str))
    ) +
    theme_bw() +
    theme(legend.justification = c(1, 1), legend.position = c(1, 1),
            legend.key = element_rect(colour = "transparent", fill = NA),
            legend.box.background = element_rect(colour = "black"),
            legend.title = element_text(size = rel(.75), hjust = 0.5),
            legend.text = element_text(size = rel(.75)),
            legend.margin = margin(c(0, 0, 0, 0), unit = "mm"),
            legend.spacing = unit(0, "mm"),
            legend.spacing.x = unit(0, "mm"),
            legend.spacing.y = unit(0, "mm"),
            panel.spacing.y = unit(3, "mm"),
            plot.title = element_text(size = rel(0.5)),
            plot.caption = element_text(hjust = 0, size = rel(0.5))
    ) +
    scale_colour_manual(values = cb_palette) +
    scale_shape_manual(values = shape_order) +
    scale_linetype_manual(values = linetype_order) +
    guides(color = g, shape = g, linetype = g)
  return(plt)
}

# 1. Comparing Lasso and SqrtLasso methods
lasso_sqrtlasso <- c("Lasso", "SqrtLasso", "PostLasso", "PostSqrtLasso")
des_plt <- c("Diagonal", "BlockDiag", "NearBand",
             "Correlated", "Heteroskedastic", "HeavyTailed")

for (thisstat in seq_along(stats_plt)) {
  which_stat <- stats_plt[thisstat]
  plt <- plot_stat_method_by_design(stats_plt[thisstat],
                                    met_plt = lasso_sqrtlasso,
                                    des_plt = des_plt)
  file_name <- paste("all_lasso_sqrtlasso_methods_",
                     stats_plt[thisstat], "_max_ell2_errors_method_by_design_",
                     nummc, "_MC_", min(nvec), "_to_", max(nvec),
                     "_n_", min(pvec), "_to_", max(pvec), "_p", sep = "")
  ggsave(plt, filename = paste0("img/", file_name, ".png"),
         width = 8, height = 12)
}

# 2. Comparing Lasso and BIC methods
lasso_biclasso <- c("Lasso", "BICLasso", "PostLasso", "PostBICLasso")
des_plt <- c("Diagonal", "BlockDiag", "NearBand",
             "Correlated", "Heteroskedastic", "HeavyTailed")

for (thisstat in seq_along(stats_plt)) {
  which_stat <- stats_plt[thisstat]
  plt <- plot_stat_method_by_design(stats_plt[thisstat],
                                    met_plt = lasso_biclasso,
                                    des_plt = des_plt)
  file_name <- paste("all_lasso_bic_methods_",
                     stats_plt[thisstat], "_max_ell2_errors_method_by_design_",
                     nummc, "_MC_", min(nvec), "_to_", max(nvec),
                     "_n_", min(pvec), "_to_", max(pvec), "_p", sep = "")
  ggsave(plt, filename = paste0("img/", file_name, ".png"),
         width = 8, height = 12)
}