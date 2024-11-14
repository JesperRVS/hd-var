# Clear
rm(list = ls(all.names = TRUE)) # will clear all (including hidden) objects.
invisible(gc()) #free up memory

# Packages for plotting
libplt <- c("devtools", "ggh4x", "ggplot2", "latex2exp", "reshape", "scales")
lapply(libplt, require, character.only = TRUE)

# devtools::install_github("stefano-meschiari/latex2exp") # to get ell in TeX

# Load workspace
load("simulations_workspace_80_MC_100_to_500_n_4_to_64_p_3_h_dot9_rho_with_num_upd_LINUX.RData")
# load("simulations_workspace_100_MC_100_to_500_n_4_to_64_p_3_h_dot9_rho_with_num_upd.RData")
# load("simulations_workspace_50_MC_100_to_500_n_4_to_64_p_3_h_dot9_rho_with_num_upd.Rdata")
# load("simulations_workspace_80_MC_100_to_500_n_16_to_128_p_3_h_LINUX_with_num_upd.RData")
# load("simulations_workspace_250_MC_200_to_1000_n_16_to_128_p_3_h_LINUX_with_num_upd.RData")

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
                   "Heteroskedastic_y" = "E: Diag., y-Hetero. Gaussian",
                   "Heteroskedastic_eta" = "E: Diag., eta-Hetero. Gaussian",
                   "HeavyTailed" = "F: Diag., IID Student's t")
  return(des_list[string])
}

# met_lab <- function(string) {
#   met_list <- list("Lasso" = "Lasso",
#                    "AICLasso" = "AIC Lasso",
#                    "BICLasso" = "BIC Lasso",
#                    "HQICLasso" = "HQIC Lasso",
#                    "SqrtLasso" = "Square-Root Lasso",
#                    "PostLasso" = "Post Lasso",
#                    "PostAICLasso" = "Post AIC Lasso",
#                    "PostBICLasso" = "Post BIC Lasso",
#                    "PostHQICLasso" = "Post HQIC Lasso",
#                    "PostSqrtLasso" = "Post Square-Root Lasso")
#   return(met_list[string])
# }

# # Rewrite met_lab using a switch instead of a list
# met_lab <- function(string) {
#   switch(string,
#          "Lasso" = "Lasso",
#          "AICLasso" = "AIC Lasso",
#          "BICLasso" = "BIC Lasso",
#          "HQICLasso" = "HQIC Lasso",
#          "SqrtLasso" = "Square-Root Lasso",
#          "PostLasso" = "Post Lasso",
#          "PostAICLasso" = "Post AIC Lasso",
#          "PostBICLasso" = "Post BIC Lasso",
#          "PostHQICLasso" = "Post HQIC Lasso",
#          "PostSqrtLasso" = "Post Square-Root Lasso",
#          stop("Invalid method"))
# }
# Rewrite met_lab to allow for vector of strings
met_lab <- function(string) {
  met_list <- c("Lasso" = "Lasso",
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
# Rewrite met_lab to only return the formatted strings




# Specify statistics to plot and methods to include
stats_plt <- c("mean", "median", "q90", "q95") # statistics to plot
# stats_plt <- c("max", "mean", "median", "q90", "sd") # statistics to plot

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

# Function providing function for statistic
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
cb_palette <- c("#E69F00", "#56B4E9", "#000000", "#009E73",
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Shape order
shape_order <- c(0, 1, 2, 5, 4, 3)
# 0=square, 1=circle, 2=triangle, 5=diamond, 4=cross, 3=plus

# Linetype order
linetype_order <- c("solid", "dashed", "dotdash", 
                    "dotted", "longdash", "twodash")

# Function to plot error statistics as method (cols) by design (rows)
t_str <- " of $\\max_{i\\in[p]}\\|\\widehat{\\beta}_i-\\beta_{{0}{i}}\\|_{\\ell_2}$ Relative to that of Lasso"

plot_rel_stat_p_by_design <- function(stat, met_plt, des_plt) {
  where_des_plt <- which_designs(des_plt)
  where_met_plt <- which_methods(met_plt)
  stat_max_ell2 <- apply(max_ell2_errors, c(2, 3, 4, 5), stat_fctn(stat))
  stat_max_ell2_plt <- stat_max_ell2[, , where_des_plt, where_met_plt]
  rel_stat_max_ell2_plt <- array(NA, dim(stat_max_ell2_plt))
  dimnames(rel_stat_max_ell2_plt) <- dimnames(stat_max_ell2_plt)
  for (i in seq_along(where_des_plt)) {
    for (j in seq_along(where_met_plt)) {
      rel_stat_max_ell2_plt[, , i, j] <- stat_max_ell2_plt[, , i, j] /
                                         stat_max_ell2_plt[, , i, 1] # <- Lasso
    }
  }
  rel_stat_max_ell2_plt <- rel_stat_max_ell2_plt[, , , -1] #  drop Lasso (all 1)
  df <- reshape::melt(rel_stat_max_ell2_plt)
  g <- guide_legend()
  plt <- ggplot(df, aes(x = n, y = value)) +
    geom_line(aes(linetype = method, color = method), linewidth = 0.5) +
    geom_point(aes(shape = method, color = method), size = 2) +
    scale_shape(solid = FALSE) +
    geom_hline(yintercept = 1, linetype = "dotted", color = "black") +
    ggh4x::facet_grid2(design ~ p,
                       labeller = labeller(p = plab, design = des_lab,
                                           method = met_lab),
                       scales = "free_y") +
    # facetted_pos_scales(
    #       y = list(
    #         design == "BlockDiag" ~
    #           scale_y_continuous(limits = c(NA, 1.4),
    #                             breaks = scales::breaks_extended(n = 5)),
    #         design != "BlockDiag" ~
    #           scale_y_continuous(breaks = scales::breaks_extended(n = 5))
    #       )
    #     ) +
    labs(
      x = TeX(r"($n$)"),
      y = TeX(paste0("Monte Carlo ", stat_string(stat), t_str))
    ) +
    theme_bw() +
    theme(#legend.justification = c(1, 1), legend.position = c(1, 1),
          legend.position = "bottom",
          # legend.key = element_rect(colour = "transparent", fill = NA),
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
    scale_colour_manual(values = cb_palette) +
    scale_shape_manual(values = shape_order) +
    scale_linetype_manual(values = linetype_order) +
    guides(color = g, shape = g, linetype = g)
  return(plt)
}

# 1. Comparing Lasso, Sqrt-Lasso, Post-Lasso and BIC-Lasso
# met_plt <- c("Lasso", "PostLasso", "SqrtLasso", "PostSqrtLasso", "BICLasso", "PostBICLasso")
# des_plt <- c("Diagonal", "Correlated")
# met_plt <- c("Lasso", "PostLasso", "SqrtLasso", "BICLasso")
met_plt <- c("Lasso", "PostLasso", "SqrtLasso", "BICLasso", "PostBICLasso")
des_plt <- c("Diagonal", "BlockDiag", "NearBand", "Correlated",
             "Heteroskedastic_y", "Heteroskedastic_eta", "HeavyTailed")
# met_plt <- c("Lasso", "PostLasso", "SqrtLasso", "PostSqrtLasso", "BICLasso")
# des_plt <- c("Diagonal", "BlockDiag", "NearBand",
#              "Correlated", "Heteroskedastic", "HeavyTailed")
stats_plt <- c("mean", "median", "q90", "q95") # statistics to plot
for (thisstat in seq_along(stats_plt)) {
  which_stat <- stats_plt[thisstat]
  plt <- plot_rel_stat_p_by_design(stats_plt[thisstat],
                                   met_plt = met_plt,
                                   des_plt = des_plt)
  file_name <- paste("rel_",
                     stats_plt[thisstat], "_max_ell2_errors_p_by_design_",
                     nummc, "_MC_", min(nvec), "_to_", max(nvec),
                     "_n_", min(pvec), "_to_", max(pvec), "_p", sep = "")
  ggsave(plt, filename = paste0("img/", file_name, ".png"),
         width = 8, height = 12)
}
