# Clear
rm(list = ls(all.names = TRUE)) # will clear all (including hidden) objects.
invisible(gc()) #free up memory

# Packages for plotting
libplt <- c("devtools", "ggh4x", "ggplot2", "latex2exp", "reshape", "scales")
lapply(libplt, require, character.only = TRUE)

# devtools::install_github("stefano-meschiari/latex2exp") # to get ell in TeX

# Load workspace
load("simulations_workspace_1000_MC_100_to_1000_n_16_to_128_p_with_num_upd.RData")
# load("simulations_workspace_1000_MC_200_to_1000_n_16_to_128_p_3_h_dot9_rho_with_num_upd_LINUX.RData")
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
  TeX(sprintf("$p = %s$", string)) # note: labeller expects a string - not integer
}

des_lab <- function(string) {
  des_list <-
    c("Diagonal" = TeX("$A$: Diagonal, IID Normal"),
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
  met_list <- c("Lasso" = "Weighted Lasso",
                "AICLasso" = "AIC-Lasso",
                "BICLasso" = "BIC-Lasso",
                "HQICLasso" = "HQIC-Lasso",
                "SqrtLasso" = "Sqrt-Lasso",
                "PostLasso" = "Post-Lasso",
                "PostAICLasso" = "Post-AIC-Lasso",
                "PostBICLasso" = "Post-BIC-Lasso",
                "PostHQICLasso" = "Post-HQIC-Lasso",
                "PostSqrtLasso" = "Post-Sqrt-Lasso")
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

# Function to plot error statistics as method (cols) by design (rows)
t_str_stat <- " of $\\max_{i\\in[p]}\\|\\widehat{\\beta}_i-\\beta_{{0}{i}}\\|_{\\ell_2}$"

plot_stat_p_by_design <- function(stat, met_plt, des_plt) {
  where_des_plt <- which_designs(des_plt)
  where_met_plt <- which_methods(met_plt)
  dimnames(max_ell2_errors)[["p"]] <- plab(dimnames(max_ell2_errors)[["p"]])
  dimnames(max_ell2_errors)[["design"]] <-
    des_lab(dimnames(max_ell2_errors)[["design"]])
  dimnames(max_ell2_errors)[["method"]] <-
    met_lab(dimnames(max_ell2_errors)[["method"]])
  stat_max_ell2 <- apply(max_ell2_errors, c(2, 3, 4, 5), stat_fctn(stat))
  stat_max_ell2_plt <- stat_max_ell2[, , where_des_plt, where_met_plt]
  df <- reshape::melt(stat_max_ell2_plt)
  g <- guide_legend(nrow = 1)
  plt <- ggplot(df, aes(x = n, y = value)) +
    geom_line(aes(linetype = method, color = method), linewidth = 0.5) +
    geom_point(aes(shape = method, color = method), size = 2) +
    scale_shape(solid = FALSE) +
    # geom_hline(yintercept = 1, linetype = "dotted", color = "black") +
    ggh4x::facet_grid2(design ~ p, labeller = label_parsed, scales = "free_y") +
    # facetted_pos_scales(
    #       y = list(
    #         design == des_lab("BlockDiag")[[1]] ~
    #           scale_y_continuous(limits = c(NA, 1.5),
    #                              breaks = scales::breaks_extended(n = 5)),
    #         design != des_lab("BlockDiag")[[1]] ~
    #           scale_y_continuous(breaks = scales::breaks_extended(n = 5))
    #       )
    #     ) +
    labs(
      x = TeX(r"($n$)"),
      y = TeX(paste0("Monte Carlo ", stat_string(stat), t_str_stat))
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
    scale_x_continuous(breaks = seq(200, 1000, 200), limits = c(100, 1000)) +
    scale_colour_manual(values = cb_palette) +
    scale_shape_manual(values = shape_order) +
    scale_linetype_manual(values = linetype_order) +
    guides(color = g, shape = g, linetype = g)
  return(plt)
}

# Function to plot error statistics as method (cols) by design (rows)
t_str_rel_stat <- " of $\\max_{i\\in[p]}\\|\\widehat{\\beta}_i-\\beta_{{0}{i}}\\|_{\\ell_2}$ Relative to that of Weighted Lasso (= 1)"

plot_rel_stat_p_by_design <- function(stat, met_plt, des_plt) {
  where_des_plt <- which_designs(des_plt)
  where_met_plt <- which_methods(met_plt)
  dimnames(max_ell2_errors)[["p"]] <- plab(dimnames(max_ell2_errors)[["p"]])
  stat_max_ell2 <- apply(max_ell2_errors, c(2, 3, 4, 5), stat_fctn(stat))
  rel_stat_max_ell2 <- array(NA, dim(stat_max_ell2))
  dimnames(rel_stat_max_ell2) <- dimnames(stat_max_ell2)
  for (des in designs) {
    for (met in methods) {
      rel_stat_max_ell2[, , des, met] <- stat_max_ell2[, , des, met] /
                                         stat_max_ell2[, , des, "Lasso"]
    }
  }
  dimnames(rel_stat_max_ell2)[["design"]] <-
    des_lab(dimnames(rel_stat_max_ell2)[["design"]])
  dimnames(rel_stat_max_ell2)[["method"]] <-
    met_lab(dimnames(rel_stat_max_ell2)[["method"]])
  rel_stat_max_ell2_plt <- rel_stat_max_ell2[, , where_des_plt, where_met_plt]
  rel_stat_max_ell2_plt <- rel_stat_max_ell2_plt[, , , -1] #  drop Lasso (all 1)
  # stat_max_ell2_plt <- stat_max_ell2[, , where_des_plt, where_met_plt]
  # rel_stat_max_ell2_plt <- array(NA, dim(stat_max_ell2_plt))
  # dimnames(rel_stat_max_ell2_plt) <- dimnames(stat_max_ell2_plt)
  # for (i in seq_along(where_des_plt)) {
  #   for (j in seq_along(where_met_plt)) {
  #     rel_stat_max_ell2_plt[, , i, j] <- stat_max_ell2_plt[, , i, j] /
  #                                        stat_max_ell2_plt[, , i, 1] # <- Lasso
  #   }
  # }
  # rel_stat_max_ell2_plt <- rel_stat_max_ell2_plt[, , , -1] #  drop Lasso (all 1)
  # levels(rel_stat_max_ell2_plt$method) <-
  #   met_lab(levels(rel_stat_max_ell2_plt$method))
  df <- reshape::melt(rel_stat_max_ell2_plt)
  g <- guide_legend()
  plt <- ggplot(df, aes(x = n, y = value)) +
    geom_line(aes(linetype = method, color = method), linewidth = 0.5) +
    geom_point(aes(shape = method, color = method), size = 2) +
    scale_shape(solid = FALSE) +
    geom_hline(yintercept = 1, linetype = "dotted", color = "black") +
    ggh4x::facet_grid2(design ~ p, labeller = label_parsed, scales = "free_y") +
    # facetted_pos_scales(
    #       y = list(
    #         design == "BlockDiag" ~
    #           scale_y_continuous(limits = c(NA, 1.5),
    #                             breaks = scales::breaks_extended(n = 5)),
    #         design != "BlockDiag" ~
    #           scale_y_continuous(breaks = scales::breaks_extended(n = 5))
    #       )
    #     ) +
    labs(
      x = TeX(r"($n$)"),
      y = TeX(paste0("Monte Carlo ", stat_string(stat), t_str_rel_stat))
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
    scale_colour_manual(values = cb_palette[-1]) +
    scale_shape_manual(values = shape_order[-1]) +
    scale_linetype_manual(values = linetype_order[-1]) +
    guides(color = g, shape = g, linetype = g)
  return(plt)
}

plot_stat_numupd_p_by_design <- function(stat, des_plt) {
  where_des_plt <- which_designs(des_plt)
  dimnames(num_upd)[["p"]] <- plab(dimnames(num_upd)[["p"]])
  dimnames(num_upd)[["design"]] <-
    des_lab(dimnames(num_upd)[["design"]])
  dimnames(num_upd)[["method"]] <-
    met_lab(dimnames(num_upd)[["method"]])
  stat_num_upd <- apply(num_upd, c(2, 3, 4, 5), stat_fctn(stat))
  stat_num_upd_plt <- stat_num_upd[, , where_des_plt, ]
  df <- reshape::melt(stat_num_upd_plt)
  g <- guide_legend(nrow = 1)
  plt <- ggplot(df, aes(x = n, y = value)) +
    geom_line(aes(linetype = method, color = method), linewidth = 0.5) +
    geom_point(aes(shape = method, color = method), size = 2) +
    scale_shape(solid = FALSE) +
    facet_grid(design ~ p, labeller = label_parsed) +
    labs(
      x = TeX(r"($n$)"),
      y = TeX(paste0("Monte Carlo ", stat_string(stat),
                     " Number of Loading Updates (Maximum 15)"))
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
    scale_y_continuous(breaks = c(1, 5, 10, 15), limits = c(1, 15)) +
    scale_colour_manual(values = cb_palette) +
    scale_shape_manual(values = shape_order) +
    scale_linetype_manual(values = linetype_order) +
    guides(color = g, shape = g, linetype = g)
}

# 1. Comparing Lasso, Post-Lasso and Sqrt-Lasso
met_plt <- c("Lasso", "PostLasso", "SqrtLasso")
des_plt <- c("Diagonal", "NearBand", "BlockDiag",
             "Correlated", "HeavyTailed",
             "Heteroskedastic", "NearUnity")
stats_plt <- c("mean", "median", "q90", "q95") # statistics to plot

# ... in terms of absolute max ell2 errors
for (thisstat in seq_along(stats_plt)) {
  which_stat <- stats_plt[thisstat]
  plt <- plot_stat_p_by_design(stats_plt[thisstat],
                               met_plt = met_plt,
                               des_plt = des_plt)
  file_name <- paste(stats_plt[thisstat], "_max_ell2_errors_p_by_design_",
                     nummc, "_MC_", min(nvec), "_to_", max(nvec),
                     "_n_", min(pvec), "_to_", max(pvec), "_p", sep = "")
  ggsave(plt, filename = paste0("img/", file_name, ".png"),
         width = 8, height = 12)
}

# ... in terms of relative max ell2 errors
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

# ... and in terms of number of loading updates (<=15; Lasso & Post-Lasso only)
for (thisstat in seq_along(stats_plt)) {
  which_stat <- stats_plt[thisstat]
  plt <- plot_stat_numupd_p_by_design(stats_plt[thisstat], des_plt = des_plt)
  file_name <- paste(stats_plt[thisstat], "_num_upd_p_by_design_",
                     nummc, "_MC_", min(nvec), "_to_", max(nvec),
                     "_n_", min(pvec), "_to_", max(pvec), "_p", sep = "")
  ggsave(plt, filename = paste0("img/", file_name, ".png"),
         width = 8, height = 12)
}