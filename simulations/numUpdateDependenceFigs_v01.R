# Clear
rm(list = ls(all.names = TRUE)) # will clear all (including hidden) objects.
invisible(gc()) #free up memory

# Packages for plotting
libplt <- c("devtools", "ggh4x", "ggplot2", "latex2exp", "reshape", "scales")
lapply(libplt, require, character.only = TRUE)

# Load workspace
load("numupd_dependence_workspace_1000_MC_500_n_16_p_15_k_Diagonal_design.RData")

met_lab <- function(string) {
  met_list <- c("Lasso" = TeX("Weighted Lasso"),
                "PostLasso" = TeX("Post-Lasso")
              )
  return(met_list[string])
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

# # (a) Monte Carlo mean of relative change in penalty loadings for Lasso
# # and Post-Lasso
# rel_change_ups <- array(NA, dim = c(nummc, k, 2))
# dimnames(rel_change_ups) <- list(mc = 1:nummc, k = 1:k,
#                                  method = c("Lasso", "PostLasso"))
# rel_change_ups[, , "Lasso"] <- lasso_rel_change_loadings
# rel_change_ups[, , "PostLasso"] <- postl_rel_change_loadings
# dimnames(rel_change_ups)[["method"]] <-
#   met_lab(dimnames(rel_change_ups)[["method"]])
# mean_rel_change_ups <- apply(rel_change_ups, c("k", "method"), mean)

# df <- reshape::melt(mean_rel_change_ups)
# g <- guide_legend(ncol = 1)
# plt_rel_change_ups <- ggplot(df, aes(x = k, y = value, color = method)) +
#   geom_line(aes(group = method, linetype = method), size = 1) +
#   geom_point(aes(shape = method), size = 3) +
#   scale_shape(solid = FALSE) +
#   scale_color_manual(values = cb_palette[1:2]) +
#   scale_linetype_manual(values = linetype_order[1:2]) +
#   labs(
#     x = "Update (k)", y = "Relative Change in Penalty Loadings (k vs. k-1)",
#     # title = "Relative Change: Loadings"
#   ) +
#   theme_bw() +
#   theme(legend.justification = c(1, 1), legend.position = c(1, 1),
#         legend.box.background = element_rect(colour = "black"),
#         legend.title = element_blank(),
#         legend.text = element_text(size = 12),
#         axis.title = element_text(size = 14),
#         axis.text = element_text(size = 12),
#         plot.title = element_text(size = 16, face = "bold"),
#         plot.subtitle = element_text(size = 14, face = "bold")
#   ) +
#   scale_x_continuous(breaks = c(1, 5, 10, 15), limits = c(1, 15)) + 
#   scale_colour_manual(values = cb_palette) +
#   scale_shape_manual(values = shape_order) +
#   scale_linetype_manual(values = linetype_order) +
#   guides(color = g, shape = g, linetype = g)

# # (b) Monte Carlo mean of relative change in estimates for Lasso and Post-Lasso
# rel_change_est <- array(NA, dim = c(nummc, k, 2))
# dimnames(rel_change_est) <- list(mc = 1:nummc, k = 1:k,
#                                  method = c("Lasso", "PostLasso"))
# rel_change_est[, , "Lasso"] <- lasso_rel_change_estimates
# rel_change_est[, , "PostLasso"] <- postl_rel_change_estimates
# dimnames(rel_change_est)[["method"]] <-
#   met_lab(dimnames(rel_change_est)[["method"]])
# mean_rel_change_est <- apply(rel_change_est, c("k", "method"), mean)

# df <- reshape::melt(mean_rel_change_est)
# # g <- guide_legend(ncol = 1)
# plt_rel_change_est <- ggplot(df, aes(x = k, y = value, color = method)) +
#   geom_line(aes(group = method, linetype = method), size = 1) +
#   geom_point(aes(shape = method), size = 3) +
#   scale_shape(solid = FALSE) +
#   scale_color_manual(values = cb_palette[1:2]) +
#   scale_linetype_manual(values = linetype_order[1:2]) +
#   labs(
#     x = "Update (k)", y = "Relative Change in Coefficient Estimates (k vs. k-1)",
#   ) +
#   theme_bw() +
#   theme(legend.position = "none",
#         axis.title = element_text(size = 14),
#         axis.text = element_text(size = 12),
#         plot.title = element_text(size = 16, face = "bold"),
#         plot.subtitle = element_text(size = 14, face = "bold")
#   ) +
#   scale_x_continuous(breaks = c(1, 5, 10, 15), limits = c(1, 15)) + 
#   scale_colour_manual(values = cb_palette) +
#   scale_shape_manual(values = shape_order) +
#   scale_linetype_manual(values = linetype_order)

# # (c) Monte Carlo mean of the 
# rel_est_error <- array(NA, dim = c(nummc, k, 2))
# dimnames(rel_est_error) <- list(mc = 1:nummc, k = 1:k,
#                                 method = c("Lasso", "PostLasso"))
# rel_est_error[, , "Lasso"] <- lasso_max_ell2_errors_refi /
#   lasso_max_ell2_errors_init
# rel_est_error[, , "PostLasso"] <- postl_max_ell2_errors_refi /
#   postl_max_ell2_errors_init
# dimnames(rel_est_error)[["method"]] <-
#   met_lab(dimnames(rel_est_error)[["method"]])
# mean_rel_est_error <- apply(rel_est_error, c("k", "method"), mean)

# df <- reshape::melt(mean_rel_est_error)
# g <- guide_legend(ncol = 1)
# plt_rel_est_error <- ggplot(df, aes(x = k, y = value, color = method)) +
#   geom_line(aes(group = method, linetype = method), size = 1) +
#   geom_point(aes(shape = method), size = 3) +
#   scale_shape(solid = FALSE) +
#   scale_color_manual(values = cb_palette[1:2]) +
#   scale_linetype_manual(values = linetype_order[1:2]) +
#   labs(
#     x = "Update (k)", y = "Relative Error in Coefficient Estimates (k vs. init)",
#     # title = "Relative Error: Estimates"
#   ) +
#   theme_bw() +
#   theme(legend.position = "none",
#         axis.title = element_text(size = 14),
#         axis.text = element_text(size = 12),
#         plot.title = element_text(size = 16, face = "bold"),
#         plot.subtitle = element_text(size = 14, face = "bold")
#   ) +
#   scale_x_continuous(breaks = c(1, 5, 10, 15), limits = c(1, 15)) + 
#   scale_colour_manual(values = cb_palette) +
#   scale_shape_manual(values = shape_order) +
#   scale_linetype_manual(values = linetype_order)

# # put the two plots side by side
# library(gridExtra)
# plts <- arrangeGrob(plt_rel_change_ups, plt_rel_change_est,
#                     plt_rel_est_error, ncol = 3)
# file_name <- paste0("numupd_dependence_",
#                     nummc, "_MC_",
#                     n, "_n_", p, "_p_", k, "_k_",
#                     design, "_design")
# ggsave(plts, filename = paste0("img/", file_name, ".png"),
#        width = 8, height = 6)

# Similar three plots but only for Lasso

rel_change_ups <- array(NA, dim = c(nummc, k, 1))
dimnames(rel_change_ups) <- list(mc = 1:nummc, k = 1:k,
                                 method = "Lasso")
rel_change_ups[, , "Lasso"] <- lasso_rel_change_loadings
mean_rel_change_ups <- apply(rel_change_ups, c("k", "method"), mean)

df <- reshape::melt(mean_rel_change_ups)
plt_rel_change_ups <- ggplot(df, aes(x = k, y = value, color = method)) +
  geom_line(aes(group = method, linetype = method), size = 1) +
  geom_point(aes(shape = method), size = 3) +
  scale_shape(solid = FALSE) +
  scale_color_manual(values = cb_palette[1]) +
  scale_linetype_manual(values = linetype_order[1]) +
  labs(
    x = "Update (k)", y = "Relative Change in Loadings (k vs. k-1)",
  ) +
  theme_bw() +
  theme(legend.position = "none") +
  scale_x_continuous(breaks = c(1, 5, 10, 15), limits = c(1, 10)) + 
  scale_colour_manual(values = cb_palette) +
  scale_shape_manual(values = shape_order) +
  scale_linetype_manual(values = linetype_order)

rel_change_est <- array(NA, dim = c(nummc, k, 1))
dimnames(rel_change_est) <- list(mc = 1:nummc, k = 1:k,
                                 method = "Lasso")
rel_change_est[, , "Lasso"] <- lasso_rel_change_estimates
mean_rel_change_est <- apply(rel_change_est, c("k", "method"), mean)

df <- reshape::melt(mean_rel_change_est)

plt_rel_change_est <- ggplot(df, aes(x = k, y = value, color = method)) +
  geom_line(aes(group = method, linetype = method), size = 1) +
  geom_point(aes(shape = method), size = 3) +
  scale_shape(solid = FALSE) +
  scale_color_manual(values = cb_palette[1]) +
  scale_linetype_manual(values = linetype_order[1]) +
  labs(
    x = "Update (k)", y = "Relative Change in Estimates (k vs. k-1)",
  ) +
  theme_bw() +
  theme(legend.position = "none") +
  scale_x_continuous(breaks = c(1, 5, 10, 15), limits = c(1, 10)) + 
  scale_colour_manual(values = cb_palette) +
  scale_shape_manual(values = shape_order) +
  scale_linetype_manual(values = linetype_order)

rel_est_error <- array(NA, dim = c(nummc, k, 1))
dimnames(rel_est_error) <- list(mc = 1:nummc, k = 1:k,
                                method = "Lasso")
rel_est_error[, , "Lasso"] <- lasso_max_ell2_errors_refi /
    lasso_max_ell2_errors_init
mean_rel_est_error <- apply(rel_est_error, c("k", "method"), mean)

df <- reshape::melt(mean_rel_est_error)
g <- guide_legend(ncol = 1)
plt_rel_est_error <- ggplot(df, aes(x = k, y = value, color = method)) +
  geom_line(aes(group = method, linetype = method), size = 1) +
  geom_point(aes(shape = method), size = 3) +
  scale_shape(solid = FALSE) +
  scale_color_manual(values = cb_palette[1]) +
  scale_linetype_manual(values = linetype_order[1]) +
  labs(
    x = "Update (k)", y = "Relative Estimation Error (k vs. init)",
    # title = "Relative Error: Estimates"
  ) +
  theme_bw() +
  theme(legend.position = "none") +
  scale_x_continuous(breaks = c(1, 5, 10, 15), limits = c(1, 10)) +
  scale_colour_manual(values = cb_palette) +
  scale_shape_manual(values = shape_order) +
  scale_linetype_manual(values = linetype_order)


plts_lasso <- arrangeGrob(plt_rel_change_ups, plt_rel_change_est,
                    plt_rel_est_error, ncol = 3)
file_name <- paste0("numupd_dependence_",
                    nummc, "_MC_",
                    n, "_n_", p, "_p_", k, "_k_",
                    design, "_design", "_Lasso_only")
ggsave(plts_lasso, filename = paste0("img/", file_name, ".png"),
       width = 8, height = 4)
