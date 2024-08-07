# Clear
rm(list = ls(all.names = TRUE)) # will clear all (including hidden) objects.
invisible(gc()) #free up memory

# Packages for plotting
libplt <- c("devtools", "ggplot2", "reshape")

# Auto-installer (checks if installed - if not, installs and loads)
using <- function(...) {
  libs <- unlist(list(...))
  req <- unlist(lapply(libs, require, character.only = TRUE))
  need <- libs[req == FALSE]
  if (length(need) > 0) {
    install.packages(pkgs = need, repos = "https://cran.us.r-project.org")
    lapply(need, require, character.only = TRUE)
  }
}
using(libplt)
# devtools::install_github("stefano-meschiari/latex2exp")
library("latex2exp")

# Load workspace
load("simulations_workspace_10_MC_100_to_500_n_4_to_16_p_3_h.RData")

# Specify statistics to plot and methods to include
stats_plt <- c("mean", "max") # statistics to plot (fed to apply)

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

# Create labels (via functions) for plots
plab <- function(string) {
  sprintf("p = %s", string) # note: labeller expects a string - not integer
}

# Function providing string for statistic (included in y-label)
stat_string <- function(stat) {
  switch(stat,
         "mean" = "Average",
         "max" = "Maximum",
         stop("Invalid statistic"))
}

# Function taking statistic and methods to plot, returning ggplot
t_str <- " of $\\max_{i\\in[p]}\\|\\widehat{\\beta}_i-\\beta_{{0}{i}}\\|_{\\ell_2}$"
plot_stat_p_by_design <- function(stat, met_plt) {
  where_met_plt <- which_methods(met_plt)
  stat_max_ell2_errors <- apply(max_ell2_errors, c(2, 3, 4, 5), stat)
  df <- reshape::melt(stat_max_ell2_errors[, , , method = where_met_plt])
  plt <- ggplot(df, aes(x = n, y = value, color = as.factor(method))) +
    geom_line(linewidth = 0.5) +
    geom_point(size = 2) +
    facet_grid(design ~ p, labeller = labeller(p = plab)) +
    labs(
      x = TeX(r"($n$)"),
      y = TeX(paste0("Monte Carlo ", stat_string(stat), t_str))
    ) +
    theme_bw() +
    theme(legend.justification = c(1, 1), legend.position = c(1, 1),
            legend.key = element_rect(colour = "transparent", fill = NA),
            legend.box.background = element_rect(colour = "black"),
            legend.title = element_blank(),
            legend.text = element_text(size = rel(0.5)),
            legend.margin = margin(c(0, 0, 0, 0), unit = "mm"),
            legend.spacing = unit(0, "mm"),
            legend.spacing.x = unit(0, "mm"),
            legend.spacing.y = unit(0, "mm"),
            panel.spacing.y = unit(3, "mm"),
            plot.title = element_text(size = rel(0.5)),
            plot.caption = element_text(hjust = 0, size = rel(0.5))
    )
  return(plt)
}

# 1. Comparisons w/o refitting (i.e. no "Post" in name)
met_plt <- c("Lasso", "BICLasso", "SqrtLasso")
# met_plt <- c("Lasso", "AICLasso", "BICLasso", "SqrtLasso")
# met_plt <- methods[!grepl("Post", methods)]

for (thisstat in seq_along(stats_plt)) {
  plt <- plot_stat_p_by_design(stats_plt[thisstat], met_plt)
  file_name <- paste("plain_methods_",
                     stats_plt[thisstat], "_max_ell2_errors_p_by_design_",
                     nummc, "_MC_", min(nvec), "_to_", max(nvec),
                     "_n_", min(pvec), "_to_", max(pvec), "_p", sep = "")
  ggsave(plt, filename = paste0("img/", file_name, ".png"),
         width = 8, height = 12)
}

# 2. Comparisons w/ refitting (i.e. "Post" in name)
met_plt <- c("PostLasso", "PostBICLasso", "PostSqrtLasso")
# met_plt <- c("PostLasso", "PostAICLasso", "PostBICLasso", "PostSqrtLasso")
# met_plt <- methods[grepl("Post", methods)]
for (thisstat in seq_along(stats_plt)) {
  plt <- plot_stat_p_by_design(stats_plt[thisstat], met_plt)
  file_name <- paste("post_methods_",
                     stats_plt[thisstat], "_max_ell2_errors_p_by_design_",
                     nummc, "_MC_", min(nvec), "_to_", max(nvec),
                     "_n_", min(pvec), "_to_", max(pvec), "_p", sep = "")
  ggsave(plt, filename = paste0("img/", file_name, ".png"),
         width = 8, height = 12)
}