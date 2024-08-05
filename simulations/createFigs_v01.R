# Clear
rm(list = ls(all.names = TRUE)) # will clear all (including hidden) objects.
invisible(gc()) #free up memory

# Load workspace
load("simulations/simulations_workspace_10_MC_100_to_500_n_4_to_16_p.RData")

# average over the MC repetitions
avg_max_ell2_errors <- apply(max_ell2_errors, c(2, 3, 4, 5), mean)

# plot the average ell_2 errors as a function of n on the first axis, and a line for each p
# for the methods Lasso, BICLasso, and SqrtLasso (horizontal facets for the methods)
# and each design as a facet (vertical facets for the designs)
met_plt <- c("Lasso", "PostLasso", "BICLasso", "SqrtLasso")
nummet_plt <- length(met_plt)
where_met_plt <- numeric(nummet_plt)
for (thismet in 1:nummet_plt) {
  where_met_plt[thismet] <- which(methods == met_plt[thismet])
}
df <- reshape::melt(avg_max_ell2_errors[, , , method = where_met_plt])
library("ggplot2")
p_mean <- ggplot(df, aes(x = n, y = value, color = as.factor(p))) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 2) +
  facet_grid(rows = vars(design), cols = vars(method)) +
  ylab("Average Maximum Rowwise Estimation Error") +
  theme_bw()
p_mean

# create a file name where I include the number of MC repetitions into the string
# Also include the range of p and n
file_name <- paste("mean_max_ell2_errors_", nummc, "_MC_", min(nvec), "_to_", max(nvec),
                   "_n_", min(pvec), "_to_", max(pvec), "_p", sep = "")
ggsave(p_mean, filename = paste0("simulations/img/", file_name, ".pdf"),
       width = 8, height = 12)


