# Clear
rm(list = ls(all.names = TRUE)) # will clear all (including hidden) objects.
invisible(gc()) #free up memory

n <- 10000
p <- 4
sigma_eps <- 0.1
h <- 0.1
nburn <- 100

source("simulations/simData.R")

y0ton <- sim_data_h(n = n, p = p, sigma_eps = sigma_eps,
                    h = h, nburn = nburn)

colMeans(y0ton)
colMeans(y0ton^2)

library("ggplot2")
library("reshape")

# I have p time series in y0ton with 1 + n observations each.
# I want to put these into a suitable data frame for ggplot2
# and then plot them all on the same plot.

# First, I need to reshape the data so that I have a column for each time series.
# I can do this with the reshape package.
y0ton_melt <- melt(y0ton, id.vars = "t")

# Now I can plot the data.
ggplot(y0ton_melt, aes(x = X1, y = value, color = as.factor(X2))) +
  geom_line() +
  geom_point(size = 2) +
  labs(title = "Time series",
       x = "Time",
       y = "Value") +
  theme_minimal()

# Show me just the first one.
ggplot(y0ton_melt[y0ton_melt$X2 == 1, ], aes(x = X1, y = value)) +
  geom_line() +
  geom_point(size = 2) +
  labs(title = "Time series",
       x = "Time",
       y = "Value") +
  theme_minimal()

# calculate the sample variance of each time series

