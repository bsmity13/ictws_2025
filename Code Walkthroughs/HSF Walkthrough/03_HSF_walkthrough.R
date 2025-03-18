#######################################################X
#----------------------ICTWS 2025----------------------X
#--------------Habitat Selection Workshop--------------X
#----------------Last updated 2025-03-17---------------X
#-----------------HSF Code Walkthrough-----------------X
#######################################################X

# Load packages ----
library(tidyverse)
library(terra)
library(amt)

# Load data ----

# ... habitat raster ----

# I've already randomly generated 4 variables:
#   - forage (g/m^2) [resource]
#   - mean annual temperature (°C) [condition]
#   - predator density (predators/100 km^2) [risk]
#   - land cover type (grassland = 1, forest = 2, wetland = 3) [condition]

# If you want to see how I did that, check out the script "01_sim_habitat.R". 
# We can load those here as a SpatRast.

hab <- rast("geo/habitat.tif")
names(hab) <- c("forage", "temp", "predator", "cover")

# Cover is a factor, let's code it that way
levels(hab[[4]]) <- data.frame(id = 1:3,
                               cover = c("grassland", "forest", "wetland"))
hab$cover
plot(hab)

# ... location data ----

# I simulated these data under the generating model of the HSF.

# If you want to see how I did that, check out the script "02_sim_data.R". 
# We can load those here from CSV file.
gps <- read.csv("sim_gps.csv")

head(gps)

# Fitting an HSF with 'amt' ----
# Let's fit an HSF.
mod_dat <- gps %>% 
  # First, we turn it into a 'track_xy' object
  make_track(x, y, crs = 32612) %>% 
  # This is technically sampling within a 100% MCP by default,
  # but that is practically the extent of our raster.
  random_points(n = nrow(gps) * 100) %>% 
  # We attach all the covariates from our raster
  extract_covariates(hab) %>% 
  # Assign large weights to available points
  mutate(weight = ifelse(case_, 1, 1e5))

# Fit a model
hsf <- glm(case_ ~ forage + temp + I(temp^2) + predator + cover,
           family = binomial(), weights = weight, data = mod_dat)

# Summary
summary(hsf)

# Compare our fitted coefficients to true values
b <- coef(hsf)
ci <- confint(hsf)

# These are the true values. Again, see "02_sim_data.R" for these.
beta_forage <- log(5)/500
beta_pred <- log(0.25)/5
beta_temp2 <- -1 * log(2)/36
beta_temp <- beta_temp2 * -26
beta_forest <- log(2)
beta_wetland <- log(1/2)

# Combine in data.frame for plotting
res <- as.data.frame(cbind(b, ci))
res$truth <- c(NA,
               beta_forage,
               beta_temp,
               beta_temp2,
               beta_pred,
               beta_forest,
               beta_wetland)
names(res) <- c("est", "lwr", "upr", "truth")

# Plot
(beta_fig <- res %>% 
    filter(!is.na(truth)) %>%
    mutate(name = row.names(.)) %>% 
    ggplot(aes(x = name, y = est)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_point(size = 1.25) +
    geom_errorbar(aes(ymin = lwr, ymax = upr),
                  width = 0.2) +
    geom_point(aes(y = truth), color = "red") +
    xlab(expression(beta)) +
    ylab(NULL) +
    theme_bw())

# Zoom in to see forage and temp^2
beta_fig +
  coord_cartesian(ylim = c(-0.025, 0.005))

# Looks great!

# Model interpretation ----
# Using RSS is the best way to understand what's going on biologically.

# We can calculate log-RSS using a function from `amt`.

?log_rss

# *Importantly*, the data.frame 'x2' should almost always have only 1-row
# to avoid any confusion with R's vector-recycling rules.

# ... scenarios ----
# One useful thing we can do with RSS is to generate some biologically 
# interesting scenarios and ask how many times more points we expect
# in one habitat (real or hypothetical) vs another habitat.

# For example, we might be interested in the tradeoff between foraging
# and predation risk. We could create a scenario with high forage and low
# predator density vs. a habitat with low forage and high predator density.
# How many times more points do we expect if:
#   - x1: 
#       - forage: 750 g/m^2
#       - predator density: 3 predators/100 km^2
#   - x2: 
#       - forage: 250 g/m^2
#       - predator density: 8 predators/100 km^2

# Note that we didn't specify temperature or land cover. We still have to
# pass these variables to 'log_rss()', but *IT DOES NOT MATTER* what values
# we choose for x1 and x2, as long as they are the same (because our model
# has no interactions).

# Define x1
x1 <- data.frame(forage = 750,
                 predator = 3,
                 temp = 15,
                 # Note factors need the same levels as the original data
                 cover = factor("grassland", 
                                levels = c("grassland", "forest", "wetland")))

# Define x2
x2 <- data.frame(forage = 250,
                 predator = 8,
                 temp = 15,
                 cover = factor("grassland", 
                                levels = levels(mod_dat$cover)))

# Calculate log-RSS and 95% CI
scenario_rss <- log_rss(hsf, x1, x2, ci = "se")

# Examine the structure of the resulting object
str(scenario_rss)

# The first object, "df", is typically what you'll want to work with.
scenario_rss$df

# The data.frame shows us the values at x1 and the resulting log-RSS.
# Exponentiate to answer our question.
scenario_rss$df %>% 
  dplyr::select(log_rss:upr) %>% 
  mutate(across(everything(), exp))

# We expect almost 20x as many points if forage is high and predation is low
# vs. if forage is low and predation is high. We can also see the 95% CI
# for that estimate.

# You can see how generating interesting scenarios with combinations of 
# different habitats can be interesting.

# Another common approach is to use RSS to look across a range of values
# for one habitat axis, with all others held constant. Let's do that for each
# of our habitat dimensions.

# *Note* that I am re-using the names 'x1' and 'x2' in each example. This
# works fine if you always run the script top-to-bottom, but it gets quite
# dangerous if you are skipping around and running the script interactively.
# Be careful if you decide not to give each x1 and x2 unique names!

# ... forage ----
# x1 is a sequence across the whole range of forage
# All other variables held constant
x1 <- data.frame(forage = seq(0, 1000, length.out = 100),
                 temp = mean(mod_dat$temp),
                 predator = mean(mod_dat$predator),
                 cover = factor("grassland", 
                                levels = c("grassland", "forest", "wetland")))

# x2 is still a single row
x2 <- data.frame(forage = mean(mod_dat$forage) + sd(mod_dat$forage),
                 temp = mean(mod_dat$temp),
                 predator = mean(mod_dat$predator),
                 cover = factor("grassland", 
                                levels = c("grassland", "forest", "wetland")))

# Calculate log-RSS
forage_rss <- log_rss(hsf, x1, x2, ci = "se")

# Plot RSS
(forage_plot <- forage_rss$df %>% 
    mutate(rss = exp(log_rss),
           exp_lwr = exp(lwr),
           exp_upr = exp(upr)) %>% 
    ggplot(aes(x = forage_x1, y = rss, ymin = exp_lwr, ymax = exp_upr)) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    geom_ribbon(linewidth = 0.5, linetype = "dashed", fill = "gray70", color = "black") +
    geom_line(linewidth = 1) +
    xlab(expression("Forage " * (g/m^2))) +
    ylab("RSS") +
    theme_bw())

# ... temperature ----
x1 <- data.frame(forage = mean(mod_dat$forage),
                 temp = seq(2, 20, length.out = 100),
                 predator = mean(mod_dat$predator),
                 cover = factor("grassland", 
                                levels = c("grassland", "forest", "wetland")))

x2 <- data.frame(forage = mean(mod_dat$forage),
                 temp = 13,
                 predator = mean(mod_dat$predator),
                 cover = factor("grassland", 
                                levels = c("grassland", "forest", "wetland")))

temp_rss <- log_rss(hsf, x1, x2, ci = "se")

(temp_plot <- temp_rss$df %>% 
    mutate(rss = exp(log_rss),
           exp_lwr = exp(lwr),
           exp_upr = exp(upr)) %>% 
    ggplot(aes(x = temp_x1, y = rss, ymin = exp_lwr, ymax = exp_upr)) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    geom_ribbon(linewidth = 0.5, linetype = "dashed", fill = "gray70", color = "black") +
    geom_line(linewidth = 1) +
    xlab("Temperature (°C)") +
    ylab("RSS") +
    theme_bw())

# ... predator density ----
x1 <- data.frame(forage = mean(mod_dat$forage),
                 temp = mean(mod_dat$temp),
                 predator = seq(0, 12, length.out = 100),
                 cover = factor("grassland", 
                                levels = c("grassland", "forest", "wetland")))

x2 <- data.frame(forage = mean(mod_dat$forage),
                 temp = mean(mod_dat$temp),
                 predator = mean(mod_dat$predator),
                 cover = factor("grassland", 
                                levels = c("grassland", "forest", "wetland")))

pred_rss <- log_rss(hsf, x1, x2, ci = "se")

(pred_plot <- pred_rss$df %>% 
    mutate(rss = exp(log_rss),
           exp_lwr = exp(lwr),
           exp_upr = exp(upr)) %>% 
    ggplot(aes(x = predator_x1, y = rss, ymin = exp_lwr, ymax = exp_upr)) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    geom_ribbon(linewidth = 0.5, linetype = "dashed", fill = "gray70", color = "black") +
    geom_line(linewidth = 1) +
    xlab(expression("Predator Density " * (pred/km^2))) +
    ylab("RSS") +
    theme_bw())

# ... land cover ----
x1 <- data.frame(forage = mean(mod_dat$forage),
                 temp = mean(mod_dat$temp),
                 predator = mean(mod_dat$predator),
                 cover = factor(c("grassland", "forest", "wetland"), 
                                levels = c("grassland", "forest", "wetland")))

x2 <- data.frame(forage = mean(mod_dat$forage),
                 temp = mean(mod_dat$temp),
                 predator = mean(mod_dat$predator),
                 cover = factor("wetland", 
                                levels = c("grassland", "forest", "wetland")))

cover_rss <- log_rss(hsf, x1, x2, ci = "se")

(cover_plot <- cover_rss$df %>% 
    mutate(rss = exp(log_rss),
           exp_lwr = exp(lwr),
           exp_upr = exp(upr)) %>% 
    ggplot(aes(x = cover_x1, y = rss, ymin = exp_lwr, ymax = exp_upr)) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    geom_errorbar(width = 0.25, color = "black") +
    geom_point(size = 2) +
    xlab("Land Cover") +
    ylab("RSS") +
    theme_bw())

# Mapping the HSF ----
# We can use the HSF to predict across space now.
# We may want to plot w(x). We can use R's generic 'predict()' function
# to do this, but recall, we don't want the intercept. Let's call the 
# linear prediction with the intercept y(x).
map_dat <- as.data.frame(hab, xy = TRUE) %>% 
  mutate(y_x = predict(hsf, newdata = .),
         # Subtract off the intercept to get g(x)
         g = y_x - coef(hsf)[1],
         # Exponentiate to get w(x)
         w = exp(g))

# Map w(x)
(map_w <- ggplot(map_dat, aes(x = x, y = y, fill = w)) +
  geom_raster() +
  coord_sf(crs = 32612) +
  scale_fill_viridis_c(name = expression(w(x))) +
  xlab(NULL) +
  ylab(NULL) +
  theme_bw()
)

# If we add the point, we can gain some intuition for what the generating
# model is:
map_w +
  geom_point(aes(x = x, y = y), data = gps, inherit.aes = FALSE,
             color = "white", size = 0.5)
# Points are sprinkled randomly, but more intensely where w(x) is higher

# Interpreting w(x) is difficult here, so it might be more intuitive to plot
# RSS vs mean conditions.

x1 <- as.data.frame(hab, xy = TRUE)

x2 <- data.frame(forage = mean(mod_dat$forage),
                 temp = mean(mod_dat$temp),
                 predator = mean(mod_dat$predator),
                 cover = factor("grassland", 
                                levels = c("grassland", "forest", "wetland")))

map_rss <- log_rss(hsf, x1, x2)

(map <- map_rss$df %>% 
    mutate(rss = exp(log_rss)) %>% 
    ggplot(aes(x = x_x1, y = y_x1, fill = rss)) +
    geom_raster() +
    coord_sf(crs = 32612) +
    scale_fill_viridis_c(name = "RSS") +
    xlab(NULL) +
    ylab(NULL) +
    theme_bw())

# The new map looks the same because the RSS is linearly proportional to w(x),
# but the values are easier to understand intuitively now.

# The problem with this is that we can't really tell what's going on below 1.
# log-RSS would be better for that

(map <- map_rss$df %>% 
    ggplot(aes(x = x_x1, y = y_x1, fill = log_rss)) +
    geom_raster() +
    coord_sf(crs = 32612) +
    scale_fill_gradient2(name = "log-RSS",
                        low = "navyblue",
                        mid = "white",
                        high = "firebrick",
                        midpoint = 0) +
    xlab(NULL) +
    ylab(NULL) +
    theme_bw())

# By using a color palette that distinguishes between negative, positive, and
# 0, we can at a glance see if a pixel has habitats that would be used more than
# expected, less than expected, or directly in proportion to their availability.