#######################################################X
#----------------------ICTWS 2025----------------------X
#--------------Habitat Selection Workshop--------------X
#----------------Last updated 2025-03-17---------------X
#-----------------SSF Code Walkthrough-----------------X
#######################################################X

# Load packages ----
library(tidyverse)
library(terra)
library(amt)
library(lubridate)
library(circular)

# Load data ----

# ... habitat variables ----
# We'll use the same habitat layers we generated for the HSF walkthrough.
hab <- rast("../HSF Walkthrough/geo/habitat.tif")
names(hab) <- c("forage", "temp", "predator", "cover")

# Cover is a factor, let's code it that way
levels(hab[[4]]) <- data.frame(id = 1:3,
                               cover = c("grassland", "forest", "wetland"))
plot(hab)

# ... location data ----
# Here we'll load the movement trajectory data we simulated under the SSF.
traj <- read.csv("sim_gps.csv")

# ... true values ----
# Here are the true values from the simulation for later.

## Movement-free selection kernel
beta_forage = log(5)/500
beta_pred = log(0.25)/5
beta_temp2 = -1 * log(2)/36
beta_temp = beta_temp2 * -26
beta_forest = log(2)
beta_wetland = log(1/2)

## Selection-free movement kernel
shp <- 4
scl <- 50

# Note that the mean of the gamma is given by shape * scale
shp * scl

kappa <- 1
mu <- 0

# Model fitting ----
# Now we can use our simulated data and fit a model.
head(traj)

# ... format data ----
# Random seed since we're generating random steps
set.seed(20240429 * 2)

issf_dat <- traj %>%
  # Format time column as a POSIX object
  mutate(t_ = as.POSIXct(t_)) %>% 
  # Format as a 'track_xyt' object
  make_track(x_, y_, t_, crs = 32612) %>% 
  # or use steps_by_burst()
  steps() %>% 
  # Default is to use gamma and von Mises
  random_steps(n_control = 100) %>% 
  # This is a useful function if you want to add time of day to each point
  time_of_day(where = "both") %>% 
  # Extract habitat covariates at the end of each step
  extract_covariates(hab, where = "end") %>% 
  # Convert tod_* to factor
  mutate(tod_start_ = factor(tod_start_),
         tod_end_ = factor(tod_end_)) %>% 
  mutate(log_sl_ = log(sl_),
         cos_ta_ = cos(ta_)) %>% 
  filter(!is.na(ta_))

# Take a look at what we've got.
print(issf_dat, n = 5, width = Inf)

# Notice that it has attributes, including the CRS, the tentative movement
# parameters, and the number of control steps.
str(issf_dat, 2)

# We have helper functions to extract the tentative movement parameters
sl_distr(issf_dat)
ta_distr(issf_dat)

# For later, we will need just the observed steps, so we can
# store them separately now.
obs <- issf_dat %>% 
  filter(case_)

# ... fit iSSF ----
# Fit the model

issf <- issf_dat %>% 
  # Make your own dummy variables -- R often makes too many levels with
  # categorical variables for clogit models.
  # Remember grassland is our reference level
  mutate(forest = as.numeric(cover == "forest"),
         wetland = as.numeric(cover == "wetland")) %>%
  # Fit the model
  fit_issf(case_ ~ 
             # Habitat
             forage + temp + I(temp^2) + predator + 
             # All the cover terms
             forest + wetland + 
             # Movement
             log_sl_ + sl_ +
             cos_ta_ +
             # Strata (steps)
             strata(step_id_), model = TRUE)

# Model summary
summary(issf)

# Let's look at the structure of the 'issf' object.
str(issf, 1)

# Notice that it is a list with 4 elements at the top level.
#   - $model: the actual fitted model
#   - $sl_: the tentative step-length distribution
#   - $ta_: the tentative turn-angle distribution
#   - $more: (currently empty) a placeholder for additional information

# Let's make a figure to see how well we did estimating the parameters. Note
# that more generic functions are available to us if we call the model object
# directly rather than the "fit_clogit" object.
b <- coef(issf)
ci <- confint(issf$model)

res <- as.data.frame(cbind(b, ci))
res$truth <- c(beta_forage,
               beta_temp,
               beta_temp2,
               beta_pred,
               beta_forest,
               beta_wetland,
               rep(NA, 3))
names(res) <- c("est", "lwr", "upr", "truth")

# Plot
(beta_fig <- res %>% 
    # Keep only habitat betas
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
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)))

# Zoom in to see forage and temp^2
beta_fig +
  coord_cartesian(ylim = c(-0.025, 0.005))

# We did very well!

# ... update movement distributions ----
# We did well with the movement-free habitat kernel. How about with the
# selection-free movement kernel?

# We need to use the movement betas from the fitted model to update the 
# tentative step-length and turn-angle distributions to the estimated 
# selection-free distributions.

# Recall that the formulas are available in Appendix C of Fieberg et al. 2021.
# https://conservancy.umn.edu/bitstream/handle/11299/218272/AppC_iSSA_movement.html

# We also have functions in `amt` to make this a bit easier for you.
?update_sl_distr

# Note that you cannot simply use 'update_sl_distr()' or 'update_ta_distr()' if
# you have interactions with your movement parameters. You need to pass the
# fully updated betas to the correct updating function. For example:
?update_gamma

# Update step-length distribution
sl_upd <- update_sl_distr(issf)

# Update turn-angle distribution
ta_upd <- update_ta_distr(issf)

# How did we do?
data.frame(parm = c("shape", "scale", "kappa"),
           est = c(sl_upd$params$shape,
                   sl_upd$params$scale,
                   ta_upd$params$kappa),
           truth = c(shp,
                     scl,
                     kappa))

# We did very well!

# Figures ----

# ... landcover ----
# Let's use RSS to express selection for landcover types
# We can use 'amt::log_rss()' just like we did for HSFs.


x1_lc <- data.frame(forage = mean(obs$forage),
                    temp = mean(obs$temp),
                    predator = mean(obs$predator),
                    cover = factor(c("grassland", "forest", "wetland"),
                                   levels = c("grassland", "forest", "wetland")),
                    sl_ = 100,
                    log_sl_ = log(100),
                    cos_ta_ = 1) %>% 
  # Make your own dummy variables
  # Remember grassland is our reference level
  mutate(forest = as.numeric(cover == "forest"),
         wetland = as.numeric(cover == "wetland"))

x2_lc <- data.frame(forage = mean(obs$forage),
                    temp = mean(obs$temp),
                    predator = mean(obs$predator),
                    # We don't have to compare to the reference level
                    # We can choose any level
                    cover = factor("grassland",
                                   # Need to pass the same factor levels as fitted model
                                   levels = c("grassland", "forest", "wetland")),
                    tod_start_ = factor("day", levels = c("day", "night")),
                    tod_end_ = factor("day", levels = c("day", "night")),
                    sl_ = 100,
                    log_sl_ = log(100),
                    cos_ta_ = 1) %>% 
  # Make your own dummy variables
  # Remember grassland is our reference level
  mutate(forest = as.numeric(cover == "forest"),
         wetland = as.numeric(cover == "wetland"))

# Calculate log-RSS
log_rss_lc <- log_rss(issf, x1 = x1_lc, x2 = x2_lc, ci = "se")

# Plot
(lr_lc_plot <- log_rss_lc$df %>% 
    ggplot(aes(x = cover_x1, y = log_rss)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    geom_point(position = position_dodge(0.25)) +
    geom_errorbar(aes(ymin = lwr, ymax = upr), 
                  position = position_dodge(0.25), width = 0.1) +
    xlab("Land Cover") +
    ylab("log-RSS vs Grassland") +
    theme_bw())

# ... temperature ----
# Let's use RSS to express selection for selection for temperature
# We can use 'amt::log_rss()' just like we did for HSFs.

x1_temp <- data.frame(forage = mean(obs$forage),
                      temp = seq(from = 3, to = 20, length.out = 100),
                      predator = mean(obs$predator),
                      cover = factor(c("grassland"),
                                     # Need all the levels from the original data
                                     levels = c("grassland", "forest", "wetland")),
                      sl_ = 100,
                      log_sl_ = log(100),
                      cos_ta_ = 1) %>% 
  # Make your own dummy variables
  # Remember grassland is our reference level
  mutate(forest = as.numeric(cover == "forest"),
         wetland = as.numeric(cover == "wetland"))

x2_temp <- data.frame(forage = mean(obs$forage),
                      # We can compare this to any value. Often I choose the
                      # mean, but we can choose 10 degrees here
                      temp = 5,
                      predator = mean(obs$predator),
                      cover = factor("grassland",
                                     levels = c("grassland", "forest", "wetland")),
                      tod_start_ = factor("day", levels = c("day", "night")),
                      tod_end_ = factor("day", levels = c("day", "night")),
                      sl_ = 100,
                      log_sl_ = log(100),
                      cos_ta_ = 1) %>% 
  # Make your own dummy variables
  # Remember grassland is our reference level
  mutate(forest = as.numeric(cover == "forest"),
         wetland = as.numeric(cover == "wetland"))

# Calculate log-RSS
log_rss_temp <- log_rss(issf, x1 = x1_temp, x2 = x2_temp, ci = "se")

# Plot
(lr_temp_plot <- log_rss_temp$df %>% 
    ggplot(aes(x = temp_x1, y = log_rss)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    geom_ribbon(aes(ymin = lwr, ymax = upr), 
                linewidth = 0.1, linetype = "dashed",
                color = "gray30", fill = "gray90") +
    geom_line(linewidth = 1) +
    xlab("Temperature (C)") +
    ylab("log-RSS vs 5 degrees") +
    theme_bw())

# ... step-length distribution ----
# We also want to plot the step-length distribution.

data.frame(sl = seq(0.1, 1000, length.out = 100)) %>% 
  mutate(
    shp = sl_upd$params$shape,
    scl = sl_upd$params$scale,
    y = dgamma(sl, shape = shp, scale = scl)
  ) %>% 
  ggplot(aes(x = sl, y = y)) +
  geom_line(linewidth = 1) +
  xlab("Step Length (m)") +
  ylab("Probability Density") +
  theme_bw()


# ... turn-angle distribution ----
# We also want to plot the turn angle distribution.

data.frame(ta = seq(-pi, pi, length.out = 100)) %>% 
  rowwise() %>% # need this because circular::dvonmises is not vectorized
  mutate(
    k = ta_upd$params$kappa,
    y = suppressWarnings(dvonmises(ta, mu = 0, kappa = k))
  ) %>% 
  ggplot(aes(x = ta, y = y)) +
  geom_line(linewidth = 1) +
  xlab("Turn Angle (radians)") +
  ylab("Probability Density") +
  theme_bw()

