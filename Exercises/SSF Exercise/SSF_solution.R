#######################################################X
#----------------------ICTWS 2025----------------------X
#--------------Habitat Selection Workshop--------------X
#----------------Last updated 2025-03-06---------------X
#-----------------SSF Exercise Solution----------------X
#######################################################X

# Using UT cougar data

# Load packages ----
library(tidyverse)
library(amt)
library(terra)

# 1. Load data ----
# Location data
dat <- read_csv("data/coyote_cougar.csv") %>% 
  # Subset to just cougar F53
  filter(id == "F53")

# Set the timezone
tz(dat$t_) <- "US/Mountain"

# Habitat data as stack
hab <- rast("data/coyote_cougar_habitat.tif")
names(hab) <- c("elevation", "trees", "biomass", "dist_to_road")

# Format as track_xyt
trk <- dat %>% 
  make_track(x_, y_, t_, crs = 32612)

# Do we have regular steps?
summarize_sampling_rate(trk)

# The median sampling rate is 4h, but the track isn't perfectly regular.
# Let's separate into bursts of 4h steps, with a tolerance of 15 minutes.
# Only keep bursts that have at least 3 relocations. Then we can turn those
# locations into steps.

stp <- track_resample(trk, rate = hours(4), tolerance = minutes(15)) %>% 
  filter_min_n_burst(min_n = 3) %>% 
  steps_by_burst()

# What do our step durations look like now?
hist(as.numeric(stp$dt_))

# 2. Consider habitat and movement variables ----

## Habitats
plot(hab)

# We have:
#   - elevation (m)
#   - tree cover (%)
#   - biomass of annual and perennial grasses & forbs (kg/ha)
#   - distance to road (m)

# I would categorize these as:
#   - elevation: risk during winter (too much snow at high elevation?)
#               deep snow would also probably slow down movement (short step lengths)
#   - tree cover: resource during warm days, provides resting shade
#   - biomass: resource -- maybe high biomass of forage attracts prey species
#   - distance to road: no effect during winter because roads are not open

# Notice the temporal dynamics in some of my a priori predictions. Some of
# my predictions depend on time of day (e.g., trees as shade resource) and some
# depend on season (e.g., only need shade during summer). For simplicity here,
# I am just going to fit a model for a single season.

## Movement variables

# We are going to use a gamma distribution to model step lengths, so we need
# to include 'sl_' and 'log(sl_)' in our model.

# We are going to use a von Mises distribution to model turn angles, so we 
# need to include 'cos(ta_)' in our model.

# How might these interact? Maybe high elevations also have a lot of snow
# during winter. So I expect that at high elevation, step lengths will be
# shorter during winter.

# There may also be daily patterns. Cougars primarily hunt at night, so perhaps
# steps are longer at night than during the day.

# 3. Fit iSSFs ----
issf_dat <- stp %>% 
  # Keep the defaults for generating available steps (gamma and von Mises distrs)
  # Increase the number of random steps
  random_steps(n_control = 50) %>% 
  # Attach habitat variables
  extract_covariates(hab) %>% 
  # Add temporal covariates
  time_of_day(where = "both") %>% 
  # Add additional movement covariates
  mutate(log_sl_ = log(sl_),
         cos_ta_ = cos(ta_))

# Subset for winter
winter_dat <- issf_dat %>% 
  filter(month(t1_) %in% c(11, 12, 1:3))

# Winter model
winter_issf <- winter_dat %>% 
  fit_issf(case_ ~ 
             # Habitat selection main effects
             elevation + trees + biomass + 
             # Habitat interactions
             trees : tod_end_ + 
             # Movement main effects
             sl_ + log_sl_ + cos_ta_ +
             # Movement interactions
             sl_ : tod_start_ + log_sl_ : tod_start_ + cos_ta_ : tod_start_ +
             sl_ : elevation + log_sl_ : elevation +
             # Don't forget the strata
             strata(step_id_),
           # And include model = TRUE so we can use 'log_rss()' later
           model = TRUE)

summary(winter_issf) 

# 4. RSS Figures ----
# We found a significant effect of trees on habitat selection and a significant
# interaction with time of day.
# Let's create a figure to examine our hypothesis that our cougar prefers
# tree cover during the day, but not at night.

# day
x1_day <- data.frame(elevation = 2500,
                           trees = seq(0, 100, length.out = 25),
                           biomass = 2000,
                           dist_to_road = 1500,
                           tod_start_ = factor("day", 
                                               levels = c("day", "night")),
                           tod_end_ = factor("day", 
                                             levels = c("day", "night")), 
                           sl_ = 900,
                           log_sl_ = log(900),
                           cos_ta_ = 1)

x2_day <- data.frame(elevation = 2500,
                           trees = 0,
                           biomass = 2000,
                           dist_to_road = 1500,
                           tod_start_ = factor("day", 
                                               levels = c("day", "night")),
                           tod_end_ = factor("day", 
                                             levels = c("day", "night")), 
                           sl_ = 900,
                           log_sl_ = log(900),
                           cos_ta_ = 1)

rss_day <- log_rss(winter_issf, x1_day, x2_day, ci = "se")

# night
x1_night <- data.frame(elevation = 2500,
                             trees = seq(0, 100, length.out = 25),
                             biomass = 2000,
                             dist_to_road = 1500,
                             tod_start_ = factor("night", 
                                                 levels = c("day", "night")),
                             tod_end_ = factor("night", 
                                               levels = c("day", "night")), 
                             sl_ = 900,
                             log_sl_ = log(900),
                             cos_ta_ = 1)

x2_night <- data.frame(elevation = 2500,
                             trees = 0,
                             biomass = 2000,
                             dist_to_road = 1500,
                             tod_start_ = factor("night", 
                                                 levels = c("day", "night")),
                             tod_end_ = factor("night", 
                                               levels = c("day", "night")), 
                             sl_ = 900,
                             log_sl_ = log(900),
                             cos_ta_ = 1)

rss_night <- log_rss(winter_issf, x1_night, x2_night, ci = "se")

# Now grab all the data.frames from the 'log_rss' objects and combine
fig_dat <- bind_rows("Day" = rss_day$df,
                     "Night" = rss_night$df,
                     .id = "time") %>% 
  # Convert log-RSS to RSS
  mutate(rss = exp(log_rss),
         rss_lwr = exp(lwr),
         rss_upr = exp(upr))

# Check
head(fig_dat)

# Plot
fig_dat %>% 
  ggplot(aes(x = trees_x1, y = log_rss, ymin = lwr, ymax = upr,
             color = time, fill = time)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1) +
  geom_ribbon(linetype = "dashed", alpha = 0.2, linewidth = 1) +
  geom_line(linewidth = 1) +
  xlab("Tree Cover (%)") +
  ylab("log-RSS") +
  theme_bw()

# You can see that the 95% confidence interval during night overlaps with 0.
# However, during the day, we see a clear preference for higher tree cover.

# 5. Movement figures ----
# Now let's make a figure to look at our hypothesis that movement slows down
# at high elevation during winter seasons. 

# Note that none of the interaction terms for the step-length distribution
# are very large, so the answer is most likely no. But let's go through the
# exercise anyway.
summary(winter_issf)

# Get betas from our winter-season model
b <- coef(winter_issf)

# Update step-length distributions for different elevations.
elev <- seq(2000, 3500, length.out = 3)


gamma_elev <- expand.grid(elev = elev,
                          tod_start_ = c("day", "night")) %>% 
  # Recall, that the betas for sl_ and for log_sl_ are now functions of elevation
  # and time of day (day is the reference).
  mutate(b_sl_ = b[["sl_"]] + b[["elevation:sl_"]] * elev + 
           b[["sl_:tod_start_night"]] * (tod_start_ == "night"), 
         b_log_sl_ = b[["log_sl_"]] + b[["elevation:log_sl_"]] * elev + 
           b[["log_sl_:tod_start_night"]] * (tod_start_ == "night"))

# Now that we have the betas for sl and log_sl as a function of elev and tod,
# we can update our gamma distribution.
# (notice that this is vectorized)
gamma_elev_distr <- update_gamma(sl_distr(winter_issf), 
                                 gamma_elev$b_sl_,
                                 gamma_elev$b_log_sl_)

# And now add the shape and scale to our data.frame
gamma_elev$shp <- gamma_elev_distr$params$shape
gamma_elev$scl <- gamma_elev_distr$params$scale

# Check
gamma_elev

# Plot mean step length
gamma_elev %>% 
  mutate(mean_sl = shp * scl) %>% 
  ggplot(aes(x = elev, y = mean_sl, color = tod_start_)) +
  geom_line(linewidth = 1) +
  xlab("Elevation (m)") +
  ylab("Mean Step Length (m)") +
  scale_color_discrete(name = "Time") +
  theme_bw()

# Mean step length decreases with elevation regardless of time of day

# Now we want to calculate the probability density under the gamma with those
# parameters for a range of step lengths

step_fig_dat <- gamma_elev %>% 
  # Convert to tibble for nested data.frame
  as_tibble() %>%
  # Add elevation as list column
  mutate(elev_list = list(tibble(sl = seq(1, 500, length.out = 100)))) %>% 
  # Unnest
  unnest(cols = elev_list) %>% 
  # Calculate probability density
  mutate(dens = dgamma(sl, shape = shp, scale = scl))


# Plot distribution
ggplot(step_fig_dat, aes(x = sl, y = dens, col = elev, group = elev)) +
  facet_wrap(~ tod_start_) +
  geom_line() +
  xlab("Step Length (m)") +
  ylab("Proability Density") +
  theme_bw()

# These distributions are all concentrated near 0 and have very long tails,
# something we might expect from a large ambush predator like a cougar.

# From this figure, it is not clear that there is a large difference in step
# length distributions by elevation or by time of day.

