# Differences between day and night?

# Load packages ----
library(tidyverse)
library(amt)
library(terra)

# 1. Load data ----
# Location data
dat <- read_csv("data/coyote_cougar.csv") %>% 
  # Subset to just cougar F53
  # filter(id == "F53") %>% 
  # Subset to just cougar F64
  # filter(id == "F64") %>% 
  # Subset to just coyote C028
  filter(id == "C028") %>%
  arrange(t_)

# Set the timezone
tz(dat$t_) <- "US/Mountain"

# Habitat data as stack
hab <- rast("data/coyote_cougar_habitat.tif")
names(hab) <- c("elevation", "trees", "biomass", "dist_to_road")

# Format as track_xyt
trk <- dat %>% 
  make_track(x_, y_, t_, crs = 32612)

stp <- track_resample(trk, rate = hours(round(summarize_sampling_rate(trk)[["median"]])), 
                      tolerance = minutes(30)) %>% 
  filter_min_n_burst(min_n = 3) %>% 
  steps_by_burst()

# Format data
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


# ToD model
m <- issf_dat %>% 
  fit_issf(case_ ~ 
             # Habitat selection main effects
             elevation + trees + biomass + 
             # Movement main effects
             sl_ + log_sl_ + cos_ta_ +
             # Movement interactions
             sl_ : tod_start_ + log_sl_ : tod_start_ + 
             # Don't forget the strata
             strata(step_id_),
           # And include model = TRUE so we can use 'log_rss()' later
           model = TRUE)

summary(m) 


dd <- issf_dat %>% 
  filter(case_)

ggplot(dd, aes(x = sl_, color = tod_start_)) +
  geom_density()

day_sl <- update_gamma(sl_distr(issf_dat),
                       beta_sl = coef(m)[["sl_"]],
                       beta_log_sl = coef(m)[["log_sl_"]])

night_sl <- update_gamma(sl_distr(issf_dat),
                         beta_sl = coef(m)[["sl_:tod_start_night"]],
                         beta_log_sl = coef(m)[["log_sl_:tod_start_night"]])

expand.grid(sl = seq(1, 8000, length.out = 100),
            tod = c("day", "night")) %>% 
  rowwise() %>% 
  mutate(shp = get(paste0(tod, "_sl"))$params$shape,
         scl = get(paste0(tod, "_sl"))$params$scale,
         y = dgamma(sl, shape = shp, scale = scl)) %>% 
  ggplot(aes(x = sl, y = y, color = tod)) +
  geom_line()
