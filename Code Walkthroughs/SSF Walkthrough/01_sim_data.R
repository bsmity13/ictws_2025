#######################################################X
#----------------------ICTWS 2025----------------------X
#--------------Habitat Selection Workshop--------------X
#----------------Last updated 2025-03-17---------------X
#-----------------SSF Code Walkthrough-----------------X
#######################################################X

# Simulate data under iSSF for code walkthrough

# Load packages ----
library(tidyverse)
library(terra)
library(amt)
library(lubridate)
library(circular)

# Generate example data ----
# We need to generate some example telemetry data to fit our iSSF.
# Generating points from an iSSF is trickier than generating points
# from an HSF because of (1) the dependence on the start location  
# (i.e., it is a time series) and (2) the addition of the movement model.

# ... habitat variables ----
# We'll use the same habitat layers we generated for the HSF walkthrough.
hab <- rast("../HSF Walkthrough/geo/habitat.tif")
names(hab) <- c("forage", "temp", "predator", "cover")

# Cover is a factor, let's code it that way
levels(hab[[4]]) <- data.frame(id = 1:3,
                               cover = c("grassland", "forest", "wetland"))
plot(hab)

# ... movement-free habitat kernel ----
# Now we have to choose betas that will give us our movement-free habitat
# selection kernel.

# We'll start with the same values for beta from the HSF walkthrough. 
# We still interpret each beta as the log-RSS for a one-unit change in 
# the covariate. But now we need to recognize that the interpretation
# is the probability of *a step* ending in that habitat. The estimated
# values also depend on scale, since they are conditional on the movement
# model and the start location.

beta_forage = log(5)/500
beta_pred = log(0.25)/5
beta_temp2 = -1 * log(2)/36
beta_temp = beta_temp2 * -26
beta_forest = log(2)
beta_wetland = log(1/2)

# Since our landscape is small we can calculate the movement-free habitat
# selection kernel for the entire landscape. 
hab_vals <- as.data.frame(hab) %>% 
  # Calculate g(x) 
  mutate(g =
           # forage
           beta_forage * forage +
           # two terms for temperature
           beta_temp * temp +
           beta_temp2 * temp^2 +
           # predator density
           beta_pred * predator +
           # landcover
           beta_forest * (cover == "forest") +
           beta_wetland * (cover == "wetland")) %>% 
  # Exponentiate to get w(x)
  mutate(w = exp(g)) %>% 
  # Normalize so they sum to 1
  mutate(w_prime = w/sum(w))

# Place in raster
# Use hab[[1]] as template
hab_kern <- hab[[1]]
# Insert values (not using the normalized values here)
values(hab_kern) <- hab_vals$w
names(hab_kern) <- "hab"
# And a vector of values
hab_kern_vals <- values(hab_kern)

# Plot
plot(hab_kern)

# ... selection-free movement kernel ----

# Now we have to specify the movement part of the model. 

# Let's assume our steps come from a gamma distribution and our turn
# angles come from a von Mises distribution.

# Also note that our landscape is only 2 km x 2 km, so our steps should be
# quite short to stay within the landscape.

# Let's choose these shape and scale parameters.
shp <- 4
scl <- 50

# Note that the mean of the gamma is given by shape * scale
shp * scl

# Plot the distributions
data.frame(sl = seq(0.1, 1000, length.out = 100)) %>% 
  mutate(y = dgamma(sl, shape = shp, scale = scl)) %>% 
  ggplot(aes(x = sl, y = y)) +
  geom_line(linewidth = 1) +
  xlab("Step Length (m)") +
  ylab("Probability Density") +
  theme_bw()

# We also need to specify our turn angle distribution.

kappa <- 1
mu <- 0

# Plot the distributions
data.frame(ta = seq(-pi, pi, length.out = 100)) %>% 
  rowwise() %>% 
  mutate(y = circular::dvonmises(ta, mu = mu, kappa = kappa)) %>% 
  ggplot(aes(x = ta, y = y)) +
  geom_line(linewidth = 1) +
  xlab("Turn Angle (radians)") +
  ylab("Probability Density") +
  coord_cartesian(ylim = c(0, NA)) +
  theme_bw()

# ... times ----

# Let's setup a data.frame with timestamps to hold our simulated data

dat <- data.frame(id = "A01",
                  time = seq(ymd_hms("2021-02-05 0:00:00", tz = "US/Mountain"),
                             ymd_hms("2021-05-05 00:00:00", tz = "US/Mountain"),
                             by = "1 hour"),
                  x = NA,
                  y = NA)

# We'll start our animal right in the middle of our map
dat$x[1] <- mean(ext(hab))[1]
dat$y[1] <- mean(ext(hab))[2]

# Check
plot(hab[[1]])
points(dat$x, dat$y, pch = 16, col = "red")

#
dat <- dat %>% 
  # Rearrange/rename columns
  dplyr::select(id, x1 = x, y1 = y, t1 = time) %>% 
  # Convert to steps
  mutate(x2 = lead(x1),
         y2 = lead(y1),
         t2 = lead(t1)) %>% 
  filter(!is.na(t2))

# Let's define our first step's endpoint as being 50m directly
# north to get us started moving. This is also the second step's start point.
dat$y2[1] <- dat$y1[2] <- dat$y1[1] + 50
dat$x2[1] <- dat$x1[2] <- dat$x1[1] + 0

# Lastly, let's add the absolute angle of the first step
dat$abs_angle <- NA
dat$abs_angle[1] <- 0 # directly north

# ... simulate movement ----
# We're ready to simulate!

# One useful thing to have pre-calculated are the xy coordinates of each
# raster cell.
coords <- xyFromCell(hab, 1:ncell(hab))

# We'll also want our jitter function from the HSF walkthrough.

# Function to jitter data
jitter <- function(x, y, min = -25, max = 25) {
  res <- data.frame(x = x + runif(1, min, max),
                    y = y + runif(1, min, max))
  return(res)
}

# I have written this code for ease of understanding, 
# not computational speed or handling general use cases. 
# I don't recommend that you try to adapt this code to simulate from an iSSA. 
# Instead, follow the procedure described in Signer et al. (2024) and use
# `amt` to develop the simulation for you.

# Note also that this only works fairly well because we are using
# a small raster and a simple habitat selection model. This code
# will not generalize well to any possible iSSA formulation.

# We already have the first step. We need to simulate the rest.
set.seed(20240429)

for (i in 2:nrow(dat)) {
  # Report status
  cat("\nStep", i, "of", nrow(dat))
  
  ## Calculate selection-free movement kernel
  # Start point
  start <- cbind(dat$x1[i], dat$y1[i])
  # Distances along x and y to every cell
  dx <- coords[, 1] - start[, 1]
  dy <- coords[, 2] - start[, 2]
  # Distance to every cell
  dists <- sqrt(dx^2 + dy^2)
  # Absolute angle to every cell
  abs <- (pi/2 - atan2(dy, dx)) %% (2*pi)
  # Relative angle difference
  rel_diff <- (abs - dat$abs_angle[i-1])
  # Relative angle
  rel_angle <- ifelse(rel_diff > pi, rel_diff - 2*pi, rel_diff)
  # Likelihood of step length
  sl_like <- dgamma(dists, 
                    shape = shp,
                    scale = scl)
  # Likelihood of turn angle 
  suppressWarnings({
    ta_like <- dvonmises(rel_angle,
                         mu = mu,
                         kappa = kappa)
    
  })
  
  # Calculate kernel values
  move_kern_vals <- sl_like * ta_like
  
  # # If you want to plot this
  # move_kern <- hab[[1]]
  # values(move_kern) <- move_kern_vals
  # plot(move_kern, main = "Movement Kernel")
  
  ## Movement-free habitat kernel
  # We pre-computed this
  
  # If you want to plot
  # plot(hab_kern, main = "Habitat Kernel")
  
  ## Combine
  step_kern_vals <- move_kern_vals * hab_kern_vals
  # Normalize
  step_kern_vals <- step_kern_vals/sum(step_kern_vals)
  
  # # If you want to visualize
  # step_kern <- hab[[1]]
  # values(step_kern) <- step_kern_vals
  # plot(step_kern, main = "Habitat x Movement Kernel")
  
  # Randomly select cell to move into based on the probabilities
  next_cell <- sample(x = 1:ncell(hab),
                      size = 1,
                      prob = step_kern_vals)
  
  # Get cell coordinates
  next_cell_coords <- xyFromCell(hab, next_cell)
  
  # If you want to plot
  # points(next_cell_coords[,"x"], next_cell_coords[,"y"], pch = 16)
  
  # Jitter
  next_coords <- jitter(next_cell_coords[, 1], next_cell_coords[, 2])
  
  # Insert into data.frame
  if (i != nrow(dat)) {
    dat$x2[i] <- dat$x1[i+1] <- next_coords[, 1]
    dat$y2[i] <- dat$y1[i+1] <- next_coords[, 2]
  } else {
    dat$x2[i] <- next_coords[, 1]
    dat$y2[i] <- next_coords[, 2]
  }
  
  # Calculate absolute angle
  dx <- dat$x2[i] - dat$x1[i]
  dy <- dat$y2[i] - dat$y1[i]
  dat$abs_angle[i] = (pi/2 - atan2(dy, dx)) %% (2*pi)
  
  # If you want to check
  # dat[i, ]
  
}

# Now that we've simulated our data, let's see what our trajectory looks like.
traj <- dat %>% 
  filter(!is.na(t2)) %>% 
  select(id, x = x2, y = y2, time = t2) %>% 
  make_track(x, y, time, id = id, crs = 32612)

plot(hab[[1]])
points(traj)
lines(traj)

# Save
write.csv(traj, "sim_gps.csv", row.names = FALSE)
