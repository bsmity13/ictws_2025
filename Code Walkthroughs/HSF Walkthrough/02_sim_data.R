#######################################################X
#----------------------ICTWS 2025----------------------X
#--------------Habitat Selection Workshop--------------X
#----------------Last updated 2025-03-17---------------X
#-----------------HSF Code Walkthrough-----------------X
#######################################################X

# Simulate point locations for HSF walkthrough

# Load packages ----
library(tidyverse)
library(terra)
library(amt)

# Generate example data ----
# We need to generate some example telemetry data to fit our HSF.
# Generating points is fairly easy since the HSF assumes all points are 
# independent. However, we should spend some time thinking about how
# we define the relationships between our study animal and habitat.

# ... habitat axes ----
# First we need habitat variables to define our e-space.

# It is often useful to break our variables into:
#   - resources (more is better)
#   - conditions (some intermediate value is best)
#   - risk (less is better)

# I've already randomly generated 4 variables:
#   - forage (g/m^2) [resource]
#   - mean annual temperature (°C) [condition]
#   - predator density (predators/100 km^2) [risk]
#   - land cover type (grassland = 1, forest = 2, wetland = 3) [condition]

# If you want to see how I did that, check out the script "habitat.R". We can 
# load those here as a SpatRast.

hab <- rast("geo/habitat.tif")
names(hab) <- c("forage", "temp", "predator", "cover")

# Cover is a factor, let's code it that way
levels(hab[[4]]) <- data.frame(id = 1:3,
                               cover = c("grassland", "forest", "wetland"))
hab$cover
plot(hab)

# ... habitat relationships ----
# Now we have to define our relationship between density and our habitat axes.
# In general:
#   - resources should have a single, positive slope parameter (more is better)
#   - risks should have a single, negative slope parameter (less is better)
#   - conditions should be defined by a concave-down parabola (intermediate is best)
#       - a parabola requires a linear and quadratic term
#       - the x-coordinate of the vertex will be at -b/2a
#           - a is the beta for the quadratic term
#           - b is the beta for the linear term

# Note: we are ignoring the intercept here, so the y-axis is just 
# *relative* density.

# We will also assume (for now) that our habitats are strictly additive
# (no interactions).

# ... ... resources ----
# Let's begin with resources.
# We know we want our beta to be positive. How big should it be?
hist(hab$forage)

# Our forage variable ranges from 0 to 1000 (g/m^2). Let's say that at 1000 
# g/m^2, the density is 5x the density at 500 g/m^2. I.e., RSS(x1, x2) = 5 if
# x1 has forage == 1000 and x2 has forage == 500. 

# The beta is just the log of RSS for a 1-unit change, and it is also the
# difference between the linear predictors. I.e.,
# log(5) = (beta_forage * 1000) - (beta_forage * 500)
#   ==>
beta_forage = log(5)/500

# Let's check our work. The RSS for 1000 vs 500 g/m^2 should be 5.
exp(beta_forage * 1000) / exp(beta_forage * 500)

# How many times more animals do we have at 500 g/m^2 than 200 g/m^2?
exp(beta_forage * 500) / exp(beta_forage * 200)

# How many times more animals do we have at 1000 g/m^2 than 0 g/m^2?
exp(beta_forage * 1000) / exp(beta_forage * 0)

# ... ... risks ----
# Now onto risks.
# We know we want our beta to be negative, but again, how big?
hist(hab$predator)

# Our predator density ranges from 0 - 12 (predators / 100 km^2).
# Let's say that we have 0.25x as many animals at predator == 10 than
# predator == 5. I.e., RSS(x1, x2) = 0.25.

# As above,
# log(0.25) = (beta_pred * 10) - (beta_pred * 5)
#   ==>
beta_pred = log(0.25)/5

# Let's check our work.
exp(beta_pred * 10) / exp(beta_pred * 5)

# How many MORE animals will we have if predator density is 0 than if predator
# density is 12?
exp(beta_pred * 0) / exp(beta_pred * 12)

# Note that the overall change in density of our critter changes ~ 25x as we
# move through the full range of the risk and resource variables. But the
# magnitude of 'beta_forage' is ~ 100x smaller than the magnitude of 
# 'beta_pred'. That's because the range of forage is about 100x greater than
# the range of predator densities. I.e., be careful not to interpret the
# magnitude of the beta if you haven't standardized your covariates.

# ... ... conditions ----
# Now onto conditions.
# This one is tricky because we are trying to parameterize a parabola, so
# we need 2 parameters. Let's start by deciding where the vertex (the best
# habitat) should be.
#
# For a parabola of the form y = ax^2 + bx + c,
# the x-coordinate of the vertex falls at -b/2a.
#
# What does our distribution of temperatures look like?
hist(hab$temp)

# Say our critter prefers a mean annual temperature of 13 degrees C, i.e.,
# -b/2a = 13, or to give our parameters better names,
# (-1 * beta_temp) / (2 * beta_temp2) = 13
#   ==>
# beta_temp/beta_temp2 = 13 * -2 = -26

# Now we need a second piece of information to make the parameters identifiable.
# So let's say our our animal density is double at 13 degrees vs 7 degrees.
# I.e., RSS(x1, x2) = 2
#   ==>
# (beta_temp2 * 13^2 + beta_temp * 13) - (beta_temp2 * 7^2 + beta_temp * 7) = log(2)

# I won't try to type out the algebra here, but if you solve this system
# of equations, you'll find that:
beta_temp2 = -1 * log(2)/36
beta_temp = beta_temp2 * -26

# Let's check our work.
exp(beta_temp2 * 13^2 + beta_temp * 13) / exp(beta_temp2 * 7^2 + beta_temp * 7)

# How many times more animals will we have at temp == 10 than temp == 5?
exp(beta_temp2 * 10^2 + beta_temp * 10) / exp(beta_temp2 * 5^2 + beta_temp * 5)

# How many times more animals will we have at temp == 13 than temp == 2.5?
exp(beta_temp2 * 13^2 + beta_temp * 13) / exp(beta_temp2 * 2.5^2 + beta_temp * 2.5)

# "Most preferred temperature"
-1*beta_temp/(2*beta_temp2)

# ... ... land cover ----
# Finally, we tackle land cover.
# One of the categories will be the reference category. I.e., it would be
# captured by the "intercept", but since we don't have one, it takes on the
# values from all our other variables. We'll make grassland our reference.

# The beta for forest is then the log-RSS for forest vs. grassland.
# The beta for wetland is then the log-RSS for wetland vs grassland.

# Let's say density is twice as high in forest as in grassland. I.e.,
beta_forest <- log(2)

# Let's say density is half as high in wetland as in grassland. I.e.,
beta_wetland <- log(1/2)

# Dummy variables
#     Forest  Grassland Wetland
# F   1       0         0
# W   0       0         1

# Okay, we have all of our betas!

# ... calculate w(x) ----
# Recall that w(x) will be proportional to our number of points. We will decide
# later how many total data points to collect, but for now, let's calculate
# w(x) for all locations in space.

# Get our raster data into a data.frame
dat <- as.data.frame(hab, xy = TRUE) %>%
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
  # Calculate w(x)
  mutate(w = exp(g))

# ... calculate expected number of points ----
# Let's say that all of our simulated data will come from a single individual
# wearing a GPS collar.
#   - Assume our GPS fixes are far enough apart that there is no autocorrelation
#   - Assume the entire raster extent is equally available.

# How many points should we have? Let's say we've deployed our collar long
# enough to have ~ 1000 locations.

# What's the expected number of points in each cell?

# Recall, under the IPP, lambda(s) is proportional to w(x(s)). Our expected
# number of points in each raster cell is 1000 * w'(x), where w'(x) is a 
# normalized version of w(x), i.e., it sums to 1.
dat$w_prime <- dat$w/sum(dat$w)

# Expected points in each cell
dat$lambda <- 1000 * dat$w_prime
head(dat)
hist(dat$lambda)

# ... draw points ----
# Finally, we can draw n, the number of locations in each cell. Recall, this
# is a Poisson IPP, so we will draw from the Poisson distribution.
set.seed(20240429)
dat$n <- rpois(n = nrow(dat), lambda = dat$lambda)

# This is all we need if we assume we would use the Poisson GLM to fit our
# HSF.
dat %>% 
  dplyr::select(x, y, forage, temp, predator, cover, n) %>% 
  head()

# However, if we really want GPS coordinates, we can take the coordinates of
# each cell with n > 1, and randomly place n points within that cell. Let's do 
# that.

# Our coordinates are the center of our raster cell, and each cell is 50m x 50m.
# So we can jitter up to 25m in any direction and still be in the same cell.

# Function to jitter data
jitter <- function(x, y, min = -25, max = 25) {
  res <- data.frame(x = x + runif(1, min, max),
                    y = y + runif(1, min, max))
  return(res)
}

# E.g.,
jitter(0, 0)

# Now we split each row with n > 1 into a list element
dat_list <- dat %>% 
  filter(n > 0) %>%
  split(1:nrow(.))

# And now we can create jittered points for each element of our list.
# We will use 'bind_rows()' (twice) to return a single data.frame
set.seed(20220126)
gps <- lapply(dat_list, function(d) {
  replicate(d$n, jitter(x = d$x, y = d$y), simplify = FALSE) %>% 
    bind_rows()
}) %>% 
  bind_rows()

# This is what the result looks like.
head(gps)

plot(hab$forage)
points(gps$x, gps$y, pch = 16, cex = 0.5)

# Save the results
write.csv(gps, "sim_gps.csv", row.names = FALSE)
