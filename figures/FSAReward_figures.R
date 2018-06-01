

pkgdir <- paste0(getwd(), "/packages") # set package directory
.libPaths(c(pkgdir, .libPaths())) # add package directory

library(tidyverse)
library(akima) # for interpolation
library(scales) # scale functions for visualization
library(mgcv) # mixed GAM computation vehicle with automatic smoothness estimation
library(gridExtra) # for graphics

################################################################################################
######################################### TOPOGRAPHIES #########################################
##################### code adapted from https://github.com/craddm/eegUtils ##################### 
################################################################################################

electrodeLocs <- read_delim(paste0(getwd(), "/figures/Antonio_BioSemi64.locs"), # load electrode locations
  "\t", # delimiter (tab)
  col_names = c("chanNo", "theta", "radius", "electrode"), # column names
  escape_double = FALSE, trim_ws = TRUE
)

electrodeLocs$radianTheta <- pi / 180 * electrodeLocs$theta # convert theta values from degrees to radians

# calculate Cartesian coordinates
electrodeLocs <- electrodeLocs %>%
  mutate(
    x = .$radius * sin(.$radianTheta),
    y = .$radius * cos(.$radianTheta)
  )

electrodeLocs <- electrodeLocs[order(electrodeLocs$electrode), ] # sort electrodes in alphabetical order

# create ggplot theme without axes and background grids
theme_topo <- function(base_size=12) {
  theme_bw(base_size = base_size) %+replace%
    theme(
      rect = element_blank(),
      line = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank()
    )
}

# draw head and nose
circleFun <- function(center=c(0, 0), diameter=1, npoints=100) {
  r <- diameter / 2
  tt <- seq(0, 2 * pi, length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}
headShape <- circleFun(c(0, 0), 1, npoints = 100)
nose <- data.frame(x = c(-.075, 0, .075), y = c(.495, .575, .495))

# # plot electrode coordinates on head
# ggplot(headShape, aes(x, y)) +
#   geom_path() +
#   geom_text(data=electrodeLocs, aes(x, y, label=electrode)) +
#   geom_line(data=nose, aes(x, y)) +
#   theme_topo() +
#   coord_equal()

# prepare topography
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")) # jet color map
gridRes <- 268 # number of points for each grid dimension, corresponding to the resolution/smoothness of the interpolation
maskRing <- circleFun(diameter = 1.42) # create a circle around the outside of the plotting area to mask the jagged edges of the interpolation

amplim <- c(0, .301) # amplitude limits
contour.binwidth <- .03 # contour width

topos <- read_csv(paste0(getwd(), "/figures/topos.csv")) %>% # load data
  bind_cols(., electrodeLocs) %>% # bind electrode locations and amplitudes
  select(c("10 Hz", "12 Hz", "x", "y"))

######################################### 10 Hz #########################################

topo_10Hz <- topos %>%
  select(-"12 Hz") %>%
  rename(amp = "10 Hz")

splineSmooth <- gam(amp ~ s(x, y, bs = "ts"), data = topo_10Hz)
GAMtopo <- data.frame(expand.grid(
  x = seq(min(topo_10Hz$x) * 2,
    max(topo_10Hz$x) * 2,
    length = gridRes
  ),
  y = seq(min(topo_10Hz$y) * 2,
    max(topo_10Hz$y) * 2,
    length = gridRes
  )
))

GAMtopo$amplitude <- predict(splineSmooth, GAMtopo, type = "response")
GAMtopo$incircle <- (GAMtopo$x)^2 + (GAMtopo$y)^2 < .7^2

# plot topography
ggplot(GAMtopo[GAMtopo$incircle, ], aes(x, y, fill = amp)) +
  geom_raster() +
  stat_contour(aes(z = amp), binwidth = contour.binwidth) +
  theme_topo() +
  scale_fill_gradientn(
    colours = jet.colors(10),
    limits = amplim,
    guide = "colourbar",
    oob = squish
  ) +
  geom_path(data = maskRing, aes(x, y, z = NULL, fill = NULL), colour = "white", size = 6) +
  geom_point(
    data = topo_10Hz,
    aes(x, y, fill = NULL)
  ) +
  geom_path(data = nose, aes(x, y, z = NULL, fill = NULL), size = 1.5) +
  geom_path(data = headShape, aes(x, y, z = NULL, fill = NULL), size = 1.5) +
  coord_quickmap()




