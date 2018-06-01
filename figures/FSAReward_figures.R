
pkgdir <- paste0(getwd(), "/packages") # set package directory
.libPaths(c(pkgdir, .libPaths())) # add package directory

library(tidyverse)
library(eegUtils)
library(viridis)
library(cowplot)

################################################################################################
######################################### TOPOGRAPHIES #########################################
################################################################################################

electrodeLocs <- read_delim(paste0(getwd(), "/figures/Antonio_BioSemi64.locs"), # load electrode locations
  "\t", # delimiter (tab)
  col_names = c("chanNo", "theta", "radius", "electrode"), # column names
  escape_double = FALSE, trim_ws = TRUE
) %>%
  mutate(
    radianTheta = pi / 180 * theta, # convert theta values from degrees to radians
    # Cartesian coordinates
    x = radius * sin(radianTheta),
    y = radius * cos(radianTheta)
  )

topos <- read_csv(paste0(getwd(), "/figures/topos.csv")) %>% # load data
  bind_cols(., electrodeLocs) %>% # bind electrode locations and amplitudes
  select(c("10 Hz", "12 Hz", "electrode", "x", "y"))

# 10 Hz
topo_10Hz <- topos %>%
  select(-"12 Hz") %>%
  rename(amplitude = "10 Hz") %>%
  topoplot(.,
    time_lim = NULL,
    limits = c(0, 1),
    chanLocs = electrodeLocs,
    method = "Biharmonic",
    r = NULL,
    grid_res = 67,
    palette = "viridis",
    interp_limit = "skirt",
    contour = TRUE,
    chan_marker = "point",
    quantity = "amplitude",
    montage = NULL
  ) + theme(plot.margin = unit(c(6, 0, 6, 0), "pt"))

# 12 Hz

topo_12Hz <- topos %>%
  select(-"10 Hz") %>%
  rename(amplitude = "12 Hz") %>%
  topoplot(.,
    time_lim = NULL,
    limits = c(0, 1),
    chanLocs = electrodeLocs,
    method = "Biharmonic",
    r = NULL,
    grid_res = 67,
    palette = "viridis",
    interp_limit = "skirt",
    contour = TRUE,
    chan_marker = "point",
    quantity = "amplitude",
    montage = NULL
  ) + theme(plot.margin = unit(c(6, 0, 6, 0), "pt"))

# arrange the plots in a single row
prow <- plot_grid(topo_10Hz + theme(legend.position = "none"),
  topo_12Hz + theme(legend.position = "none"),
  align = "vh",
  labels = c("A", "B"),
  hjust = -1,
  nrow = 1
)

legend <- get_legend(topo_10Hz) # extract legend
plot_grid(prow, legend, rel_widths = c(3, .3)) # plots & legend






