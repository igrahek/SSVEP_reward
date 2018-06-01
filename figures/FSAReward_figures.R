
pkgdir <- paste0(getwd(), "/packages") # set package directory
.libPaths(c(pkgdir, .libPaths())) # add package directory

library(Rmisc) # must be loaded before 'tidyverse' or it will cause compatibility issues
library(tidyverse)
library(viridis)
library(cowplot)
# devtools::install_github('craddm/eegUtils', dependencies = TRUE, lib = pkgdir)
library(eegUtils)

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

################################################################################################
########################################### SPECTRA ############################################
################################################################################################

# adapted from:
# https://craddm.github.io/2016/11/28/erp-visualization-within-subject-confidence-intervals/
spectra <- read_csv(paste0(getwd(), "/figures/spectra.csv")) %>% # load data
  gather(
    key = frequency, # grouping column
    value = amplitude, # value column
    "0":"511.8181818" # columns that need gathering
  ) %>%
  mutate(
    participant = as.factor(participant), # convert participant as factor
    condition = recode(
      factor(condition), # convert condition as factor and rename levels
      "1" = "BslnRedAttended", "2" = "BslnBlueAttended",
      "3" = "AcqRedAttended", "4" = "AcqBlueAttended",
      "5" = "ExtRedAttended", "6" = "ExtBlueAttended"
    ),
    frequency = as.numeric(frequency) # convert frequency as numeric
  ) %>%
  filter(frequency <= 16) # filter out frequencies not in graph

# summarized data from each time point, including within-subject 95% CIs
spectra.pointsummary <- spectra %>%
  split(.$frequency) %>%
  map(~ summarySEwithin(
    data = .,
    measurevar = "amplitude",
    withinvars = "condition",
    idvar = "participant"
  ))

# extract within-subject 95% CIs
spectra..withinssjCI <- map_df(spectra.pointsummary, magrittr::extract) %>%
  mutate(frequency = rep(unique(spectra$frequency), each = length(unique(spectra$condition))))

# plot (frequency range: 0-16 Hz)
ggplot(spectra, aes(frequency, amplitude)) + # basic plot
  stat_summary(aes(colour = condition), # lines: means
    fun.y = mean,
    geom = "line",
    size = 1.3
  ) +
  geom_ribbon( # ribbons: 95% CI
    data = spectra..withinssjCI,
    aes(ymin = amplitude - ci, ymax = amplitude + ci, fill = condition),
    linetype = "solid",
    alpha = .2
  ) +
  scale_color_viridis(discrete = TRUE) +
  scale_fill_viridis(discrete = TRUE) +
  geom_vline(xintercept = seq(0, 16, 2), linetype = "dotted", colour = "#999999", size = .8, alpha = .5) +
  geom_hline(yintercept = seq(0, 1.6, .2), linetype = "dotted", colour = "#999999", size = .8, alpha = .5) +
  scale_y_continuous(breaks = seq(0, 1.6, .2)) +
  scale_x_continuous(breaks = seq(0, 16, 2)) +
  labs(
    title = "Amplitude (0 - 16 Hz)",
    x = "frequency (Hz)",
    y = expression(paste("amplitude (", mu, "V)"))
  ) +
  guides(fill = "none", color = guide_legend(title = NULL)) +
  theme_classic(base_size = 18) +
  theme(
    legend.position = c(.2, .8),
    plot.title = element_text(size = 24, hjust = .5)
  )

# plot (frequency range: 9-11 Hz)
filter(spectra, frequency >= 9 & frequency <= 11) %>%
  ggplot(., aes(frequency, amplitude)) +
  stat_summary(aes(colour = condition),
    fun.y = mean,
    geom = "line",
    size = 1.3
  ) +
  geom_ribbon(
    data = filter(spectra..withinssjCI, frequency >= 9 & frequency <= 11),
    aes(ymin = amplitude - ci, ymax = amplitude + ci, fill = condition),
    linetype = "solid",
    alpha = .2
  ) +
  scale_color_viridis(discrete = TRUE) +
  scale_fill_viridis(discrete = TRUE) +
  geom_vline(xintercept = seq(9, 11, 1), linetype = "dotted", colour = "#999999", size = .8, alpha = .5) +
  geom_hline(yintercept = seq(0, 1.6, .2), linetype = "dotted", colour = "#999999", size = .8, alpha = .5) +
  scale_y_continuous(breaks = seq(0, 1.6, .2)) +
  scale_x_continuous(breaks = seq(9, 11, 1)) +
  coord_cartesian(
    ylim = c(.6, 1.6), # zoom in
    xlim = c(9.5, 10.5)
  ) +
  labs(
    title = "Amplitude (10 Hz)",
    x = "frequency (Hz)",
    y = expression(paste("amplitude (", mu, "V)"))
  ) +
  guides(fill = "none", color = guide_legend(title = NULL)) +
  theme_classic(base_size = 18) +
  theme(
    legend.position = c(.2, .8),
    plot.title = element_text(size = 24, hjust = .5)
  )

# plot (frequency range: 11-13 Hz)
filter(spectra, frequency >= 11 & frequency <= 13) %>%
  ggplot(., aes(frequency, amplitude)) +
  stat_summary(aes(colour = condition),
    fun.y = mean,
    geom = "line",
    size = 1.3
  ) +
  geom_ribbon(
    data = filter(spectra..withinssjCI, frequency >= 11 & frequency <= 13),
    aes(ymin = amplitude - ci, ymax = amplitude + ci, fill = condition),
    linetype = "solid",
    alpha = .2
  ) +
  scale_color_viridis(discrete = TRUE) +
  scale_fill_viridis(discrete = TRUE) +
  geom_vline(xintercept = seq(11, 13, 1), linetype = "dotted", colour = "#999999", size = .8, alpha = .5) +
  geom_hline(yintercept = seq(.5, 1, .1), linetype = "dotted", colour = "#999999", size = .8, alpha = .5) +
  scale_y_continuous(breaks = seq(.5, 1, .1)) +
  scale_x_continuous(breaks = seq(11, 13, 1)) +
  coord_cartesian(
    ylim = c(.5, 1), # zoom in
    xlim = c(11.5, 12.5)
  ) +
  labs(
    title = "Amplitude (12 Hz)",
    x = "frequency (Hz)",
    y = expression(paste("amplitude (", mu, "V)"))
  ) +
  guides(fill = "none", color = guide_legend(title = NULL)) +
  theme_classic(base_size = 18) +
  theme(
    legend.position = c(.2, .8),
    plot.title = element_text(size = 24, hjust = .5)
  )

################################################################################################
################################################################################################
################################################################################################
