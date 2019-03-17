
######################################### SETUP #########################################

.libPaths("E:/R_workspace/packages") # add package directory

setwd(paste0(getwd(), "/figures/")) # when using R projects
# setwd("E:/Experiments/Grahek_Ivan/FSAReward/repo/figures")

### install packages
# install.packages(c("tidyverse", "Rmisc",
#                    "parallel", "rstan", "brms",
#                    "viridis", "cowplot"),
#                  dependencies = TRUE)
# devtools::install_github('craddm/eegUtils', dependencies = TRUE)

### load packages
library(Rmisc) # must be loaded before 'tidyverse' or it will cause compatibility issues
library(tidyverse)
library(viridis)
library(cowplot)
library(eegUtils)

######################################### TOPOGRAPHIES #########################################

# load electrode locations
electrodeLocs <- read_delim(paste0(getwd(), "/BioSemi64.locs"),
  "\t",
  col_names = c("chanNo", "theta", "radius", "electrode"),
  escape_double = FALSE, trim_ws = TRUE
) %>%
  mutate(
    radianTheta = pi / 180 * theta, # convert theta values from degrees to radians
    # Cartesian coordinates
    x = radius * sin(radianTheta),
    y = radius * cos(radianTheta)
  )

topos <- read_csv(paste0(getwd(), "/grandAverage_topos.csv")) %>%
  bind_cols(., electrodeLocs) %>% # bind electrode locations and amplitudes
  select(c("10 Hz", "12 Hz", "electrode", "x", "y"))

# 10 Hz
topo_10Hz <-
  topos %>%
  select(-"12 Hz") %>%
  rename(amplitude = "10 Hz") %>%
  topoplot(.,
    limits = c(0, 1), # min-max amplitude
    chanLocs = electrodeLocs,
    method = "Biharmonic",
    palette = "viridis",
    interp_limit = "skirt",
    contour = TRUE,
    chan_marker = "point",
    quantity = "amplitude",
    highlights = c("IZ", "OZ", "POZ", "O2")
  ) +
  ggtitle("10 Hz") +
  guides(fill = guide_colorbar(title = "amplitude (a.u.)",
                               title.position = "right",
                               barwidth = rel(1),
                               barheight = rel(6),
                               title.theme = element_text(angle = 270))) +
  theme(
    plot.title = element_text(size = 22, hjust = .5, face = "bold"), # modify plot title
    plot.margin = unit(c(6, 0, 6, 0), "pt") # decrease plot margins
  )

# 12 Hz
topo_12Hz <-
  topos %>%
  select(-"10 Hz") %>%
  rename(amplitude = "12 Hz") %>%
  topoplot(.,
    limits = c(0, 1),
    chanLocs = electrodeLocs,
    method = "Biharmonic",
    palette = "viridis",
    interp_limit = "skirt",
    contour = TRUE,
    chan_marker = "point",
    quantity = "amplitude",
    highlights = c("IZ", "OZ", "POZ", "O2")
  ) +
  ggtitle("12 Hz") +
  guides(fill = guide_colorbar(title = "amplitude (a.u.)",
                               title.position = "right",
                               barwidth = rel(1),
                               barheight = rel(6),
                               title.theme = element_text(angle = 270))) +
  theme(
    plot.title = element_text(size = 22, hjust = .5, face = "bold"),
    plot.margin = unit(c(6, 0, 6, 0), "pt")
  )

########################################### SPECTRA ############################################

# adapted from:
# https://craddm.github.io/2016/11/28/erp-visualization-within-subject-confidence-intervals/
spectra <- read_csv(paste0(getwd(), "/grandAverage_spectra.csv")) %>% # load data
  dplyr::select(-c(`16.0455`:`30`)) %>% # keep frequencies until 16 Hz for plotting
  gather(
    key = frequency,
    value = amplitude,
    "0":"16"
  ) %>%
  mutate(
    participant = as.factor(participant),
    frequency = as.numeric(levels(as.factor(frequency)))[as.factor(frequency)],
    phase = recode(
      factor(condition),
      "1" = "baseline", "2" = "baseline", "3" = "rewarded", "4" = "rewarded", "5" = "non-rewarded", "6" = "non-rewarded"
    ),
    attended = recode(
      factor(condition),
      "1" = "red", "2" = "blue", "3" = "red", "4" = "blue", "5" = "red", "6" = "blue"
    ),
    condition = as.factor(paste(phase, attended, sep = " "))
  ) %>%
  dplyr::select(participant, frequency, condition, phase, attended, amplitude)

# extract within-subject 95% CIs
spectra.withinssjCI <- spectra %>%
  split(.$frequency) %>%
  # summarize data from each frequency
  map_df(~ summarySEwithin(
    data = .,
    measurevar = "amplitude",
    withinvars = c("phase", "attended", "condition"),
    idvar = "participant"
  )) %>%
  mutate(frequency = as.numeric(rep(unique(spectra$frequency), each = length(unique(spectra$condition))))) %>%
  dplyr::select(frequency, condition, phase, attended, N, amplitude, ci)

# plot (frequency range: 0-16 Hz)
spectra_all <-
  ggplot(
    spectra.withinssjCI,
    aes(
      x = frequency,
      y = amplitude,
      color = condition,
      linetype = phase
    )
  ) +
  geom_line(size = 1.3) +
  geom_ribbon( # ribbons: 95% CI
    data = spectra.withinssjCI,
    aes(ymin = amplitude - ci, ymax = amplitude + ci, fill = attended),
    linetype = "blank",
    alpha = .2
  ) +
  scale_color_manual(values = rep(c("blue", "red"), 3)) +
  scale_fill_manual(values = rep(c("red", "blue"), 3)) +
  guides(fill = "none", color = "none") +
  geom_vline(xintercept = seq(0, 16, 2), linetype = "dotted", color = "#999999", size = .8, alpha = .5) +
  geom_hline(yintercept = seq(0, 1.5, .3), linetype = "dotted", color = "#999999", size = .8, alpha = .5) +
  scale_y_continuous(breaks = seq(0, 1.5, .3)) +
  scale_x_continuous(breaks = seq(0, 16, 2)) +
  coord_cartesian(xlim = c(0, 16), ylim = c(0, 1.5)) +
  labs(
    x = "frequency (Hz)",
    y = "amplitude (a.u.)"
  ) +
  theme_classic(base_size = 18) +
  theme(
    legend.position = c(.2, .7),
    plot.title = element_text(size = 24, hjust = .5)
  )
  









# plot (frequency range: 0-16 Hz)
spectra_all <-
  ggplot(spectra, aes(frequency, amplitude)) +
  stat_summary(aes(color = attended, linetype = phase),
    fun.y = mean,
    geom = "line",
    size = 1.3
  ) +
  geom_ribbon( # ribbons: 95% CI
    data = spectra.withinssjCI,
    aes(ymin = amplitude - ci, ymax = amplitude + ci, fill = attended),
    linetype = "solid",
    alpha = .2
  ) +
  scale_color_manual(values = c("red", "blue")) +
  scale_fill_manual(values = c("red", "blue")) +
  geom_vline(xintercept = seq(0, 16, 2), linetype = "dotted", color = "#999999", size = .8, alpha = .5) +
  geom_hline(yintercept = seq(0, 1.5, .3), linetype = "dotted", color = "#999999", size = .8, alpha = .5) +
  scale_y_continuous(breaks = seq(0, 1.5, .3)) +
  scale_x_continuous(breaks = seq(0, 16, 2)) +
  coord_cartesian(xlim = c(0, 16), ylim = c(0, 1.5)) +
  labs(
    x = "frequency (Hz)",
    y = "amplitude (a.u.)"
  ) +
  theme_classic(base_size = 18) +
  theme(
    legend.position = c(.2, .8),
    plot.title = element_text(size = 24, hjust = .5)
  )

# zoom in to relevant frequencies
# plot (frequency range: 9-11 Hz)
spectra_10Hz <- 
  filter(spectra, frequency >= 9 & frequency <= 11) %>%
  ggplot(., aes(frequency, amplitude)) +
  stat_summary(aes(color = attended, linetype = phase),
    fun.y = mean,
    geom = "line",
    size = 1.3
  ) +
  geom_ribbon(
    data = filter(spectra.withinssjCI, frequency >= 9 & frequency <= 11),
    aes(ymin = amplitude - ci, ymax = amplitude + ci, fill = attended),
    linetype = "solid",
    alpha = .2
  ) +
  scale_color_manual(values = c("red", "blue")) +
  scale_fill_manual(values = c("red", "blue"))
  
  # geom_vline(xintercept = seq(9, 11, 1), linetype = "dotted", color = "#999999", size = .8, alpha = .5) +
  geom_hline(yintercept = seq(0, 1.5, .1), linetype = "dotted", color = "#999999", size = .8, alpha = .5) +
  scale_y_continuous(breaks = seq(0, 1.5, .1)) +
  scale_x_continuous(breaks = seq(9, 11, 1)) +
  coord_cartesian( # zoom in
    ylim = c(.5, 1.5),
    xlim = c(9.5, 10.5)
  ) +
  labs(
    title = "10 Hz",
    x = "",
    y = expression(paste("amplitude (", mu, "V)"))
  ) +
  guides(fill = "none", color = guide_legend(title = NULL)) +
  theme_classic(base_size = 18) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 24, hjust = .5)
  )

# plot (frequency range: 11-13 Hz)
spectra_12Hz <- 
  filter(spectra, frequency >= 11 & frequency <= 13) %>%
  ggplot(., aes(frequency, amplitude)) +
  stat_summary(aes(color = condition),
    fun.y = mean,
    geom = "line",
    size = 1.3
  ) +
  geom_ribbon(
    data = filter(spectra.withinssjCI, frequency >= 11 & frequency <= 13),
    aes(ymin = amplitude - ci, ymax = amplitude + ci, fill = condition),
    linetype = "solid",
    alpha = .2
  ) +
  scale_color_viridis(discrete = TRUE) +
  scale_fill_viridis(discrete = TRUE) +
  geom_vline(xintercept = seq(11, 13, 1), linetype = "dotted", color = "#999999", size = .8, alpha = .5) +
  geom_hline(yintercept = seq(.5, 1, .1), linetype = "dotted", color = "#999999", size = .8, alpha = .5) +
  scale_y_continuous(breaks = seq(.5, 1, .1)) +
  scale_x_continuous(breaks = seq(11, 13, 1)) +
  coord_cartesian( # zoom in
    ylim = c(.5, 1.5),
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

######################################## ARRANGE PLOTS #########################################

# arrange topographies in a single row
topo.row <- plot_grid(topo_10Hz + theme(legend.position = "none"),
  topo_12Hz + theme(legend.position = "none"),
  align = "vh",
  hjust = -1,
  nrow = 1,
  scale = .8 # slightly smaller size
)

legend <- get_legend(topo_10Hz) # extract legend

plot_grid(topo.row,
  legend,
  spectra_all,
  rel_widths = c(3, .3)
) %>% # the second column (containing the legend) is 10 times smaller than the first one (containing the two topographies)
  save_plot(paste0(getwd(), "/topos_spectra.jpg"),
    .,
    base_height = 10,
    base_aspect_ratio = 1.1
  )
