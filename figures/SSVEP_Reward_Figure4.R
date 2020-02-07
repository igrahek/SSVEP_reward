
######################################### SETUP #########################################

### install packages
# install.packages("here")
# install.packages("Rmisc")
# install.packages("tidyverse")
# install.packages("viridis")
# install.packages("cowplot")
# install.packages("devtools")
# devtools::install_github("r-lib/conflicted")
# devtools::install_github('craddm/eegUtils')
# devtools::install_github("mikabr/ggpirate")

### load packages
library(conflicted)
library(here)
library(Rmisc)
library(tidyverse)
library(viridis)
library(cowplot)
library(eegUtils)
library(ggpirate)

# avoid function conflicts
conflict_prefer("here", "here")
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("mutate", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("summarize", "dplyr")

########################################### FIGURE 4 ###########################################
######################################### TOPOGRAPHIES #########################################

# load electrode locations
electrodeLocs <- read_delim(here("figures/BioSemi64.locs"),
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

topos <- read_csv(here("figures/grandAverage_topos.csv")) %>%
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
    highlights = c("Iz", "Oz", "POz", "O2")
  ) +
  ggtitle("10 Hz") +
  guides(fill = guide_colorbar(
    title = expression(paste("amplitude (", mu, "V)")),
    title.position = "right",
    barwidth = rel(1),
    barheight = rel(6),
    title.theme = element_text(angle = 270)
  )) +
  theme(
    plot.title = element_text(size = 22, hjust = .5), # modify plot title
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
    highlights = c("Iz", "Oz", "POz", "O2")
  ) +
  ggtitle("12 Hz") +
  guides(fill = guide_colorbar(
    title = expression(paste("amplitude (", mu, "V)")),
    title.position = "right",
    barwidth = rel(1),
    barheight = rel(6),
    title.theme = element_text(angle = 270)
  )) +
  theme(
    plot.title = element_text(size = 22, hjust = .5),
    plot.margin = unit(c(6, 0, 6, 0), "pt")
  )

########################################### SPECTRA ############################################

# adapted from:
# https://craddm.github.io/2016/11/28/erp-visualization-within-subject-confidence-intervals/
spectra <- read_csv(here("figures/grandAverage_spectra.csv")) %>% # load data
  select(-c(`16.0455`:`30`)) %>% # keep frequencies until 16 Hz for plotting
  gather(
    key = frequency,
    value = amplitude,
    "0":"16"
  ) %>%
  mutate(
    participant = as.factor(participant),
    frequency = as.numeric(levels(as.factor(frequency)))[as.factor(frequency)],
    Phase = recode(
      factor(condition),
      "1" = "Baseline", "2" = "Baseline", "3" = "Training", "4" = "Training", "5" = "Test", "6" = "Test"
    ),
    attended = recode(
      factor(condition),
      "1" = "red", "2" = "blue", "3" = "red", "4" = "blue", "5" = "red", "6" = "blue"
    ),
    condition = as.factor(paste(Phase, attended, sep = " "))
  ) %>%
  select(participant, frequency, condition, Phase, attended, amplitude) %>%
  filter(participant != "14") %>% # participant eliminated from analyses due to noisy signal
  droplevels()

# extract within-subject 95% CIs
spectra_withinssjCI <- spectra %>%
  split(.$frequency) %>%
  # summarize data from each frequency
  map_df(~ summarySEwithin(
    data = .,
    measurevar = "amplitude",
    withinvars = c("Phase", "attended", "condition"),
    idvar = "participant"
  )) %>%
  mutate(frequency = as.numeric(rep(unique(spectra$frequency), each = length(unique(spectra$condition))))) %>%
  select(frequency, condition, Phase, attended, N, amplitude, ci)

# plot (frequency range: 0-16 Hz)
spectra_all <-
  ggplot(
    spectra_withinssjCI,
    aes(
      x = frequency,
      y = amplitude,
      color = condition,
      linetype = Phase
    )
  ) +
  geom_line(size = 1.3) +
  geom_ribbon( # ribbons: 95% CI
    data = spectra_withinssjCI,
    aes(ymin = amplitude - ci, ymax = amplitude + ci, fill = attended, linetype = Phase),
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
    plot.title = element_text(size = 22, hjust = .5)
  )

########################################### RDI PLOT ###########################################

ssVEP_RDI <-
  spectra %>%
  filter(frequency %in% c(10, 12)) %>%
  mutate(frequency = recode(factor(frequency),
    "10" = "10 Hz", "12" = "12 Hz"
  )) %>%
  ggplot(aes(
    x = Phase,
    y = amplitude,
    color = attended,
    fill = attended
  )) +
  geom_pirate(
    bars = FALSE,
    cis = TRUE,
    lines = TRUE, lines_params = list(color = "black"),
    points = TRUE, points_params = list(shape = 21, size = 5, alpha = .1),
    violins = TRUE, violins_params = list(size = 1),
    show.legend = TRUE
  ) +
  scale_color_manual(values = rep(c("red", "blue"), 3)) +
  scale_fill_manual(values = rep(c("red", "blue"), 3)) +
  scale_y_continuous(
    limits = c(0, 5),
    breaks = seq(0, 5, 1)
  ) +
  ylab("amplitude (a.u.)") +
  facet_wrap(~frequency) +
  theme_minimal(base_size = 16) +
  theme(
    strip.text = element_text(
      hjust = .5,
      size = 20
    ),
    plot.title = element_text(size = 26, hjust = .5),
    legend.box.background = element_rect(color = "transparent"),
    legend.position = "right"
  )

######################################## ARRANGE PLOTS #########################################

# panel A: spectra
cowplot_spectra_all <-
  plot_grid(spectra_all,
    labels = "A"
  )

# panel B: topographies
cowplot_topos <-
  plot_grid(topo_10Hz + theme(legend.position = "none"),
    topo_12Hz + theme(legend.position = "none"),
    align = "vh",
    hjust = -1,
    nrow = 2,
    labels = "C"
  )

# extract legend
cowplot_legend <-
  get_legend(topo_10Hz + theme(legend.box.margin = margin(0, 30, 0, 0)))

# unite plot and legend
cowplot_topos_legend <-
  plot_grid(cowplot_topos,
    cowplot_legend,
    rel_widths = c(1, .1)
  )

# panel C: RDI plot
cowplot_ssVEP_RDI <-
  plot_grid(ssVEP_RDI,
    labels = "B"
  )

# all plots in one figure
figure_4 <-
  plot_grid(cowplot_spectra_all,
    plot_grid(cowplot_ssVEP_RDI,
      cowplot_topos_legend,
      ncol = 2,
      rel_widths = c(1.3, .7)
    ),
    ncol = 1
  )

# save as tiff
ggsave(
  plot = figure_4,
  filename = here("figures/Figure4.tiff"),
  units = "in",
  width = 12, height = 10,
  dpi = 600, compression = "lzw"
)

###################################################################################################
