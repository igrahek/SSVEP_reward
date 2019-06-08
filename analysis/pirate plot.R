EmoSSR.data %>%
  filter(experiment == "pilot1") %>%
  ggplot(aes(
    x = cond4plot, #Exp phase X Reward (all of the combinations)
    y = CS, # amplitude
    color = cond4plot, #attention
    fill = cond4plot #attention
  )) +
  geom_pirate(
    bars = FALSE,
    cis = TRUE,
    lines = TRUE, lines_params = list(color = "black"),
    points = TRUE, points_params = list(shape = 16, color = "black", size = 5, alpha = .2),
    violins = TRUE, violins_params = list(size = 1),
    show.legend = FALSE
  ) +
  scale_y_continuous(
    limits = c(-.1, 1),
    breaks = seq(0, 1, .1)
  ) +
  coord_cartesian(ylim = c(-.06, .4)) +
  geom_hline(
    yintercept = seq(0, 1, .1),
    linetype = "dotted",
    colour = "#999999",
    size = .8,
    alpha = .5
  ) +
  labs(
    x = "",
    y = "CS"
  ) +
  ggtitle("pilot 1") +
  theme_EmoSSR