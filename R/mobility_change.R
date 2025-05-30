source("R/lib.R")

source("R/functions.R")

# load, cache, and format the mobility data
mobility <- all_mobility() %>%
  append_google_data()

write_mobility_dates(mobility)

saveRDS(mobility, file = "outputs/cached_mobility.RDS")
# mobility <- readRDS("outputs/cached_mobility.RDS")

n_weeks_ahead <- 6
first_date <- min(mobility$date)
last_date <- max(mobility$date)

mobility_fitted <- mobility %>%
  rename(
    state_long = state,
  ) %>%
  filter(
    !is.na(state_long) &
      !is.na(trend)
  ) %>%
  mutate(
    state = abbreviate_states(state_long)
  ) %>%
  group_by(state, datastream) %>%
  do(
    predict_mobility_trend(
      .,
      min_date = first_date,
      max_date = last_date + 7 * n_weeks_ahead
    )
  ) %>%
  ungroup()

all_states <- na.omit(unique(mobility_fitted$state_long))

mobility_ticks_labels <- split_ticks_and_labels(
  data = mobility_fitted,
  tick_freq = "2 month",
  label_freq = "4 months",
  label_format = "%b %y"
)

for (this_state in all_states) {
  
  mobility_fitted %>%
    filter(state_long == this_state) %>%
    ggplot() +
    aes(date, fitted_trend) +
    geom_hline(
      yintercept = 0,
      colour = "grey80",
      linetype = 3
    ) +
    geom_vline(
      aes(xintercept = date),
      data = interventions() %>%
        filter(state == abbreviate_states(this_state)),
      colour = "grey80"
    ) +
    geom_vline(
      aes(xintercept = last_date),
      colour = "grey80",
      linetype = 2
    ) +
    facet_wrap(
      ~datastream,
      ncol = 3,
      scales = "free"
    ) +
    geom_ribbon(
      aes(
        ymin = fitted_trend_lower,
        ymax = fitted_trend_upper
      ),
      fill = grey(0.9),
      colour = grey(0.8),
      size = 0.1
    ) +
    # fitted trend
    geom_line(
      aes(date, fitted_trend),
      colour = "gray40"
    ) +
    geom_point(
      aes(date, trend),
      size = 0.2,
      col = "purple"
    ) +
    # # predicted trend
    # geom_line(
    #   aes(date, predicted_trend),
    #   size = 1
    # ) +
    coord_cartesian(
      xlim = c(as.Date("2020-03-01"), last_date)# + 7 * n_weeks_ahead)
    ) +
    scale_y_continuous(position = "right") +
    scale_x_date(
      breaks = mobility_ticks_labels$ticks,
      labels = mobility_ticks_labels$labels,
      limits = range(mobility_fitted$date)
    ) +
    xlab("") +
    ylab("") +
    ggtitle(
      sprintf(
        "Percentage change in selected mobility datastreams up to %s, %s",
        format(last_date, format = "%B %d"),
        format(last_date, format = "%Y")
      )
    ) +
    cowplot::theme_cowplot() +
    cowplot::panel_border(remove = TRUE) +
    theme(legend.position = "none",
          strip.background = element_blank(),
          strip.text = element_text(hjust = 0, face = "bold"),
          axis.text.x = element_text(size = 10),
          axis.title.y.right = element_text(vjust = 0.5, angle = 90),
          panel.spacing = unit(1.2, "lines"))
  
  dpi <- 150
  ggsave(
    filename = sprintf(
      "outputs/figures/%s_datastream_model_fit_%s.png",
      this_state,
      last_date
    ),
    width = 1500 / dpi,
    height = 1250 / dpi,
    dpi = dpi,
    scale = 1.2
  )
}

# save predictions in correct format for macro and mobility models
mobility_fitted %>%
  filter(
    grepl("^Google: ", datastream)
  ) %>%
  mutate(
    change = 1 + (predicted_trend / 100),
    state_datastream = paste(state_long, datastream)
  ) %>%
  select(
    state_datastream,
    state = state_long,
    datastream,
    change,
    date
  ) %>%
  saveRDS("outputs/google_change_trends.RDS")

# output 3-column plot
target_datastreams <- c("Google: time at workplaces",
                        "Google: time at retail and recreation",
                        "Google: time at transit stations")

mobility_fitted %>%
  filter(
    datastream %in% target_datastreams
  ) %>%
  mutate(
    datastream = factor(datastream, levels = target_datastreams)
  ) %>%
  ggplot() +
  aes(date, fitted_trend) +
  geom_hline(
    yintercept = 0,
    colour = "grey80",
    linetype = 3
  ) +
  geom_vline(
    aes(xintercept = last_date),
    colour = "grey80",
    linetype = 2
  ) +
  facet_grid(
    rows = vars(state),
    cols = vars(datastream),
    switch = "y",
    scales = "free"
  ) +
  geom_vline(
    aes(xintercept = date),
    data = interventions(),
    colour = "grey80"
  ) +
  geom_ribbon(
    aes(
      ymin = fitted_trend_lower,
      ymax = fitted_trend_upper
    ),
    fill = grey(0.9),
    colour = grey(0.8),
    size = 0.1
  ) +
  # fitted trend
  geom_line(
    aes(date, fitted_trend),
    colour = "gray40",
    size = 0.25
  ) +
  geom_point(
    aes(date, trend),
    size = 0.2,
    col = "purple"
  ) +
  coord_cartesian(
    xlim = c(as.Date("2020-03-01"), last_date)
  ) +
  scale_y_continuous(
    position = "right",
    breaks = seq(-100, 500, by = 25)
  ) +
  scale_x_date(
    date_breaks = "2 months",
    date_labels = "%m",
    limits = range(mobility_fitted$date)
  ) +
  xlab("") +
  ylab("") +
  ggtitle(
    sprintf(
      "Percentage change in selected mobility datastreams up to %s, %s",
      format(last_date, format = "%B %d"),
      format(last_date, format = "%Y")
    )
  ) +
  cowplot::theme_cowplot() +
  cowplot::panel_border(remove = TRUE) +
  theme(legend.position = "none",
        strip.background = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        strip.text.y = element_text(size = 18),
        strip.text.x = element_text(size = 12),
        plot.title = element_text(size = 20, face = "plain"),
        panel.border = element_rect(colour = grey(0.4), fill = NA),
        panel.spacing = unit(1.2, "lines")) 

dpi <- 300
ggsave(
  filename = sprintf(
    "outputs/figures/multistate_model_fit_%s.png",
    last_date
  ),
  width = 2481 / dpi,
  height = 3507 / dpi,
  dpi = dpi,
  scale = 1.2
)


mobility_date_six_month <- last_date - months(6)

mobility_fitted %>%
  filter(
    datastream %in% target_datastreams,
    date >= mobility_date_six_month
  ) %>%
  mutate(
    datastream = factor(datastream, levels = target_datastreams)
  ) %>%
  ggplot() +
  aes(date, fitted_trend) +
  geom_hline(
    yintercept = 0,
    colour = "grey80",
    linetype = 3
  ) +
  geom_vline(
    aes(xintercept = last_date),
    colour = "grey80",
    linetype = 2
  ) +
  facet_grid(
    rows = vars(state),
    cols = vars(datastream),
    switch = "y",
    scales = "free"
  ) +
  geom_vline(
    aes(xintercept = date),
    data = interventions(),
    colour = "grey80"
  ) +
  geom_ribbon(
    aes(
      ymin = fitted_trend_lower,
      ymax = fitted_trend_upper
    ),
    fill = grey(0.9),
    colour = grey(0.8),
    size = 0.1
  ) +
  # fitted trend
  geom_line(
    aes(date, fitted_trend),
    colour = "gray40",
    size = 0.25
  ) +
  geom_point(
    aes(date, trend),
    size = 0.2,
    col = "purple"
  ) +
  coord_cartesian(
    xlim = c(mobility_date_six_month, last_date)
  ) +
  scale_y_continuous(
    position = "right",
    breaks = seq(-100, 500, by = 25)
  ) +
  scale_x_date(
    date_breaks = "1 months",
    date_labels = "%b %y",
    limits = range(mobility_fitted$date)
  ) +
  xlab("") +
  ylab("") +
  ggtitle(
    sprintf(
      "Percentage change in selected mobility datastreams up to %s, %s",
      format(last_date, format = "%B %d"),
      format(last_date, format = "%Y")
    )
  ) +
  cowplot::theme_cowplot() +
  cowplot::panel_border(remove = TRUE) +
  theme(legend.position = "none",
        strip.background = element_blank(),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(angle = 270),
        #strip.text.x = element_text(size = 10),
        strip.text.y = element_text(size = 18),
        strip.text.x = element_text(size = 12),
        plot.title = element_text(size = 20, face = "plain"),
        panel.border = element_rect(colour = grey(0.4), fill = NA),
        panel.spacing = unit(1.2, "lines")) 

dpi <- 300
ggsave(
  filename = sprintf(
    "outputs/figures/multistate_model_fit_six_month_%s.png",
    last_date
  ),
  width = 2481 / dpi,
  height = 3507 / dpi,
  dpi = dpi,
  scale = 1.2
)
