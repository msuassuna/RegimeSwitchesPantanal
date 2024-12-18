plot(df_wide$hydroYear, df_wide$`66825000`, type = "l", bty = "n", col = "grey",
     xlab = "Year", ylab = "Water level (cm)",
     main = "Annual minimum water levels at Ladario stream gauge")



plotData <- data.frame("Year" = df_wide$hydroYear,
                       "Amf" = df_wide$`66825000`,
                       "Class" = prstates_total)
plotData %>%
  ggplot(aes(x = factor(Class),
             y = Amf,
             fill = factor(Class))) +
  ggdist::stat_halfeye(
    adjust = 0.5
  )

plotData %>%
  ggplot(aes(x = factor(Class),
             y = Amf,
             fill = factor(Class))) +
  ggdist::stat_halfeye(
    adjust = 0.5,
    justification = -.1,
    .width = 0,
    point_colour = NA) +
  geom_boxplot(
    width = 0.12,
    alpha = 0.5) +
  ggdist::stat_dots(
    side = "left",
    justification = 1.1,
    binwidth = 5) +
  scale_fill_manual(values=c("#FF8033CC", "#33CC33CC", "#3333CCCC")) +
  theme_tq() +
  labs(
    title = "Annual minimum flows at the Ladário streamgauge",
    subtitle = "Classification on different drought classes was performed using a multisite HMM",
    x = "Drought Classes",
    y = "Annual mimimum flows (m³/s)",
    fill = "Drought classes"
  )
  