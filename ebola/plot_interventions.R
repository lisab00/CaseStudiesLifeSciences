library(tidyverse)
library(ggplot2)
library(geomtextpath)
library(latex2exp)

par_shift <- function(t,par_old,par_new, trans_speed){
  if (t < 5){
    return(par_old)
  } else {
    return(par_new + (par_old - par_new) * exp(- trans_speed * (t - 5)))
  }
}

par_shift_plot <- function(t){
    return(par_shift(t, 1,0,1))
}
    t <- seq(1,12, by = 0.001)
par_shift_vals <- data.frame(time = t, value =  sapply(t, par_shift_plot))

par_shift_vals %>%
  ggplot(aes(x = time, y = value)) +
  geom_line(linewidth = 1.5) +
  geom_textvline(xintercept = 5, color = "orange", linetype = "dashed", linewidth = 1, label = "time of intervention",
                 hjust = 0.99, vjust = -0.5) +
  labs(x="time",
       y="parameter value") +
  theme_gray(base_size=15) +
  scale_x_continuous(breaks=1:12) +
  ylim(0.0, 1.0) -> plot1

ggsave(plot = plot1, "ebola_par_shift_decline.pdf", width = 8, height = 5)
