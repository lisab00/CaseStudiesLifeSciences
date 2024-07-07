# =============================================================================
# TO DO: Plots for Report
# Author: Andrea Pinke
# =============================================================================
library(tidyverse)
library(ggplot2)
library(geomtextpath)
library(latex2exp)

# Fluctuating functions adapted from Plague paper

# f.seas --> seasonality
f.seas <- function(t){
  A <- 1.15
  B <- 0.08
  C <- 0.1
  return(A + B*sin(2*pi/12*t) + C*cos(2*pi/12*t))
}

# Plot seasonality by month
# highlight outbreak period since for us only that interval is relevant 
# --> months 8,9,10,11 (122 days in total, data for first 118 days)
t <- seq(1,12, by = 0.001)
f.seas_vals <- data.frame(time = t, value = f.seas(t))

f.seas_vals %>%
  ggplot(aes(x = time, y = value)) +
  geom_line(linewidth = 1.5) +
  geom_textvline(xintercept = 11.87, color = "orange", linetype = "dashed", linewidth = 1, label = "26/11/2017",
                 hjust = 0.99, vjust = -0.5) +
  geom_textvline(xintercept = 8, color = "orange", linetype = "dashed", linewidth = 1, label = "01/08/2017",
                 hjust = 0.99, vjust = -0.5) +
  labs(x="Month",
       y="Normalized Temperature (\u00B0C)") +
  theme_gray(base_size=15) +
  scale_x_continuous(breaks=1:12) +
  ylim(0.9, 1.5) -> plot1

ggsave(plot = plot1, "seasonality.pdf", width = 8, height = 5)


# Plot intervention effect starting at the relevant outbreak period (intervention when sudden increase in numbers)
f.int <- function(t,theta){
  return(1-1/(1+exp(theta-t)))
}

delay <- 26

t <- seq(1, 66, by = 0.001)
f.int_vals <- data.frame(time = t, value = f.int(t, delay))

f.int_vals %>%
  ggplot(aes(x = time)) +
  geom_line(aes(y = value), linewidth = 1.5) +
  geom_textvline(xintercept = delay, color = "orange", linetype = "dashed", linewidth = 1, label = "half-maximum effect",
                 hjust = 0.99, vjust = 1.5) +
  labs(x="Days since Intervention",
       y="Change in Transmission") +
  annotate("text", x = delay, y = -0.09, label = TeX("$\\tau$"), color = "orange", size = 6) +  
  coord_cartesian(ylim = c(0,1), clip = "off") +
  theme_gray(base_size=15) -> plot2



ggsave(plot = plot2, "interventions_vector.pdf", width = 8, height = 5)


# Plot adapted delay for detection
f.det <- function(t,tau){
  return(1/(1+exp(tau-t)))
}

delay <- 26

t <- seq(1, 66, by = 0.001)
f.det_vals <- data.frame(time = t, value = f.det(t, delay))

f.det_vals %>%
  ggplot(aes(x = time)) +
  geom_line(aes(y = value), linewidth = 1.5) +
  geom_textvline(xintercept = delay, color = "orange", linetype = "dashed", linewidth = 1, label = "half-maximum effect",
                 hjust = 0.99, vjust = -0.5) +
  labs(x="Days since Start of Large Epidemic Wave",
       y="Effect of Detection") +
  annotate("text", x = delay, y = -0.09, label = TeX("$\\tau$"), color = "orange", size = 6) +  
  coord_cartesian(ylim = c(0,1), clip = "off") +
  theme_gray(base_size=15) -> plot2

ggsave(plot = plot2, "delay_detection.pdf", width = 8, height = 5)



# Data
df_plague_relevant %>%
  ggplot(aes(x = time, y = I1)) +
  geom_point(color="lightseagreen", size=2) +
  labs(x="days since outbreak",
       y="infected individuals") +
  theme_gray(base_size=15) 


# Fitted curve

out_df <- data.frame(out_fit) #out_fit result from fit R file

df_plague_relevant %>%
  ggplot(aes(x = time, y = I1)) +
  geom_point(color="lightseagreen", size=2) +
  labs(x="Days since Start of Large Epidemic Wave",
       y="Number of Infected Individuals") +
  #geom_vline(xintercept = delay, color = "orange", linetype = "dashed", linewidth = 1) +
  #annotate("text", x = delay, y = -5, label = TeX("$\\tau$"), color = "orange", size = 6) +  
  #coord_cartesian(ylim = c(0,70), clip = "off") +
  theme_gray(base_size=15) +
  geom_line(aes(x=time, y=out_df$I1), size=1) -> plot3

ggsave(plot = plot3, "fitted_curve.pdf", width = 6, height = 5)