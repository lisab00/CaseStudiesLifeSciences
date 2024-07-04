'
This script contains all plots needed for the norovirus modelling.
'


## load required data and functions
library(ggplot2)
library(readxl)
source('norovirus_model_definition.R')
df_noro <- as.data.frame(read_excel("data/relevantData.xlsx", sheet="outbreak11"))
ti <- 6


## Plot the relevant data set
ggplot(df_noro, aes(x = time, y = I1)) +
  geom_vline(xintercept = ti, color = 'black', linetype='dashed', size=1) + 
  geom_point(color = 'lightseagreen', size = 2) +  
  annotate("text", x = ti, y = max(df_noro$I1), label = "Intervention", 
           color = "black", hjust = -0.1, vjust = 1.5) +  
  labs(title = 'Norovirus outbreak data',
       x = 'Days since outbreak', y = 'Number of newly infected') +
  scale_x_continuous(breaks = seq(min(df_noro$time), max(df_noro$time), by = 1)) +
  theme_gray()


## plot intervention functions
t=seq(df_noro$time[1], df_noro$time[length(df_noro$time)], 0.001)
df_signals <- data.frame(
  t=t,
  signalP = sapply(t, signalP, ti = ti, effP = 0.7),
  signalD = sapply(t, signalD, ti = ti, q = 1)
)

# plot for signalD
ggplot(df_signals, aes(x = t)) +
  geom_line(aes(y = signalD), size = 1, col="lightseagreen") +
  geom_vline(xintercept = ti, color = 'black', linetype='dashed', size=1) +
  annotate("text", x = ti, y = 1, label = "Health intervention", 
           color = "black", hjust = -0.1, vjust = 1.5) +  
  labs(title = 'Intervention function for D1 compartment', x = 'Days since outbreak') +
  scale_x_continuous(breaks = seq(min(df_noro$time), max(df_noro$time), by = 1)) +
  theme_gray() +
  ylim(0, 1)

#plot for signalP
ggplot(df_signals, aes(x = t)) +
  geom_line(aes(y = signalP), size = 1, col="lightseagreen") +
  geom_vline(xintercept = ti, color = 'black', linetype='dashed', size=1) +
  annotate("text", x = ti, y = 1, label = "Health intervention", 
           color = "black", hjust = -0.1, vjust = 1.5) +  
  labs(title = 'Intervention function for P compartment', x = 'Days since outbreak') +
  scale_x_continuous(breaks = seq(min(df_noro$time), max(df_noro$time), by = 1)) +
  theme_gray() +
  ylim(0, 1)


## plot model simulation against data
plot_against_data <- function(df_obs, df_out_fit){
  ggplot() +
    geom_line(data = df_out_fit, aes(x = time, y = I1), color = 'black', size = 1) +
    geom_point(data = df_noro, aes(x = time, y = I1), color = 'lightseagreen', size = 2) +
    labs(title = 'Fitted curve to Norovirus outbreak',
         x = 'Days since outbreak', y = 'Number of newly infected') +
    scale_x_continuous(breaks = seq(min(df_noro$time), max(df_noro$time), by = 1)) +
    theme_gray()
}


## plot residuals
plot_residuals <- function(df_res){
  ggplot(df_res, aes(x = t, y = residuals)) +
    geom_point(color = "lightseagreen") +
    labs(title = "Residuals Plot",
         x = "Time",
         y = "Residuals") +
    theme_gray()
}
