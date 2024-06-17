# =============================================================================
# TO DO: Fitting Plague 2017 Madagascar
# Author: Andrea Pinke
# =============================================================================

# Load required packages
library(tidyverse)
library(deSolve)
library(FME)

# Set working directory
setwd("C:/Users/Andrea/OneDrive/Asztali gép/TUM/Mathematik/Master/Case Studies/Project/R")

# Load data
df_plague <- read.csv(file = "Data/data_plague_2017_madagascar.csv")
df_plague

# Introduce total infected number as pneumonic + bubonic cases
df_plague %>%
  mutate(total = pneumonic+bubonic) -> df_plague

df_plague_relevant <- select(df_plague, c("date", "total"))

# Convert date to proper format for plotting
df_plague_relevant$date <- as.Date(df_plague_relevant$date, format = "%d/%m/%Y")

# as.numeric gives number of days since "1970-01-01"
df_plague_relevant$date <- as.numeric(df_plague_relevant$date)

# rename columns for modCost function --> same name as ode output
df_plague_relevant <- rename(df_plague_relevant, c(time = date, I1 = total))

# re-scale time for ode solver to start with day 1
df_plague_relevant$time <- df_plague_relevant$time - 17378

df_plague_relevant


##
# fit from 53th day on since the beginning of the outbreak --> same as in paper
df_plague_relevant <- filter(df_plague_relevant, time>52)

# re-scale time for ode solver to start with day 1
df_plague_relevant$time <- df_plague_relevant$time - 52
df_plague_relevant

# =============================================================================
### Fitting

# fluctuations with function and estimated pars from Plague paper
# A = 1.15, B = 0.08, C = 0.1
# A + B*sin(2*pi/12*t) + C*cos(2*pi/12*t)
fluct1 <- function(t){
  return(1.15 + 0.08*sin(2*pi/12*t) + 0.1*cos(2*pi/12*t))
}

# fluctuation only in the outbreak period --> months 8,9,10,11 by day (122 days in total)
fluct1_outbreak <- fluct1(seq(8,12, length.out = 122))
plot(fluct1_outbreak, xlab = "Number of days since outbreak")

# intervention effects
fluct2 <- function(t,theta,phi){
  return(1-1/(1+phi+exp(theta-t)))
}

# basic numbers for weighting
np <- sum(df_plague$pneumonic)
nb <- sum(df_plague$bubonic)
fraction_p <- np/(np+nb)
fraction_b <- nb/(np+nb)
np; nb; fraction_b; fraction_p

# estimated parameter values from Plague paper
phi_est <- 0.11
theta.b <- 17.93
theta.p <- 8.89
theta.weighted <- fraction_b*theta.b + fraction_p*theta.p
theta.weighted

# outbreak lasts 66 days
fluct2_outbreak <- fluct2(1:66, theta.weighted, phi_est)
plot(fluct2_outbreak, xlab = "Number of days since relevant start of outbreak")


## Model --> adjust parameter values first
# set fI1 = 1 since only one stage of infection
# set sigmas 0 or 1 based on whether type of transmission is possible
# total population number from Plague paper: 25 570 895
# fix constant birth term: https://www.macrotrends.net/global-metrics/countries/MDG/madagascar/birth-rate
# in 2017: 33.078 births per 1000 people --> total births per day= (25570.895*33.078)/365 = 2317
# fix death rate: https://www.macrotrends.net/global-metrics/countries/MDG/madagascar/death-rate
# in 2017: 6.282 deaths per 1000 people --> time-independent per capita death rate= (6.282/1000)/365 = 0.00001721

# weighted muI: separate bubonic/pneumonic values from Plague paper
weighted.muI <- fraction_p*0.34 + fraction_b*0.26
weighted.muI

weighted.a2 <- fraction_p*0.29 + fraction_b*0.23
weighted.a2

SEIf_model <- function (time, state, parms, signal1, signal2, ...) {
  f1 <- signal1[time+52] # relevant period for fitting is from 09/22/2017 = day 53
  f2 <- signal2(time, parms["theta"], parms["phi"])
  #f <- signal[time]
  with(as.list(c(parms, state)), {
    dS <- 2317 - a1*f1*f2*S*(I1/(S+E+I1)) - 0.00001721*S
    dE <- a1*f1*f2*S*(I1/(S+E+I1)) - weighted.a2*E - 0.00001721*E
    dI1 <- weighted.a2*E - (weighted.muI + 0.00001721)*I1
    list(c(dS, dE, dI1))
  })
}
# Note: without vector compartments the transition rate a1 implicitly includes the vector population density



## Choose initial parameter values, time interval, initial values for ode
parms <- c(a1=2,theta=theta.weighted,phi=phi_est)
# --> gives estimated parms: c(a1=2.27982943, theta=9.17602432, phi=0.09526503) # SSR=2836.42, RSE=6.657256
state0 <-  c(S=25570895,E=0,I1=1)
times <- df_plague_relevant$time


## Simulate and plot the model and the data
out <- ode(state0, times , SEIf_model, parms, signal1 = fluct1_outbreak, signal2 = fluct2)
plot(out) # plot S,E,I1 solutions

plot(out, obs=df_plague_relevant, obspar=list(pch=16, col="red")) # plot I1 against data

# Check cost value for the model
cost <- modCost(model = out[ , c(1,4)], obs = df_plague_relevant, weight = "none")
cost
cost_df <- cost["residuals"]$residuals

## Plot residuals --> should be cca. normal distributed for a good fit
cost_df %>%
  ggplot(aes(x = x, y = res.unweighted)) +
  geom_path(color = 4) +
  geom_point() +
  xlab("Number of Days") +
  ylab("Unweighted Residuals")

n <- length(times) # number of data points
SSR <- sum(cost_df$res.unweighted^2) # sum of squared residuals
RSE <- sqrt(SSR/(n-2))# residual standard error
SSR; RSE

## Fit

# Define a cost function
cost <- function(p) {
  pp <- p[c("a1", "theta", "phi")]
  out <- ode(state0, times , SEIf_model, pp, signal1 = fluct1_outbreak, signal2 = fluct2)
  modCost(model = out[ , c(1,4)], obs = df_plague_relevant, method = "Marq")
}

# Set boundary constraints (all positive) for parameters and fit the model
fit <- modFit(f = cost, p = parms, lower = c(0,0,0))
summary(fit)
coef(fit)

## Plot fitted model and data
pars_est <- coef(fit)[c("a1", "theta", "phi")]
out_fit <- ode(state0, times , SEIf_model, pars_est, signal1 = fluct1_outbreak, signal2 = fluct2)

plot(out_fit, obs=df_plague_relevant, obspar=list(pch=16, col="red"), xlab = "Number of days since relevant start of outbreak", ylab = "Total number of infected cases")


# Plot for poster
out_df <- data.frame(out_fit)
df_plague_relevant %>%
  ggplot(aes(x = time, y = I1)) +
  geom_point(color="lightseagreen", size=2) +
  labs(title="Fit curve to data",
       subtitle="Best parameter fit",
       x="days since outbreak",
       y="infected individuals") +
  theme_gray(base_size=15) +
  geom_line(aes(x=time, y=out_df$I1), size=1) -> plot3

ggsave(plot = plot3, "fitted_curve.png", width = 6, height = 5)
