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

# fit from 53th day on since the beginning of the outbreak --> same as in paper
df_plague_relevant <- filter(df_plague_relevant, time>53)

# re-scale time for ode solver to start with day 1
df_plague_relevant$time <- df_plague_relevant$time - 53
df_plague_relevant

# =============================================================================
### Fitting

# fluctuations with function and estimated pars from Plague paper
# A = 1.15, B = 0.08, C = 0.1
# A + B*sin(2*pi/12*t) + C*cos(2*pi/12*t)
fluct <- function(t){
  return(1.15 + 0.08*sin(2*pi/12*t) + 0.1*cos(2*pi/12*t))
}

# fluctuation only in the outbreak period --> months 8,9,10,11 by day (122 days in total)
fluct_outbreak <- fluct(seq(8,12, length.out = 122))
plot(fluct_outbreak, xlab = "Number of days since outbreak")


## Model --> adjust parameter values first
# set fI1 = 1 since only one stage of infection
# set sigmas 0 or 1 based on whether type of transmission is possible
# total population number from Plague paper: 25 570 895
# fix constant birth term: https://www.macrotrends.net/global-metrics/countries/MDG/madagascar/birth-rate
# in 2017: 33 078 births per 1000 people --> total births per day= (25570.895*33078)/365 = 2317
# fix death rate: https://www.macrotrends.net/global-metrics/countries/MDG/madagascar/death-rate
# in 2017: 6.282 deaths per 1000 people --> time-independent per capita death rate= 6.282/25570.895 = 0.00024567
# other possibilities: 0.017 (death rate per day) or 0.00000067 (death rate per one person per day) 

SEIf_model <- function (time, state, parms, signal, ...) {
  f <- signal[time+53] # relevant period for fitting is from 09/23/2017 = day 54
  with(as.list(c(parms, state)), {
    dS <- 2317 - a1*f*(I1/(S+E+I1))*S - 0.00024567*S
    dE <- a1*f*S*(I1/(S+E+I1)) - a2*E - 0.00024567*E
    dI1 <- a2*E - (muI + 0.00024567)*I1
    list(c(dS, dE, dI1))
  })
}
# Note: without vector compartments the transition rate a1 implicitly includes the vector population density


## Choose initial parameter values, time interval, initial values for ode

# good found parameters with SSR: 7893.71, RSE: 11.19
parms <- c(E=120, a1=0.01, a2=0.1, muI=0.07)
state0 <-  c(S=25570895,E=120,I1=5)
times <- df_plague_relevant$time


## Simulate and plot the model and the data
out <- ode(state0, times , SEIf_model, parms[2:4], signal = fluct_outbreak)
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
  #state0[2] <- p["E"]
  pp <- p[c("a1","a2","muI")]
  out <- ode(state0, times, SEIf_model, pp, signal = fluct_outbreak)
  modCost(model = out[ , c(1,4)], obs = df_plague_relevant, method = "Marq")
}

# Set boundary constraints for parameters and fit the model
# incubation time from Plague paper: 1-12 days --> a2(=inverse incubation time) in [1/12,1]

#fit <- modFit(f = cost, p = parms, lower=c(0,0,0,0), upper=c(25570895,1e5, 1e5, 1e5))
fit <- modFit(f = cost, p = parms[2:4], lower=c(0,1/12,0))
summary(fit)
coef(fit)

## Plot fitted model and data
state0_est <- state0
#state0_est[2] <- coef(fit)["E"]
pars_est <- coef(fit)[c("a1","a2","muI")]
out_fit <- ode(state0_est, times, SEIf_model, pars_est, signal = fluct_outbreak)

# scale back time
#out_df <- data.frame(out_fit)
#out_df$time <- as.Date(out_df$time+17378+53, origin = "1970-01-01")
#out_df

plot(out_fit, obs=df_plague_relevant, obspar=list(pch=16, col="red"), xlab = "Number of days since relevant start of outbreak", ylab = "Total number of infected cases")
