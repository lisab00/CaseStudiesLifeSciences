# =============================================================================
# TO DO: Fitting Plague 2017 Madagascar outbreak data -- SEI1D1SvIv model
# Author: Andrea Pinke
# =============================================================================

# Load required packages
library(tidyverse)
library(deSolve)
library(FME)

# Load datasets
df_plague <- read.csv(file = "Data/data_plague_2017_madagascar.csv")
df_plague_relevant <- read.csv(file = "Data/relevant_data_plague_2017_madagascar.csv")


## Model definition

# define time-dependent parameter functions - adapted from plague paper

# seasonality
func1 <- function(t){
  return(1.15 + 0.08*sin(2*pi/12*t) + 0.1*cos(2*pi/12*t))
}

# consider only outbreak period since for us only that interval is relevant 
# --> months 8,9,10,11 (122 days in total, data for first 118 days)
func1_outbreak <- fluct1(seq(8,12, length.out = 122))


# delay for detected compartment
func2 <- function(t,tau){
  return(1/(1+exp(tau-t)))
}

# vector intervention
func3 <- function(t,theta){
  return(1-1/(1+exp(theta-t)))
}


# define terms for weighting
np <- sum(df_plague$pneumonic) # total number of pneumonic cases
nb <- sum(df_plague$bubonic) # total number of bubonic cases
fraction_p <- np/(np+nb) # fraction of pneumonic cases
fraction_b <- nb/(np+nb) # fraction of bubonic cases

# weighted muI, a2
weighted.muI <- fraction_p*0.34 + fraction_b*0.26
weighted.a2 <- fraction_p*0.29 + fraction_b*0.23

# Model
SEIDIvSv_model <- function (time, state, parms, signal1, signal2, signal3, ...) {
  f1 <- signal1[time+52] # seasonality function only for large epidemic wave - from 09/22/2017 = day 53
  tau <- 26 # value for best parameter fit
  f2 <- signal2(time, tau)
  f3 <- signal3(time, tau)
  with(as.list(c(parms, state)), {
    dS <- 2317 - a1*S*(((1*I1+0*D1)/(S+E+I1+D1))+(Iv/(Iv+Sv))) - 0.00001721*S
    dE <- a1*S*(((1*I1+0*D1)/(S+E+I1+D1))+(Iv/(Iv+Sv))) - weighted.a2*E - 0.00001721*E
    dI1 <- weighted.a2*E - f2*a13*I1 - (weighted.muI + 0.00001721)*I1
    dD1 <- f2*a13*I1 - (weighted.muI + 0.00001721)*D1
    dSv <- 10 - f1*f3*a18*(Iv/(Iv+Sv))*Sv - 0.01*Sv
    dIv <- f1*f3*a18*(Iv/(Iv+Sv))*Sv - 0.01*Iv
    list(c(dS, dE, dI1, dD1, dSv, dIv))
  })
}


## Choose initial parameter values, time interval, initial values for ode
parms <- c(a1=0.000005,a18=0.4,a13=1) 
state0 <-  c(S=25570895,E=0,I1=1,D1=0,Sv=990,Iv=10)
times <- df_plague_relevant$time

## Simulate and plot the model and the data with initials
out <- ode(state0, times, SEIDIvSv_model, parms, signal1 = func1_outbreak, signal2 = func2, signal3 = func3)
plot(out) # plot model outputs for state variables
plot(out, obs=df_plague_relevant, obspar=list(pch=16, col="red")) # plot I1 against data


## Fit

# Define a cost function
cost <- function(p) {
  pp <- p[c("a1","a18","a13")]
  out <- ode(state0, times , SEIDIvSv_model, pp, signal1 = func1_outbreak, signal2 = func2, signal3 = func3)
  modCost(model = out[ , c(1,4)], obs = df_plague_relevant, method = "Marq")
}

# Set boundary constraints (all positive) for parameters and fit the model
fit <- modFit(f = cost, p = parms, lower = c(0,0,0))
summary(fit)
coef(fit)

# Plot fitted model with data
pars_est <- coef(fit)[c("a1","a13","a18")]
out_fit <- ode(state0, times , SEIDIvSv_model, pars_est, signal1 = func1_outbreak, signal2 = func2, signal3 = func3)
plot(out_fit, obs=df_plague_relevant, obspar=list(pch=16, col="red"), xlab = "Number of days since relevant start of outbreak", ylab = "Total number of infected cases")


# Check cost value for the model
cost <- modCost(model = out_fit[ , c(1,4)], obs = df_plague_relevant, weight = "none")
cost_df <- cost["residuals"]$residuals

# Plot residuals --> should be cca. normal distributed for a good fit
cost_df %>%
  ggplot(aes(x = x, y = res.unweighted)) +
  geom_path(color = 4) +
  geom_point() +
  xlab("Number of Days") +
  ylab("Unweighted Residuals")

n <- length(times) # number of data points
SSR <- sum(cost_df$res.unweighted^2) # sum of squared residuals
RSE <- sqrt(SSR/(n-2)) # residual standard error
SSR; RSE
