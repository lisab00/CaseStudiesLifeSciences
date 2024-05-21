# =============================================================================
# TO DO: Fitting Plague 2017 Madagascar
# Author: Andrea Pinke, Lisa Beer
# =============================================================================


# Load required packages
library(tidyverse)
library(deSolve)
library(FME)

#set working directory and load data
#setwd("C:/Users/Andrea/OneDrive/Asztali g?p/TUM/Mathematik/Master/Case Studies/Project/R")
setwd("C:/Users/lihel/Documents/master/so24/case_studies/CaseStudiesLifeSciences/")
df_plague_relevant <- read.csv(file = "data_plague_relevant.csv")


# =============================================================================
### Preparations

# fluctuations with function and estimated pars from Plague paper
# A = 1.15, B = 0.08, C = 0.1
# A + B*sin(2*pi/12*t) + C*cos(2*pi/12*t)
fluct <- function(t){
  return(1.15 + 0.08*sin(2*pi/12*t) + 0.1*cos(2*pi/12*t))
}

## Model --> adjust parameter values first
# set fI1 = 1 since only one stage of infection
# set sigmas 0 or 1 based on whether type of transmission is possible
# total population number from Plague paper: 25 570 895
# fix constant birth term: https://www.macrotrends.net/global-metrics/countries/MDG/madagascar/birth-rate
# in 2017: 33 078 births per 1000 people --> total births per day= (25570.895*33078)/365 = 2317
# https://www.macrotrends.net/global-metrics/countries/MDG/madagascar/death-rate
# in 2017: 6.282 deaths per 1000 people
# 0.017 (death rate per day) or 0.00000067 (death rate per one person per day) ?
# other idea (Lisa): 6.3 / 1000 = 0.0063 death rate per person (time-independent)

SEIf_model <- function (time, state, parms, ...) {
  f <- fluct(time)
  with(as.list(c(parms, state)), {
    dS <- 2317 - a1*f*(I1/(S+E+I1))*S - 0.0063*S
    dE <- a1*f*S*(I1/(S+E+I1)) - a2*E - 0.0063*E
    dI1 <- a2*E - (muI + 0.0063)*I1
    list(c(dS, dE, dI1))
  })
}
# Note: without vector compartments the transition rate a1 implicitly includes the vector population density


## Choose initial parameter values, time interval, initial values for ode
# found parameters with SSR: 7686, RSE: 10.86
times <- df_plague_relevant$time
state0 <-  c(S=25570895,E=150,I1=1)
parms <- c(a1=0.01, a2= 0.1, muI=0.08) #geben gutes Resultat mit E = 150
#parms <- c(a1=0.09, a2=0.02, muI=0.08)

## Simulate and plot the model and the data
out <- ode(state0, times , SEIf_model, parms)
plot(out)
plot(out, obs=df_plague_relevant, obspar=list(pch=16, col="red"))


# =============================================================================
### Fitting

#define cost function
cost <- function(p) {
  pp <- p[c("a1","a2","muI")]
  out <- ode(state0, times, SEIf_model, pp)
  modCost(model = out[ , c(1,4)], obs = df_plague_relevant, weight = "none"
          , method = "Marq")
}

#define weighted cost function
#weight the residuals for fitting: make fitting more sensitive to mode
#algo weights with 1/res_weights_inv
df_plague_relevant %>%
  mutate(res_weights_inv=1/out[,"I1"]^1) -> df_plague_relevant_weighted #power can be adjusted

cost_weighted <- function(p) {
  pp <- p[c("a1","a2","muI")]
  out <- ode(state0, times, SEIf_model, pp)
  modCost(model = out[ , c(1,4)], obs = df_plague_relevant_weighted
          , method = "Marq", err="res_weights_inv") #res are weighted by 1/err
}

# set boundary constraints for parameters and fit the model
#incubation time from Plague paper: 1-12 days --> a2(=inverse incubation time) in [1/12,1]
control_list <- list(nprint=1 #print iterations
                     , factor=1 #determine initial step bound
                     #, maxiter=25 #set max number of iterations performed
)

#fit with normal cost
fit <- modFit(f = cost, p = parms, lower=c(0,0.08,0), upper=c(1e5, 1, 1e5), 
              control=control_list)
summary(fit)

#fit with weighted cost
fit_weighted <- modFit(f = cost_weighted, p = parms, lower=c(0,0.08,0), 
                       upper=c(1e5, 1, 1e5), control=control_list)
summary(fit_weighted)
coef(fit_weighted)
coef(fit)

# =============================================================================
### Plots

#simulation with normal fit
out_fit <- ode(state0, times, SEIf_model, coef(fit), rtol=10e-5, atol=10e-3)

#weighted with weighted fit
out_fit_weighted <- ode(state0, times, SEIf_model, coef(fit_weighted),
                        rtol=10e-5, atol=10e-3)

#plot
plot(out_fit, obs=df_plague_relevant, obspar=list(pch=16, col="red"), col="blue",
     lwd=3, xlab = "Number of days since outbreak",
     ylab = "Total number of infected cases")
lines(out[,"I1"])
lines(out_fit_weighted[,"I1"], col="green", lwd=2)
legend("topright", legend=c("Initial", "Fit", "Weighted Fit"),
       col=c("black", "blue", "green"), lty=c(1,1,1), lwd=c(2,2,2))

# =============================================================================
### Compare

#compare unweighted
modCost(model = out_fit[ , c(1,4)], obs = df_plague_relevant, weight = "none")$var
modCost(model = out[ , c(1,4)], obs = df_plague_relevant, weight = "none")$var

#compare weighted
modCost(model = out_fit_weighted[ , c(1,4)], obs = df_plague_relevant_weighted, 
        weight = "none", err="res_weights_inv")$var
modCost(model = out[ , c(1,4)], obs = df_plague_relevant_weighted, 
        weight = "none", err="res_weights_inv")$var

