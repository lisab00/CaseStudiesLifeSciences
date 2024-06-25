setwd("C:/Users/lihel/Documents/master/so24/case_studies/CaseStudiesLifeSciences/")

library(readxl)
library(tidyverse)
library(deSolve)
library(FME)
library(ggplot2)

# =============================================================================
# Load Data
# consider outbreak 11 (with pathogen compartment)
df_noro <- as.data.frame(read_excel("norovirus/relevantData.xlsx", sheet="outbreak11"))
plot(df_noro)
abline(v=6, col="red")
text(6, max(df_noro$I1), label ="Health Intervention", pos=4, col="red")


## functions for health interventions
signalP <- function(t, effP){
  if (t < 6){
    return(1)
  } else {
    return(effP)
  }
}

## set up model
model <- function(time, state, parms, signalP, ...){
  signalP <- signalP(time, parms["effP"])
  with(as.list(c(parms, state)), {
    dS <- -a1 * (fI1 * I1 + fP * P) / N * S
    dE <- a1 * (fI1 * I1 + fP * P) / N * S - a2 * E
    dI1 <- a2 * E - a6 * I1
    dR <- a6 * I1
    dP <- alpha_p * signalP * fI1 * I1 - mu_p * P
    list(c(dS, dE, dI1, dR, dP))
  })
}

## set (fixed) parameters
fI1 <- 1                                # only 1 stage of infections
a2 <- 1                                 # cf omega (W)
a6 <- 0.3*0.03846 + 0.7*0.3333          # recovery rate weighted (symp, asymp)
mu_p <- 0.1                             # cf epsilon

# rmk: this leaves parms c(alpha_p, a1, fP, effP) to be fitted

# state0
N <- 1751
state0 <- c(S=N-df_noro$I1[1],
            E=df_noro$I1[1]/(1-0.3), #0.3 is proportion asymptomatic
            I1=df_noro$I1[1],
            R=0, P=0)

# initial parms
parms <- c(a1=0.1, fP=0.3, alpha_p=0.3, effP=0.8)

## simulation and plot
out <- ode(state0, df_noro$time, model, parms, signalP=signalP)
modCost(model = out[ ,c("time","I1")], obs = df_noro, method = "Marq")

plot(out)
df_out <- as.data.frame(out)

ggplot() +
  geom_line(data = df_out, aes(x = time, y = I1), color = 'blue', size = 1) +
  geom_point(data = df_noro, aes(x = time, y = I1), color = 'red', size = 2) +
  labs(title = 'Infected', x = 'Days', y = 'Number of Individuals')


## define cost
cost <- function(p) {
  out <- ode(state0, df_noro$time, model, p, signalP=signalP)
  modCost(model = out[ ,c("time","I1")], obs = df_noro, method = "Marq")
}


## fitting
fit <- modFit(f=cost, p=parms, lower=c(0,0,0,0), upper=c(1,1,1,1))
summary(fit)
coef(fit)


## plot fitted model and data
pars_est <- coef(fit)
out_fit <- ode(state0, df_noro$time, model, pars_est, signalP=signalP)
df_out_fit <- as.data.frame(out_fit)

ggplot() +
  geom_line(data = df_out_fit, aes(x = time, y = I1), color = 'blue', size = 1) +
  geom_point(data = df_noro, aes(x = time, y = I1), color = 'red', size = 2) +
  labs(title = 'Infected', x = 'Days', y = 'Number of Individuals')

