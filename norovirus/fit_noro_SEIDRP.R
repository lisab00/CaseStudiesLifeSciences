setwd("C:/Users/lihel/Documents/master/so24/case_studies/CaseStudiesLifeSciences/")

library(readxl)
library(tidyverse)
library(deSolve)
library(FME)
library(ggplot2)

# =============================================================================
# Load Data

ti <- 6   # time of health interventions, outbreak11
#ti <- 4    # outbreak2
#ti <- 3    # outbreak6

sheet <- "outbreak11"    # adapt to correct outbreak

# consider outbreaks with pathogen compartment
df_noro <- as.data.frame(read_excel("norovirus/relevantData.xlsx", sheet=sheet))
plot(df_noro$I1)
abline(v=ti, col="red")
text(ti, max(df_noro$I1), label ="Health Intervention", pos=4, col="red")


# add column with weights to df
'df_noro <- mutate(df_noro, weights_inv=1/df_noro$I1)
df_noro$weights_inv[is.infinite(df_noro$weights_inv)] <- 1  # replace Inf values
'

## functions for health interventions
# cleaning measures at time t=ti reduce virus appearance in water to effP
signalP <- function(t, effP){
  if (t < ti){
    return(1)
  } else {
    return(effP)
  }
}

# quarantine measures, signalD resembles the isolation coefficient q
signalD <- function(t, q){
  if (t < ti){
    return(0)
  } else {
    return(q)
  }
}


# plot intervention functions
t <- df_noro$time
plot(t, sapply(t, effP=0.3, signalP), ylim=c(0,1), type="l")
lines(t, sapply(t, q=1, signalD), col="red")
legend("topright", legend=c("signalP", "signalD"), col=c("black", "red"), lwd=2)


## set up model
model <- function(time, state, parms, signalP, signalD, ...){
  signalP <- signalP(time, effP)
  signalD <- signalD(time, q)
  with(as.list(c(parms, state)), {
    dS <- -a1 * (fI1 * I1 + fP * P) / N * S
    dE <- a1 * (fI1 * I1 + fP * P) / N * S - a2 * E
    dI1 <- a2 * E - ((1-signalD) * a6 + signalD * a13) * I1
    dD1 <- signalD * a13 * I1 - a15 * D1
    dR <- (1-signalD) * a6 * I1 + a15 * D1
    dP <- alpha_p * signalP * fI1 * I1 - mu_p * P
    list(c(dS, dE, dI1, dD1, dR, dP))
  })
}

## set fixed parameters from paper
a1 <- 1                        
a2 <- 1                         # cf omega (W)
a6 <- 0.3*0.03846 + 0.7*0.3333  # recovery rate weighted (0.3*symptomatic+0.7*asymptomatic)
a13 <- 2                        # isolation rate
a15 <- 0.3333                   # recovery rate of detected (only symptomatic)
mu_p <- 0.1                     # cf epsilon

# fix more parameters (can be varied)
q <- 1                          # rate of infected being sent to quarantine (0.6)
effP <- 0.1                     # fraction of pathogen getting into water after disinfection (0.1)
alpha_p <- 1                    # rate of pathogen appearance in water
#fP <- 1

# rmk: this leaves parms c(fI1, fP) to be fitted

# state0
N <- 1751      # outbreak11
#N <- 14500     # outbreak2
#N <- 3142      # outbreak6

state0 <- c(S=N-df_noro$I1[1],
            E=df_noro$I1[1]/(1-0.3), # 0.3 is proportion asymptomatic
            #E <- 8,                  # try different E0
            I1=df_noro$I1[1],
            D1=0, R=0, P=0)

# initial parms
parms <- c(fI1=0.6          # 0.6 (outb11)
           , fP=0.9         # 0.9 (outb11)
           )

## simulation and plot
timeline <- seq(df_noro$time[1], df_noro$time[length(df_noro$time)], 0.001)
out <- ode(state0, timeline, model, parms, signalP=signalP, signalD=signalD)
plot(out)
df_out <- as.data.frame(out)
ggplot() +
  geom_line(data = df_out, aes(x = time, y = I1), color = 'blue', size = 1) +
  geom_point(data = df_noro, aes(x = time, y = I1), color = 'red', size = 2) +
  labs(title = 'Infected', x = 'Days', y = 'Number of Individuals')


## define cost
cost <- function(p) {
  out <- ode(state0, timeline, model, p, signalP=signalP, signalD=signalD)
  modCost(model = out[ , c("time","I1")]
          , obs = df_noro, method = "Marq"
          #, err="weights_inv"
          )
}


## fitting
fit <- modFit(f=cost, p=parms, lower=c(0,0)
              , upper=c(1,1)
              )
summary(fit)
coef(fit)


## plot fitted model and data
pars_est <- coef(fit)
out_fit <- ode(state0, timeline, model, pars_est, signalP=signalP, signalD=signalD)
#plot(out_fit)
df_out_fit <- as.data.frame(out_fit)
ggplot() +
  geom_line(data = df_out_fit, aes(x = time, y = I1), color = 'blue', size = 1) +
  geom_point(data = df_noro, aes(x = time, y = I1), color = 'red', size = 2) +
  labs(title = 'Infected', x = 'Days', y = 'Number of Individuals') 

