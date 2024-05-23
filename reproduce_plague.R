# =============================================================================
# TO DO: Reproduce Plague 2017 Madagascar
# Author: Lisa Beer
# =============================================================================

#load required packages
library(deSolve)
library(tidyverse)

#set working directory
setwd("C:/Users/lihel/Documents/master/so24/case_studies/CaseStudiesLifeSciences/")
df_plague <- read.csv(file = "data_plague_2017_madagascar.csv")
df_plague_relevant <- read.csv(file="data_plague_relevant.csv")

#set params, state0 and time scale
#copy estimated parms from paper
N <- 25570895
p <- 1e-4 #Sb = p*N, proportion of N exposed to ratfleas

parms <- c(alpha=1.9e-3, beta=2.23, gammab=0.23, gammap=0.29, deltab=0.26,
           deltap=0.34, epsilon=0.03)
times <- df_plague_relevant$time
state0 <- c(Sb=N*p, Sp=N*(1-p), Eb=0, Ep=0, Ib=1, Ip=1)

#define fluctuant functions
theta <- 0.11
taoh <- 8.89
taof <- 17.93
A <- 1.15
B <- 0.08
C <- 0.1

f_irf <- function(t){return(A + B*sin(2*pi/12*t) + C*cos(2*pi/12*t))}
f_itvf <- function(t){return(1-1/(1+theta+exp(taof-t)))}
f_itvh <- function(t){return(1-1/(1+theta+exp(taoh-t)))}

#set up model
model <- function (time, state, parms, ...) {
  f_irf <- f_irf(time)
  f_itvf <- f_itvf(time)
  f_itvh <- f_itvh(time)
  with(as.list(c(parms, state)), {
    dSb <- -Sb*alpha*f_irf*f_itvf - beta*Sb*(Ip/N)*f_itvh
    dSp <- -beta*Sp*(Ip/N)*f_itvh
    dEb <- Sb*alpha*f_irf*f_itvf - gammab*Eb
    dEp <- beta*(Sb+Sp)*(Ip/N)*f_itvh - gammap*Ep
    dIb <- gammab*Eb - epsilon*Ib - deltab*Ib
    dIp <- gammap*Ep + epsilon*Ib - deltap*Ip
    list(c(dSb, dSp, dEb, dEp, dIb, dIp))
  })
}

#simulation
out <- ode(state0, times , model, parms)

#plot
plot(out)

#I1=Ib+Ip
I1 <- out[,"Ib"]+out[,"Ip"]
plot(df_plague_relevant)
lines(I1)
