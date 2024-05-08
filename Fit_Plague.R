# =============================================================================
# TO DO: Fitting Plague 2017 Madagascar
# Author: Andrea Pinke
# =============================================================================

# Load required packages
library(tidyverse)
library(fitode)

# Set working directory
setwd("C:/Users/Andrea/OneDrive/Asztali gép/TUM/Mathematik/Master/Case Studies/Project/R")

# Load data
df_plague <- read.csv(file = "Data/data_plague_2017_madagascar.csv")
df_plague

# Introduce total infected number as pneumonic + bubonic cases
df_plague %>%
  mutate(total = pneumonic+bubonic) -> df_plague

df_plague

# Convert date to proper format for plotting
df_plague$date <- as.Date(df_plague$date, format = "%d/%m/%Y")
df_plague

# Plot total infected number by date
df_plague %>%
  ggplot(aes(x = date, y = total)) +
  geom_path(color = 4) +
  geom_point()


# Set up model --> from General Model code
fconst <- function(t){
  return(1)
}

model <- function(t, x, parameters, f=fconst){
  
  S <- x[1]
  E <- x[2]
  I1 <- x[3]
  I2 <- x[4]
  R <- x[5]
  V1 <- x[6]
  V2 <- x[7]
  D1 <- x[8]  
  D2 <- x[9]
  SV <- x[10]
  IV <- x[11]
  P <- x[12]
  
  f1 <- f(t)
  e1 <- parameters$e1
  e2 <- parameters$e2
  f14_3 <- parameters$f14_3
  f15_6 <- parameters$f15_6
  f17_4 <- parameters$f17_4
  fE <- parameters$fE
  fI1 <- parameters$fI1
  fI2 <- parameters$fI2
  fD <- parameters$fD
  fP <- parameters$fP
  fMuI2 <- parameters$fMuI2
  
  delta <- parameters$delta
  mu <- parameters$mu
  muI <- parameters$muI
  muI1 <- muI + mu
  muI2 <- fMuI2 * muI + mu
  deltaV <- parameters$deltaV
  muV <- parameters$muV
  alphaP <- parameters$alphaP
  muP <- parameters$muP
  a1 <- parameters$a1
  a2 <- parameters$a2
  a3 <- parameters$a3
  a4 <- parameters$a4
  a5 <- parameters$a5
  a6 <- parameters$a6
  b6 <- parameters$b6
  a7 <- parameters$a7
  a8 <- parameters$a8
  a9 <- (1-e2)*a1
  a10 <- (1-e1)*a1
  a11 <- (1-e1)*a5
  a12 <- (1-e2)*a5
  a13 <- parameters$a13
  a14 <- f14_3*a3
  a15 <- f15_6*a6
  a16 <- parameters$a16
  a17 <- f17_4*a4
  a18 <- parameters$a18
  a19 <- parameters$a19
  transHH <- parameters$transHH
  transHV <- parameters$transHV
  transVV <- parameters$transVV
  transVH <- parameters$transVH
  transPH <- parameters$transPH
  
  N <- S+E+I1+I2+R+V1+V2+D1+D2 #total human population N
  infH <- fE*E+fI1*I1+fI2*I2+fD*(fI1*D1+fI2*D2) #infective Population infH weighted by relative infectivities fE, FI1, fI2, fD
  phiH <- max(0,infH/N, na.rm = T)
  NV <- SV+IV #total vector population
  phiV <- max(0,IV/NV, na.rm = T)
  phiP <- max(0,fP*P/N, na.rm = T)
  
  #model code
  dSdt <- delta + a19*R - a1*(transHH*phiH+transVH*phiV+transPH*phiP)*S - a5*S*(transHH*phiH+transVH*phiV+transPH*phiP) - a7*S - mu*S
  dEdt <- a1*S*(transHH*phiH+transVH*phiV+transPH*phiP) + a10*V1*(transHH*phiH+transVH*phiV+transPH*phiP) + a9*V2*(transHH*phiH+transVH*phiV+transPH*phiP) - a2*E - mu*E
  dI1dt <- a2*E + a5*S*(transHH*phiH+transVH*phiV+transPH*phiP) + a11*V1*(transHH*phiH+transVH*phiV+transPH*phiP) + a12*V2*(transHH*phiH+transVH*phiV+transPH*phiP) - (a3+a6+a13)*I1 + b6*R  - (muI1)*I1
  dI2dt <- a3*I1 - (a4+a16)*I2 - (muI2)*I2
  dRdt <- a4*I2 + a6*I1 + a15*D1 + a17*D2 - b6*R - a19*R - mu*R
  dV1dt <- a7*S - a8*V1 - (a10+a11)*V1*(transHH*phiH+transVH*phiV+transPH*phiP) - mu*V1
  dV2dt <- a8*V1 - (a9+a12)*V2*(transHH*phiH+transVH*phiV+transPH*phiP) - mu*V2
  dD1dt <- a13*I1 -(a14+a15)*D1 - (muI1)*D1
  dD2dt <- a14*D1 + a16*I2 - a17*D2 - (muI2)*D2
  dSVdt <- deltaV - a18*SV*(transHV*phiH + transVV*phiV) - muV*SV
  dIVdt <- a18*SV*(transHV*phiH + transVV*phiV) - muV*IV
  dPdt <- alphaP*infH - muP*P
  
  
  #return list of diffs
  list(c(dSdt, dEdt, dI1dt, dI2dt, dRdt, dV1dt, dV2dt, dD1dt, dD2dt, dSVdt, dIVdt, dPdt))
}