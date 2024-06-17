# =============================================================================
# TO DO: Fitting Plague 2017 Madagascar
# Author: Lisa Beer
# =============================================================================

# Load required packages
library(tidyverse)
library(fitode)
library(deSolve)
library(data.table)

# Set working directory
setwd("C:/Users/lihel/Documents/master/so24/case_studies/CaseStudiesLifeSciences/")

# Load data
df_plague <- read.csv(file = "data_plague_2017_madagascar.csv")
df_plague

# Introduce total infected number as pneumonic + bubonic cases
mutate(df_plague, total = pneumonic+bubonic) -> df_plague

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


#initial params estimated
delta <- 32.48 #google birth rate
mu <- 0.0063 #per capita death rate S
muI <- 0.08 #per capita death rate I (calculated from data)
a1 <- 1 #pick anything
a2 <- 1
a5 <- 1
a6 <- 1
a19 <- 1
b6 <- 1
fI1 <- 1
transHH <- T
S <- 26000000 #Madagascar population 2017
E <- 0
I1 <- 1
R <- 0

initial_params <- list(delta=delta, mu=mu, muI=muI, fMuI2=0, alphaP=0, muP=0, deltaV=0, muV=0, a1=a1, a2=a2, a3=0,
                       a4=0, a5=a5, a6=a6, b6=b6, a7=0, a8=0, a13=0, a16=0, a18=0, a19=a19, e1=0, e2=0,
                       f14_3=0, f15_6=0, f17_4=0, fE=0, fI1=fI1, fI2=0, fD=0, fP=0, transHH=T, transHV=F, transVV=F, transVH=F, transPH=F)
initial_sir <- c(S=S, E=E, I1=I1, I2=0, R=0, V1=0, V2=0, D1=0, D2=0, SV=0, IV=0, P=0)
params_to_optim <- which(initial_params != 0) #define which params i want to fit

#define time line of epidemic
df_plague$date <- as.POSIXct(df_plague$date)
start_date <- min(df_plague$date)
df_plague %>% 
  mutate(days = as.numeric(difftime(df_plague$date, start_date, units = "days"))) -> df_plague
times <- df_plague$days

#first simulate with initial params
simulation <- as.data.table(
  ode(func = model,
      y = initial_sir,
      times = times,
      parms = initial_params)
)

df_sim <- gather(simulation, key = Compartment, value = count, c("S", "E", "I1", "R"))
df_sim$Compartment <- factor(df_sim$Compartment, levels = c("S", "E", "I1", "R"))
ggplot(df_sim, aes(x = time, y=count, group = Compartment, color = Compartment)) + geom_line(size=1)


#try to fit

#define objective function for residual sum of squares
objective <- function(initial_params, observed_data){
  initial_params <- as.list(initial_params) #convert initial_params back to list format (optim forces them to numeric)
  simulation <- ode(func = model,
                    y = initial_sir,
                    times = times,
                    parms = initial_params)
  simulated_data <- simulation[,"I1"]
  residual_sum_sq <- sum((simulated_data - observed_data)^2)
  return(residual_sum_sq)
}

objective <- function(params_to_optim, observed_data){
  simulation_params <- #somehow merge initial params with optimized new ones #convert initial_params back to list format (optim forces them to numeric)
  simulation <- ode(func = model,
                    y = initial_sir,
                    times = times,
                    parms = simulation_params)
  simulated_data <- simulation[,"I1"]
  residual_sum_sq <- sum((simulated_data - observed_data)^2)
  return(residual_sum_sq)
}

#optimize objective
results <- optim(par = initial_params[params_to_optim], fn = objective, observed_data = data.table(df_plague$total))
#results <- optim(par = initial_params, fn = objective, observed_data = data.table(df_plague$total))

#plot results
simulation <- as.data.table(
  ode(func = model,
      y = initial_sir,
      times = times,
      parms = as.list(results$par))
)

df_sim <- gather(simulation, key = Compartment, value = count, c("S", "E", "I1", "R"))
df_sim$Compartment <- factor(df_sim$Compartment, levels = c("S", "E", "I1", "R"))
ggplot(df_sim, aes(x = time, y = count, group = Compartment, color = Compartment)) + geom_line(size=1)
