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

# Introduce total infected number as pneumonic + bubonic cases
mutate(df_plague, total = pneumonic+bubonic) -> df_plague

# Convert date to proper format for plotting
df_plague$date <- as.Date(df_plague$date, format = "%d/%m/%Y")

#define time line of epidemic
df_plague$date <- as.POSIXct(df_plague$date)
start_date <- min(df_plague$date)
df_plague %>% 
  mutate(days = as.numeric(difftime(df_plague$date, start_date, units = "days"))) -> df_plague
times <- df_plague$days

# Set up model
model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS <- delta - a1*(1*fI1*I1/(S+E+I1))*S - a5*S*(1*fI1*I1/(S+E+I1)) - mu*S + a19*R
    dE <- a1*S*(1*fI1*I1/(S+E+I1)) - a2*E - mu*E
    dI1 <- a2*E + a5*S*(1*fI1*I1/(S+E+I1)) - (muI)*I1 + b6*R - a6*I1
    dR <- a6*I1 - b6*R - a19*R - mu*R
    return(list(c(dS, dI1, dE, dR)))
  })
}

# #define intitial parameters
initial_params <- list(delta=32.48, mu=0.0063, muI=0.08, a1=1, a2=1, a5=1, a6=1,
                       b6=1, a19=1, fI1=1)
initial_sir <- c(S=26000, E=0, I1=1, R=0)

#first simulate with initial params
simulation <- ode(func = model, y = initial_sir, times = times, parms = initial_params)

# simulation <- as.data.table(
#   ode(func = model,
#       y = initial_sir,
#       times = times,
#       parms = initial_params)
# )

df_sim <- gather(as.data.table(simulation), key = Compartment, value = count, c("S", "E", "I1", "R"))
df_sim$Compartment <- factor(df_sim$Compartment, levels = c("S", "E", "I1", "R"))
ggplot(df_sim, aes(x = time, y=count, group = Compartment, color = Compartment)) + geom_line(size=1)
