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
df_plague_relevant

# as.numeric gives number of days since "1970-01-01"
df_plague_relevant$date <- as.numeric(df_plague_relevant$date)

# rename columns for modCost function --> same name as ode output
df_plague_relevant <- rename(df_plague_relevant, c(time = date, I1 = total))
df_plague_relevant

# Plot total infected number by date
df_plague %>%
  ggplot(aes(x = date, y = total)) +
  geom_path(color = 4) +
  geom_point()

# =============================================================================
# Fitting

# Model
SEI_model <- function (time, state, parms, ...) {
  with(as.list(c(parms, state)), {
    dS <- delta - a1*(1*fI1*I1/(S+E+I1))*S - a5*S*(1*fI1*I1/(S+E+I1)) - mu*S
    dE <- a1*S*(1*fI1*I1/(S+E+I1)) - a2*E - mu*E
    dI1 <- a2*E + a5*S*(1*fI1*I1/(S+E+I1)) - (muI1)*I1
    list(c(dS, dE, dI1))
  })
}

# re-scale time for ode solver to start with day 1
df_plague_relevant$time <- df_plague_relevant$time - 17378
df_plague_relevant

# choose initial parameter values, time interval, initial values for ode --> decide based on simulation
# intial total population number from Plague paper
parms <- c(delta=50000,a1=0.6,a2=0.6,a5=0.7,fI1=0.7,mu=0.1,muI1=0.4)
times <- df_plague_relevant$time
state0 <-  c(S=25570895,E=10,I1=1)

# simulate and plot the model and the data --> try to have shape of real data (unimodal)
out <- ode(state0, times, SEI_model, parms)
plot(out, obs=df_plague_relevant, obspar=list(pch=16, col="red"))

# Fitting
# define a cost function
cost <- function(p) {
  out <- ode(state0, times, SEI_model, p)
  modCost(model = out[ , c(1,4)], obs = df_plague_relevant, weight = "none")
}

# fit the model --> initial parameters taken from above, all positive
fit <- modFit(f = cost, p = parms, lower = 0)
summary(fit)
coef(fit)

# Plot fitted models and data
out_fit <- ode(state0, times, SEI_model, coef(fit))

# scale back time
#out_df <- data.frame(out_fit)
#out_df$time <- as.Date(out_df$time+17378, origin = "1970-01-01")

plot(out_fit, obs=df_plague_relevant, obspar=list(pch=16, col="red"), xlab = "Number of days since outbreak", ylab = "Total number of infected cases")
