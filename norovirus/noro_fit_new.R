library(deSolve)
library(readxl)
library(tidyverse)
library(FME)
library(ggplot2)

# =============================================================================
# Load Data

#ti <- 6   # time of health interventions, outbreak11
#ti <- 4    # outbreak2
ti <- 3    # outbreak6

sheet <- "outbreak6"    # adapt to correct outbreak

# consider outbreaks with pathogen compartment
df_noro <- as.data.frame(read_excel("norovirus/relevantData.xlsx", sheet=sheet))
plot(df_noro$I1)
abline(v=ti, col="red")
text(ti, max(df_noro$I1), label ="Health Intervention", pos=4, col="red")

# Define the ODE system
model <- function(time, state, parms) {
  signalP_val <- signalP(time, parms["effP"], ti)
  signalD_val <- signalD(time, parms["q"], ti)
  with(as.list(c(parms, state)), {
    fI1 <- 1  # Example value, adjust as necessary
    fP <- 1   # Example value, adjust as necessary
    dS <- -a1 * (fI1 * I1 + fP * P) / N * S
    dE <- a1 * (fI1 * I1 + fP * P) / N * S - a2 * E
    dI1 <- a2 * E - ((1 - signalD_val) * a6 + signalD_val * a13) * I1
    dD1 <- signalD_val * a13 * I1 - a15 * D1
    dR <- (1 - signalD_val) * a6 * I1 + a15 * D1
    dP <- alpha_p * signalP_val * fI1 * I1 - mu_p * P
    list(c(dS, dE, dI1, dD1, dR, dP))
  })
}

# Define the signal function for quarantine measures
signalD <- function(t, q, ti) {
  if (t < ti) {
    return(0)
  } else {
    return(q)
  }
}

# Define the signal function for health interventions
signalP <- function(t, effP, ti) {
  if (t < ti) {
    return(1)
  } else {
    return(effP)
  }
}

# Adjusted initial parameters to create a higher peak earlier
initial_parms <- c(a1 = 18, a2 = 0.05, a6 = 0.05, a13 = 0.7, a15 = 0.2, mu_p = 0.01, q = 0.4, effP = 0.2, alpha_p = 0.1, N = 3142)

# Initial state
initial_state <- c(S = 3100, E = 30, I1 = 5, D1 = 0, R = 0, P = 0.1)

# Time points
times <- df_noro$time

# Observed data for I1
observed_I1 <- df_noro$I1

# Define cost function
cost_function <- function(parms) {
  out <- ode(y = initial_state, times = times, func = model, parms = parms)
  predictions <- out[, "I1"]
  sum((observed_I1 - predictions)^2)
}

# Optimization with wider parameter ranges and tighter bounds
fit <- optim(par = initial_parms, fn = cost_function, method = "L-BFGS-B",
             lower = c(a1 = 0.01, a2 = 0.01, a6 = 0.01, a13 = 0.01, a15 = 0.01, mu_p = 0.001, q = 0.1, effP = 0.01, alpha_p = 0.01, N = 1000),
             upper = c(a1 = 20, a2 = 10, a6 = 1, a13 = 10, a15 = 1, mu_p = 1, q = 1, effP = 1, alpha_p = 1, N = 5000))

# Extract fitted parameters
fitted_parameters <- fit$par

# Print fitted parameters
print(fitted_parameters)

# Solve ODE with fitted parameters
out_fitted <- ode(y = initial_state, times = times, func = model, parms = fitted_parameters)

# Plot results with increased y-axis range
plot(times, observed_I1, type = "b", col = "red", main = "ODE Model Fit", xlab = "Time", ylab = "I1", ylim = c(0, max(observed_I1, out_fitted[, "I1"]) + 10))
lines(times, out_fitted[, "I1"], col = "blue")
legend("topright", legend = c("Observed", "Model"), col = c("red", "blue"), lty = 1:1)

