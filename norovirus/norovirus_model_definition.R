'
This script defines the model for norovirus togehter with the intervention
functions and the cost function needed to perform the fit.
'

## set up the norovirus model
model <- function(time, state, parms, signalP, signalD, ...){
  signalP <- signalP(time, ti, effP)
  signalD <- signalD(time, ti, q)
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


##define the cost function
cost <- function(p) {
  out <- ode(state0, timeline, model, p, signalP=signalP, signalD=signalD)
  modCost(model = out[ , c("time","I1")] , obs = df_noro, method = "Marq")
}


## functions for health interventions
# cleaning measures at time t=ti reduce virus appearance in water by effP
signalP <- function(t, ti, effP) {
  ifelse(t < ti, 1, 1-effP)
}


# quarantine measures, signalD resembles the isolation coefficient q
signalD <- function(t, ti, q) {
  ifelse(t < ti, 0, q)
}
