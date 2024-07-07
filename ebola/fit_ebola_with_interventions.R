# =============================================================================
# TO DO: Fitting Ebola 2014-16 Liberia
# Author: Valentin Rauscher
# =============================================================================

# Load required packages
library(tidyverse)
library(deSolve)
library(FME)



df_ebola <- read.csv(file = "./ebola/ebola_data/country_timeseries.csv")
df_ebola <- select(df_ebola, c("Date", "Day", "Cases_Guinea", "Cases_Liberia", "Cases_SierraLeone"))

plot(df_ebola$Day, df_ebola$Cases_Liberia)
relevant_data_liberia <- select(df_ebola, c("Day", "Cases_Liberia"))
relevant_data_liberia <- rename(relevant_data_liberia, c(time = Day, C = Cases_Liberia))

na_indices <- which(is.na(relevant_data_liberia$C))
relevant_data_liberia <- relevant_data_liberia[-na_indices,]
# =============================================================================
### Fitting



### definition of parameter shift
par_shift <- function(t,par_old,par_new, trans_speed){
  if (t < 196){
    return(par_old)
  } else {
    return(par_new + (par_old - par_new) * exp(- trans_speed * (t - 196)))
  }
}






#model definition
model <- function(t, x, parameters, ...){   
    S <- x[1]
    E <- x[2]
    I1 <- x[3]
    I2 <- x[4]
    R <- x[5]
    D1 <- x[6]  
    D2 <- x[7]
    C <- x[8]

    f14_3 <- parameters$f14_3
    f15_6 <- parameters$f15_6
    f17_4 <- parameters$f17_4
    fI1 <- parameters$fI1
    fI2 <- parameters$fI2
    fD <- parameters$fD    
    #delta <- parameters$delta
    #mu <- parameters$mu
    a1 <- parameters$a1
    a2 <- parameters$a2
    a3 <- parameters$a3
    a4 <- parameters$a4
    a6 <- parameters$a6
    a13 <- parameters$a13
    a14 <- f14_3*a3
    a15 <- f15_6*a6
    a17 <- f17_4*a4    
    a4_intervent_speed <- parameters$a4_intervent_speed
    fI1_intervent_speed <- parameters$fI1_intervent_speed
    fD_intervent_speed  <- parameters$fD_intervent_speed
    a13_new <- parameters$a4_new
    #parameter shifts for intervention
    a4_shift = par_shift(t, a4, 0.2040, a4_intervent_speed)
    fI1_shift = par_shift(t, fI1, 0.072, fI1_intervent_speed)
    fD_shift = par_shift(t, fD, 0.1881, fD_intervent_speed)
    fI2_shift = par_shift(t, fI2, 0.1593, log(2) / 1.115)
    a13_shift = par_shift(t, a13, a13_new, log(2) / 27.1)


    #total human population N
    N <- S+E+I1+I2+R+D1+D2
    #infective Population infH weighted by relative infectivities fE, FI1, fI2, fD
    infH <- fI1_shift * I1 + fI2_shift * I2 + fD_shift * (fI1_shift * D1 + fI2 * D2) 
    phiH <- max(0, infH / N, na.rm = T) 
    #model code
    dSdt  <- - a1 * phiH * S 
    dEdt  <- a1 * S * phiH  - a2 * E 
    dI1dt <- a2 * E - (a3 + a6 + a13_shift) * I1 
    dI2dt <- a3 * I1 - a4_shift * I2 
    dRdt  <- a4_shift * I2 + a6 * I1 + a15 * D1 + a17 * D2 
    dD1dt <- a13_shift * I1 - (a14 + a15) * D1 
    dD2dt <- a14 * D1 - a17 * D2
    dCdt  <- a2 * E
    #return list of diffs
    list(c(dSdt, dEdt, dI1dt, dI2dt, dRdt, dD1dt, dD2dt, dCdt))
  }
# Note: without vector compartments the transition rate a1 implicitly includes the vector population density



## Choose initial parameter values, time interval, initial values for ode
##from Dènes et al: delta = 1694.57, mu = 0.000077, a2 = 0.124245 (1/(1.498*7), a3= ebola_induced death rate = 0.797/7 = )

fictitous <- list(delta=1694.57, mu=0.000077 ,a1=1, a2=1/(1.498*7), a3=0.0702,
                         a4=0.1852, a6=0.0260, a13=0.01,
                         f14_3=1, f15_6=1, f17_4=1, fI1=0.076, fI2=0.303, fD=1, 
                         a4_intervent_speed=log(2) / 21.2,fI1_intervent_speed=log(2) / 12.5, fD_intervent_speed=log(2), a4_new=0.2439)
ivFictitous <- c(S = 4000000, E = 0, I1 = 8, I2 = 0, R = 0, D1 = 0, D2 = 0, C = 0)
times = rev(df_ebola$Day)

## Simulate and plot the model and the data



simulation <- ode(func = model,
      y = ivFictitous,
      times = rev(relevant_data_liberia$time),
      parms = fictitous, method="lsoda")
plot(simulation[, c(1)], simulation[, c(9)])
lines(relevant_data_liberia$time, relevant_data_liberia$C, col="blue")

# Define a cost function
cost <- function(p) {
  fI2 <- p[1]
  a4_new <- p[2]
  #pp <- list(delta=1694.57, mu=0.000077 ,a1=1, a2=1/(1.498*7), a3=0.0702, a4=0.2040, a6=0.0260, a13=0.1, f14_3=1.0, f15_6=1.0, f17_4=1.0, fI1=0.072, fI2=0.3037566, fD=0.1881,a4_intervent_speed=a4_intervent_speed,fI1_intervent_speed=fI1_intervent_speed, fD_intervent_speed= log(2))
  pp <- list(delta=1694.57, mu=0.000077 ,a1=1, a2=1/(1.498*7), a3=0.0702, a4=0.1852, a6=0.0260, a13=0.1, f14_3=1.0, f15_6=1.0, f17_4=10.0, fI1=0.076, fI2=fI2, fD=0.6165, a4_intervent_speed=log(2) / 21.2,fI1_intervent_speed=log(2) / 12.5, fD_intervent_speed=log(2), a4_new=a4_new)
  out <- ode(func = model, y = ivFictitous, times =  rev(relevant_data_liberia$time), parms = pp, method = "lsoda")
  modCost(model = out[, c(1, 9)], obs = relevant_data_liberia, weight = "none",method = "Marq")
}

starterparm <- c(0.303,0.241)
# Set boundary constraints (all positive) for parameters and fit the model
fit <- modFit(f = cost, p = starterparm, lower = c(0.1, 0.141))
summary(fit)
coef(fit)

## Plot fitted model and data
pp_beforeintervention <- list(delta=1694.57, mu=0.000077 ,a1=1, a2=1/(1.498*7), a3=0.0702, a4=0.1852, a6=0.0260, a13=0.1, f14_3=1.0, f15_6=1.0, f17_4=10.0, fI1=0.076, fI2=coef(fit)[1], fD=0.6165, a4_intervent_speed=log(2) / 21.2,fI1_intervent_speed=log(2) / 12.5, fD_intervent_speed=log(2), a4_new = coef(fit)[2])
#pp_afterintervention <- list(delta=1694.57, mu=0.000077 ,a1=1, a2=1/(1.498*7), a3=0.0702, a4=0.2040, a6=0.0260, a13=0.1, f14_3=10.0, f15_6=1.0, f17_4=1.0, fI1=0.072, fI2=coef(fit)[1], fD=0.1881)



fitted_model <- ode( ivFictitous,  rev(relevant_data_liberia$time), model, pp_beforeintervention, method = "lsoda")
# Plot for report

out_df <- data.frame(fitted_model)
out_df$C <- rev(out_df$C)
relevant_data_liberia %>%
  ggplot(aes(x = time, y = C)) +
  geom_point(color="lightseagreen", size=2) +
  labs(x="Days since Start of Large Epidemic Wave",
       y="Cumulative Number of Infected Individuals") +
  #geom_vline(xintercept = delay, color = "orange", linetype = "dashed", linewidth = 1) +
  #annotate("text", x = delay, y = -5, label = TeX("$\\tau$"), color = "orange", size = 6) +  
  #coord_cartesian(ylim = c(0,70), clip = "off") +
  theme_gray(base_size=15) +
  geom_line(aes(x=time, y=out_df$C), size=1) -> plot3

ggsave(plot = plot3, "fitted_curve.pdf", width = 6, height = 5)
