'
This script implements a grid search to find the best fit depending on initial
parameters. Multiple fits are performed and the value of the cost function (SSR)
of the fitted parameters can be compared.
'


## Load required
library(deSolve)
library(FME)
source('norovirus_model_definition.R')
source('norovirus_plots.R')
set.seed(123)

## load data of respective outbreak (outbreak11)
df_noro <- as.data.frame(read_excel("data/relevantData.xlsx", sheet="outbreak11"))
timeline <- seq(df_noro$time[1], df_noro$time[length(df_noro$time)], 0.001)
ti <- 6


## set fixed parameters obtained from reference/ set by assumption
a2 <- 1                         # cf omega (W)
a6 <- 0.3*0.03846 + 0.7*0.3333  # recovery rate weighted (0.3*symptomatic+0.7*asymptomatic)
a13 <- 2                        # isolation rate cf gamma''
a15 <- 0.3333                   # recovery rate of detected (only symptomatic), cf gamma
mu_p <- 0.1                     # cf epsilon
a1 <- 1                         # fixed constant
q <- 1                          # rate of infected being sent to quarantine (assumption)
effP <- 0.7                     # efficiency of cleaning measures (assumption)
alpha_p <- 1                    # rate of pathogen appearance in water (assumption)


## set state0
N <- 1751 # N can be found in the Chinese version of the data set (outbreak-specific)
state0 <- c(S=N-df_noro$I1[1],
            E=df_noro$I1[1]/(1-0.3), # 0.3 is proportion asymptomatic
            I1=df_noro$I1[1], D1=0, R=0, P=0)


## set parameter ranges that are supposed to be covered by grid search
range_fI1 <- seq(0.5, 1, 0.1)
range_fP <- seq(0.5, 1, 0.1)


## create parameter grid
grid <- expand.grid(fI1=range_fI1, fP=range_fP)


## perform grid search
fits <- list()

for (i in 1:length(grid[,1])){
  parms <- c(fI1=grid[i,]$fI1, fP=grid[i,]$fP)
  fits[[i]] <- modFit(f=cost, p=parms, lower=c(0,0))
  if (i%%5==0){print(i)} # progress control for larger grids
}


## evaluate results of grid search
df_fits <- data.frame("fI1"=numeric(), "fP"=numeric(), "ssr"=numeric())
for (i in 1:length(fits)){
  df_fits[i,] <- c(fits[[i]]$par[1], fits[[i]]$par[2], fits[[i]]$ssr)
}

# look at outputs and find a subset of suitable parameters
# trade off between low ssr and reasonable parameter scale
df_filtered <- df_fits[df_fits$fP>0.1,]
df_filtered <- df_filtered[df_filtered$ssr<1002,]

# I had a closer look at indices 14,20,23 and decided afterwards
# pick fit 14
# fit 14: fI1=0.826, fP=0.169
par_index <- 14
pars_est <- coef(fits[[par_index]])


# simulation of model with fitted parameters
out <- ode(state0, timeline, model, pars_est, signalP=signalP, signalD=signalD)
plot(out)
plot_against_data(df_noro, as.data.frame(out))
plot_residuals(data.frame(t=df_noro$time, residuals=fits[[par_index]]$residuals))
#plot(density(fits[[par_index]]$residuals))
