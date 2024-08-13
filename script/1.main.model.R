# original fitting procedure was written by Stefan Flasche. 
# this code has been modified by Deus Thindwa to run in fitmodel framework using fitR package
# optimal pneumococcal vaccine strategies to reduce the burden of IPD in HIV-infected adults
# age-structured mathematical model to capture pneumococcal transmission dynamics and impact of vaccination policies
# 01/03/2022

# ================================================================

# load the require packages
if (!require(pacman)){ #load packages
  install.packages("pacman")
  install.packages("RcmdrPlugin.KMggplot2")
  #install_github("sbfnk/fitR") # required gfortran MacOS Sierra installed prior
}

pacman::p_load(char = c("tidyverse", "remotes","dplyr", "here", "plyr","lubridate", "reshape2","deSolve", "adaptivetau","data.table",
                        "scales","readr","MASS","rootSolve","binom","coda","Rcpp", "gmm","fitR", "RcppArmadillo", "devtools", "lattice"))

# ================================================================

# create a example data frame for pneumococcal carriage
mwsp <- 
  data.frame(agegp = c("<3y", "3-5y", "6-17y", "18+y"),
             agegpupper = c(3, 6, 18, 75),
             popn = c(2094, 15239, 40324, 153082),
             vtprev = c(0.266, 0.245, 0.018, 0.008), 
             nvtprev = c(0.074, 0.202, 0.092, 0.013),
             sampsize = c(41, 63, 55, 360),
             clearVT = 1/c(72, 34.8, 18, 17),
             clearNVT = 1/c(72, 34.8, 18, 17))

# create a pneumo carriage dataset to to fit to the model
mwspfit <- 
  data.frame(
    time = c(1, 2, 3, 4),
    obsS = c(27, 35, 49, 352),
    obsV = c(11, 15, 1, 3),
    obsN = c(3, 13, 5, 5),
    sampsize = c(41, 63, 55, 360))

# ================================================================


SIS <- list(
  
  # model name
name = c("SIS with constant population size"),

# names of states
state.names = c("S", "V", "N", "B"),

# names of parameters
theta.names = c("scaleV", "susN1", "susN2", "susN3", "susN4"),

  # ---------------------------------------------------------------

# helps to run the model to generate data (simulations)
simulate =
function(theta, init.state, times) {

  SIS_ode <- function(time, state, parameters) {
    
    # constant parameters
    no.agegps = 4
    agegp.l = mwsp$agegpupper - c(0, mwsp$agegpupper[-no.agegps])
    ageing = 1/(agegp.l*365)
    age.in = 0 #agein = c(1/sum(1/ageing) , ageing[-no.agegps])
    age.out = 0 #ageing
    pop = mwsp$popn
    comp = rep(1/2, 4)
    clearV = mwsp$clearVT
    clearN = mwsp$clearNVT
    
    # contact matrix
    M = matrix(0, nrow = 4, ncol = 4)
    M[1,1] = 1.0; M[1,2] = 0.1; M[1,3] = 0.6; M[1,4] = 0.5
    M[2,1] = 0.5; M[2,2] = 0.5; M[2,3] = 0.1; M[2,4] = 0.6
    M[3,1] = 0.9; M[3,2] = 0.9; M[3,3] = 1.0; M[3,4] = 0.5
    M[4,1] = 0.5; M[4,2] = 0.1; M[4,3] = 0.1; M[4,4] = 0.7
    Mix_mat = M
    
    # age-specific susceptibility parameters
    susN = c(parameters[["susN1"]], parameters[["susN2"]], parameters[["susN3"]], parameters[["susN4"]])
    susV = c(parameters[["scaleV"]])*susN
    
    # combines age-specific susceptibility parameters and average rates of social contacts (between two individuals)
    betaV = sweep(Mix_mat, 1, susV, "*") # for transmissibility replace with 2 (for every row, operation)
    betaN = sweep(Mix_mat, 1, susN, "*") # for transmissibility replace with 2 (for every row, operation)
    
    # age-structured model states (<3, 3-5, 6-17, 18+)
    S <- state[1:no.agegps]
    V <- state[(no.agegps+1):(2*no.agegps)]
    N <- state[(2*no.agegps+1):(3*no.agegps)]
    B <- state[(3*no.agegps+1):(4*no.agegps)]
    
    # total population
    Total <- S + V + N + B
    
    # sum across the columns of contact matrix beta*I*S to compute the force of infection
    foiV = betaV %*% ((V+B)/Total) %>% as.vector()
    foiN = betaN %*% ((N+B)/Total) %>% as.vector() 
    
    # differential equations
    dS <- -foiV*S - foiN*S + clearV*V + clearN*N - age.out*S + age.in*c(sum(Total), S[-no.agegps])
    dV <- foiV*S - comp*foiN*V - clearV*V + clearN*B - age.out*V + age.in*c(0, V[-no.agegps])
    dN <- foiN*S - comp*foiV*N - clearN*N + clearV*B - age.out*N + age.in*c(0, N[-no.agegps])    
    dB <- comp*foiV*N + comp*foiN*V - clearV*B - clearN*B - age.out*B + age.in*c(0, B[-no.agegps])
    
    return(list(cbind(dS, dV, dN, dB)))
  }
  
  # contribution of dual carriage to either vaccine-type (btov) or non-vaccine-type (bton) carriage
  btov = 0.35
  bton = 0.65
  
  # simulated data
  trajectory <- data.frame(ode(y = init.state,
                               times = times,
                               func = SIS_ode,
                               parms = theta,
                               method = "lsoda")) %>%
    
   dplyr::filter(time == which.max(time))
      mutate(V1 = (V1 + (btov*B1)),
             V2 = (V2 + (btov*B2)),
             V3 = (V3 + (btov*B3)),
             V4 = (V4 + (btov*B4)),
             N1 = (N1 + (bton*B1)),
             N2 = (N2 + (bton*B2)),
             N3 = (N3 + (bton*B3)),
             N4 = (N4 + (bton*B4))) %>%
        
        dplyr::mutate(time = 1) %>%
        dplyr::select(time, V1:V4) %>%
        pivot_longer(cols = c(V1:V5), names_to = "serotype", values_to = "I", values_drop_na = TRUE) %>%
        mutate(time = seq(1L:5L), sampsize = c(41, 63, 55, 262, 98)) %>%
        dplyr::select(time, I, sampsize)
  #   
  # trajectory <- 
  #   cbind(
  #     trajectory %>% 
  #       dplyr::select(time, S1:S5) %>% 
  #       pivot_longer(cols = c(S1:S5), names_to = "serotype", values_to = "modS", values_drop_na = TRUE) %>%
  #       mutate(time = seq(1L:5L)),
  #     
  #     trajectory %>% 
  #       dplyr::select(time, V1:V5) %>%
  #       pivot_longer(cols = c(V1:V5), names_to = "serotype", values_to = "modV", values_drop_na = TRUE) %>%
  #       mutate(time = seq(1L:5L)),
  #     
  #     trajectory %>% 
  #       dplyr::select(time, N1:N5) %>%
  #       pivot_longer(cols = c(N1:N5), names_to = "serotype", values_to = "modN", values_drop_na = TRUE) %>%
  #       mutate(time = seq(1L:5L))
  #     ) %>%
  #     
  #    mutate(sampsize = c(41, 63, 55, 360)) %>%
  #    dplyr::select(time, modS, modV, modN, sampsize)


  return(trajectory)
},

# ================================================================

# calculates the prior density of each parameter
dprior =
  function(theta, log = FALSE) {
    log.prior.scaleV <- dunif(theta[["scaleV"]], min = 0, max = 2, log = TRUE)
    log.prior.susN1 <- dunif(theta[["susN1"]], min = 0, max = 2, log = TRUE)
    log.prior.susN2 <- dunif(theta[["susN2"]], min = 0, max = 2, log = TRUE)
    log.prior.susN3 <- dunif(theta[["susN3"]], min = 0, max = 2, log = TRUE)
    log.prior.susN4 <- dunif(theta[["susN4"]], min = 0, max = 2, log = TRUE)
    
    log.sum <- log.prior.scaleV + log.prior.susN1 + log.prior.susN2 + log.prior.susN3 + log.prior.susN4
    
    return(ifelse(log, log.sum, exp(log.sum)))
  },
#dprior(c(susV1=0.02, susV2=0.02, susV3=0.02, susV4=0.02, susV5=0.02, 
#         susN1=0.02, susN2=0.02, susN3=0.02, susN4=0.02, susN5=0.02), log=TRUE)

# ================================================================

# generates an observation from the model run (prevalence is observed through a binomial process)
rPointObs =
  function(model.point, theta){
    obs.point <- rbinom(n = 1, size = model.point[["sampsize"]], prob = (model.point[["I"]])/50)
    return(c(obs = obs.point))
  },
#rPointObs(model.point = c(sampsize = 20, I = 31), theta)

# ================================================================

# computes likelihood of a given observed data point w.r.t model data point (prevalence is observed through a binomial process)
dPointObs =
  function(data.point, model.point, theta, log = FALSE){
    return(dmultinom(x = c(data.point[["obsS"]], data.point[["obsV"]], data.point[["obsN"]]),
                     size = data.point[["sampsize"]], 
                     prob = c(model.point[["modS"]], model.point[["modV"]], model.point[["modN"]]), 
                     log = TRUE))
  }

#dPointObs(data.point = c(27, 11, 3),  41, model.point = c(0.006201376, 1.418082, 2.544851), log = TRUE)

)

# assign SIS model structure to a 'fitmodel' class
class(SIS) <- "fitmodel"

# ================================================================

# EXAMPLE
theta <- list(scaleV = 1,
              susN1 = 0.90,
              susN2 = 0.90,
              susN3 = 0.90,
              susN4 = 0.90)

init.state <- c(S = c(48, 50, 50, 50), 
                V = c(1, 0, 0, 0), 
                N = c(1, 0, 0, 0), 
                B = c(0, 0, 0, 0))

times <- 1:1000

# simulate the model functions above
traj <- SIS$simulate(theta, init.state, times)

# plot the simulated data from model run
plotTraj(traj)

# total model population
traj <- traj %>% left_join(traj %>% dplyr::select(S1:B4) %>% mutate(tot = rowSums(.)))

# ================================================================

# computes trajectory log-likelihood of theta by summing point likelihoods 
source(here::here("script", "2.fun.likelihood.calc.R"))

# ================================================================

# generates a single random observation from a single point in a model trajectory
source(here::here("script", "3.fun.traj.sim.R"))

# ================================================================

# evaluates posterior distribution at theta and returns a result suitable for mcmcMH
source(here::here("script", "4.fun.posterior.eval.R"))

# ================================================================

# wrapper for a fun that evaluates posterior distribution at theta and returns a result suitable for mcmcMH
source(here::here("script", "5.fun.posterior.eval.wrap.R"))

# ================================================================

# runs mcmcMH with adaptive proposal distribution
source(here::here("script", "6.fun.aMCMC.MH.R"))

# ================================================================

# runs mcmcMH with adaptive proposal distribution
source(here::here("script", "7.fun.model.fitobs.R"))

# ================================================================

# conduct trace analysis

# save the image to HHD for future use
save.image(here("data", "trace.Rdata"))

# using coda R package convert the trace to mcmc object
mcmc.trace <- mcmc(trace$trace)

# combine trace chains if more than one
mcmc.trace <- mcmc(trace$trace)
#mcmc.trace1 <- mcmc(trace1$trace)
#mcmc.trace2 <- mcmc(trace2$trace)
#mcmc.trace <- mcmc.list(list(mcmc.trace1, mcmc.trace2))

# check the acceptance rate (it should usually be between 0.1 to 0.6)
1 - rejectionRate(mcmc.trace)

# plot the trace using lattice package
xyplot(x = mcmc.trace)

# determine the suitable burning objectively using effective sample size
# A good optimal burn-in length would be when the ESS has hit its maximum for all parameters
plotESSBurn(mcmc.trace)

# create new trace without the burning 
# checks for flat bits (chain stays in the same state for too long) or too many consecutive steps in one direction by sampler
burn.trace <- burnAndThin(mcmc.trace, burn = 1000)
xyplot(x = burn.trace)

# ESS an estimate for the number of independent samples (taking into account autocorrelations) generated by the MCMC run
# check the effective sample size
effectiveSize(burn.trace)

# check the autocorrection of samples where it drops most and use that for thinning
# the lag-k autocorrelation is the correlation between every sample and the sample k steps before. 
acfplot(x = burn.trace)

# conduct thinning using appropriate value from eyeballing the autocorrelation plot
thin.burn.trace <- burnAndThin(burn.trace, thin = 30)
xyplot(x = thin.burn.trace)
acfplot(x = thin.burn.trace)
effectiveSize(thin.burn.trace)

# compare the posterior estimate for unthinned and thinned traces
# naive standard error is the standard error for the mean adjusting for the sample size
# time-series standard error corrects the “naive” standard error for autocorrelations
summary(burn.trace)
summary(thin.burn.trace)

# check the density plots for unthinned and thinned traces
plotPosteriorDensity(list(unthinned = burn.trace, thinned = thin.burn.trace))

# density plot of parameters
densityplot(thin.burn.trace)

# assess model fit to the observed data using 100 thetas sampled from posterior distribution
SIS_plots(trace = burn.thin.trace, 
                 fitmodel = SIS, 
                 init.state = c(S = c(47, 50, 50, 50, 50), 
                                V = c(1, 0, 0, 0, 0), 
                                N = c(1, 0, 0, 0, 0), 
                                B = c(1, 0, 0, 0, 0)),
                 data = mwspfit %>% dplyr::select(time, obs))
