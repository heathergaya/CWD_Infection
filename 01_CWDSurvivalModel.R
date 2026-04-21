################# Start Header ####
## Code for: Host-pathogen dynamics and the hidden hazards of chronic wasting disease
## Heather E. Gaya, Marcelo H. Jorge, Michael J. Chamberlain, Amy V. Nalls, 
## Nathaniel D. Denkers, Candace K. Mathiason, Mark G. Ruder, 
## Gino J. D’Angelo, and Richard B. Chandler
## Code written by Heather Gaya 
## Last modified April 21, 2026
################# End Header ####

#### Load libraries ####
library(coda)
library(nimble)
library(dplyr)
library(MCMCvis)
library(parallel)

##### *Model ####
surv_vonBert <- nimbleCode({
  for(i in 1:n.ind){ #for all individuals
    infect[i,first_live[i]] ~ dbern(p[site[i]]) #was I infected the day of capture?
    pre_daysinfect[i] ~ dunif(0, 1095) #at least 3 years
    daysinfected[i, first_live[i]] <- infect[i,first_live[i]]*pre_daysinfect[i] #status on first day of monitoring 
    for(t in (first_live[i]+1):last_live[i]){ #from 2nd to last day of monitoring
      hazard[i, t] <- p2*(1-exp(-p1[Sex[i]+1]*Prion_dens[pix[i,t-1],yr[i,t-1]]))/365 #daily hazard of CWD
      p.eff[i, t] <- (infect[i, t-1] == 1)*1 + (infect[i, t-1] == 0)*hazard[i,t] #once infected, stay infected
      infect[i, t] ~ dbern(p.eff[i,t]) #am I infected? 
      daysinfected[i, t] <- (infect[i, t] + daysinfected[i, t-1])*infect[i, t] #count days infected 
      
      logit(me.phi[i,t-1]) <- B0 + B1[infect[i,t-1]+1]*TAvg[i,t-1] + B2*Age[i,t-1] + B3*Sex[i] + B4*Rou[i,t-1] + B5*field[i,t-1]+ B6*daysinfected[i, t-1]/100 #survival
      #make sure days infected stays a relatively small number or B6 will mix poorly (hence the /100)
      survive[i,t] ~ dbern(me.phi[i,t-1]) #did I survive?
    } #end t 
  } #end i  
  
  ## Priors
  B0 ~ dnorm(2, sd= 1) #intercept
  B1[1] ~ dnorm(0, 1) #Temperature CWD negative
  B1[2] ~ dnorm(0, 1) #Temperature CWD positive
  B2 ~ dnorm(0, 1) #Age
  B3 ~ dnorm(0, 1) #Sex
  B4 ~ dnorm(0, 1) #Roughness
  B5 ~ dnorm(0, 1) #Distance to Field
  B6 ~ dnorm(0, 1) #CWD status
  
  p1[1] ~ dnorm(0,1) #effect of prion load females
  p1b ~ dnorm(0, 1) #effect of being male
  p1[2] <- p1[1]+ p1b #effect of prion load males
  p2 ~ dbeta(1,1) #max infection prob
  for(q in 1:3){
    p[q] ~ dbeta(1,1) #cwd status the day before initial capture (~site level prevalence)
  }
})

##### *Run the model ####
## grab the 'stuff' object, which includes all information needed for nimble to run the model
nimstuff <- readRDS('TVCdat4_Nov14.rds')
nimdat <- nimstuff$dat
nimconsts <- nimstuff$consts
nim.inits <- nimstuff$inits
params <- nimstuff$params


cl <- makeCluster(3)
clusterExport(cl = cl, varlist = c("nimdat", "nimconsts", "nim.inits", "params", "surv_vonBert"))
system.time(
  nim.out <- clusterEvalQ(cl = cl,{
    library(nimble) #reload packages inside clusterEvalQ
    library(coda)
    prepnim <- nimbleModel(code = surv_vonBert, constants = nimconsts,
                           data = nimdat, inits = nim.inits, calculate = T)
    prepnim$calculate()
    mcmcnim <- configureMCMC(prepnim, monitors = params, print = T)
    nimMCMC <- buildMCMC(mcmcnim) #actually build the code for those samplers
    Cmodel <- compileNimble(prepnim) #compiling the model itself in C++;
    Compnim <- compileNimble(nimMCMC, project = prepnim) # compile the samplers next
    Compnim$run(niter = 10000, nburnin = 8000, thin = 1)
    return(as.mcmc(as.matrix(Compnim$mvSamples)))
  })
)
#take about 32 minutes

saveRDS(nim.out, 'UpdatedSurvival_wSex_Jan25.rds')
MCMCtrace(nim.out, pdf = F, Rhat = T, n.eff = T)               
MCMCsummary(nim.out)
