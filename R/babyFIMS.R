library(RTMB)
library(dplyr) # sry Nathan! :D
library(tidyr)

load("data/am2022.RData")
source("R/helper.R")

head(input$obsdf, 5) # long format with all observations
# found myself wanting a unique identifier for each data set... anyone else?
# obs_type # 0=catch, 1=index, 2=agecom, 3=lencomp
# nll_type # 0=dnorm, 1=dmultinom
# fit_data # 1/0=TRUE/FALSE
# fleet    # 1=fishery, 2=survey
# obs      # transformed appropriately for nll_type (becomes keep vec)
# obserror # if nll_type obs error is an input (note this is Neff for dmultinom)

#Remove length comps for now till we implement them
input$obsdf <- input$obsdf[input$obsdf$obs_type!=3,]

# data list ----
dat <- list()
dat$obs <- input$obsdf$obs
dat$aux <- input$obsdf
dat$aux <- get_id(dat$aux)
dat$aux <- get_likelihood_index(dat$aux)
dat$year <- input$years
dat$minYear <- min(dat$year)
dat$age <-  input$ages
dat$minAge <- min(dat$age)
dat$sampleTimes <- input$srv_frac
dat$spawnTimes <- input$sp_frac
dat$waa <- input$waa
dat$mature <- input$maturity
# maybe obsType instead? not sure how this will all be indexed yet
dat$fleetTypes <- unique(input$obsdf$fleet) 
dat$srmode <- 1

# prediction data frame
dat$pred <- get_pred(dat$aux)

idx <- which(!is.na(dat$y)) 
sum(dnorm(y[idx], mean=lam[idx], sd=sdObs, log=TRUE))

# parameter ----
par <- list()
par$logsigR <- log(input$sigr)
par$logQ <- log(input$q)
# is M a constant in FIMS or by year/age?
par$logM <- matrix(log(input$natmort), nrow=length(dat$year), ncol=length(dat$age))
par$rickerpar <- if(dat$srmode==1){c(1,1)}else{numeric(0)}
par$bhpar <- if(dat$srmode==2){c(1,1)}else{numeric(0)}
par$logN <- matrix(0, nrow=length(dat$year), ncol=length(dat$age))
par$logFmort <- matrix(0, nrow=length(dat$year), ncol=1)
par$logfshslx <- log(input$fsh_slx) # need parametric selectivity
par$logsrvslx <- log(input$fsh_slx)

# todo
calc_ssb <- function(Naa, Faa, M, waa, mature, spawnTimes){
  rowSums(Naa*exp((-Faa-M)*spawnTimes)*mature*waa)
}

# model ----
f<-function(par){ # note dat isn't an argument in the fxn
  getAll(par, dat) # RTMB's attach
  obs <- OBS(obs) # access simulation, OSA residuals
  
  nobs <- length(obs) 
  nyear <- length(year)
  nage <- length(age)
  
  sigR <- exp(logsigR)
  M <- exp(logM)
  Q <- exp(logQ)
  Fmort <- exp(logFmort)
  fshslx <- exp(logfshslx)
  srvslx <- exp(logsrvslx)
  Faa <- matrix(data = 1, nrow = nyear, ncol = nage) 
  for(y in 1:nyear) Faa[y,] = Fmort[i,] * fshslx 
  
  Faa+M
  #ssb <- calc_ssb(exp(logN),Faa,natmat,waa,mature)
  
  # Z = M + Fmort
  # S = exp(-Z)
  # Naa[y,a] = exp(mean_log_rec + rec_dev[y])
  # Naa[y+1,2] = Naa[y,a-1]*S[y, a-1]
  # Naa[y+1,nage] = Naa[y,nage]*S[y,nage]
  # pred catch, sum over ages
  # logPred[i] <- logN[y,a]-log(Z)+log(1-exp(-Z))+log(Faa[y,a])
  # 
  # # survey biomass index, sum over ages
  # pred = logQ + logsrvslx + log(exp(logN) * waa) - Z * sampleTimes
  # jnll <- jnll - dnorm(obs, pred, obserror, 1)
  # # age comp
  # tmp = exp(logsrvslx + logN - Z * sampleTimes)
  # pred = tmp / sum(tmp) # this is where ageing error would be applied
  # jnll <- jnll - dmultinom(obserror*obs, pred, 1)
  
  #section to calculate population abundance
  
  
  #section to calculate predicted values
  dat_aux_pred 
  
  #section to calculate likelihoods
  for(i in unique(dat$aux$likelihood_index)){
    
  }
  
  jnll <- 0
  
  for(y in 2:nyear){
    thisSSB <- ifelse((y-minAge-1)>(-.5),ssb[y-minAge],ssb[1]) 
    if(srmode==0){ #RW
      pred <- logN[y-1,1]
    }
    if(srmode==1){ #Ricker
      pred <- rickerpar[1]+log(thisSSB)-exp(rickerpar[2])*thisSSB
    }
    if(srmode==2){ #BH
      pred <- bhpar[1]+log(thisSSB)-log(1.0+exp(bhpar[2])*thisSSB)
    }
    if(!(srmode%in%c(0,1,2))){
      stop(paste("srmode", srmode,"not implemented yet"))
    }      
    jnll <- jnll - dnorm(logN[y,1],pred,sigR,log=TRUE)
  }  
  for(y in 2:nrow){
    for(a in 2:ncol){
      pred <- logN[y-1,a-1]-F[y-1,a-1]-M[y-1,a-1]
      if(a==ncol){
        pred <- log(exp(pred)+exp(logN[y-1,a]-F[y-1,a]-M[y-1,a]))
      }
      jnll <- jnll - dnorm(logN[y,a],pred,sdS,log=TRUE)
    }
  }
  logPred <- numeric(nobs)  
  for(i in 1:nobs){
    y <- aux[i,1]-minYear+1
    f <- aux[i,2]
    a <- aux[i,3]-minAge+1
    Z <- F[y,a]+M[y,a]
    if(fleetTypes[f]==0){
      logPred[i] <- logN[y,a]-log(Z)+log(1-exp(-Z))+log(F[y,a])
    }
    if(fleetTypes[f]==2){  
      logPred[i] <- logQ[keyQ[f,a]]+logN[y,a]-Z*sampleTimes[f]
    }
    if(!(fleetTypes[f]%in%c(0,2))){  
      stop("This fleet type is has not been implemented yet")
    }
  }
  
  REPORT(logPred)
  logssb<-log(ssb)
  ADREPORT(logssb)
  jnll
}    

obj <- MakeADFun(f, par, 
                 #random=c("logN", "logF", "missing"), 
                 map=list(logsdF=as.factor(rep(0,length(par$logsdF)))), 
                 silent=FALSE)
opt <- nlminb(obj$par, obj$fn, obj$gr, control=list(eval.max=1000, iter.max=1000))
opt$objective

sdr <- sdreport(obj)
plr <- as.list(sdr,report=TRUE, "Est")
plrsd <- as.list(sdr,report=TRUE, "Std")
lines(dat$year, exp(plr$logssb), lwd=3, col="darkred")
lines(dat$year, exp(plr$logssb-2*plrsd$logssb), lwd=3, col="darkred", lty="dotted")
lines(dat$year, exp(plr$logssb+2*plrsd$logssb), lwd=3, col="darkred", lty="dotted")
