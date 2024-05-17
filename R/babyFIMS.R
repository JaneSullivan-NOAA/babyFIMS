library(RTMB)
library(dplyr) 
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
dat$srmode <- 0

# prediction data frame
dat$predObs <- get_pred(dat$aux)

# parameter ----
par <- list()
par$logsigR <- log(input$sigr)
par$logsigN <- log(0.5)
par$logQ <- log(input$q)
# is M a constant in FIMS or by year/age?
par$logM <- matrix(log(input$natmort), nrow=length(dat$year), ncol=length(dat$age))
par$rickerpar <- if(dat$srmode==1){c(1,1)}else{numeric(0)}
par$bhpar <- if(dat$srmode==2){c(1,1)}else{numeric(0)}
par$logN <- matrix(0, nrow=length(dat$year), ncol=length(dat$age))
par$logFmort <- matrix(0, nrow=length(dat$year), ncol=1)
par$logfshslx <- log(input$fsh_slx) # need parametric selectivity
par$logsrvslx <- log(input$fsh_slx)

# assumes vectors at age are supplied to function for N, F, M, W, and mature
calc_ssb <- function(Naa, Faa, M, waa, mature, spawnTimes){
  sum(Naa*exp((-Faa-M)*spawnTimes)*mature*waa)
}

# model ----
f<-function(par){ # note dat isn't an argument in the fxn
  getAll(par, dat) # RTMB's attach
  obs <- OBS(obs) # access simulation, OSA residuals
  pred_Obs_vec <- rep(0, length(predObs[,1]))
  nobs <- length(obs) 
  nyear <- length(year)
  nage <- length(age)
  
  sigR <- exp(logsigR)
  sigN <- exp(logsigN)
  M <- exp(logM)
  Q <- exp(logQ)
  Fmort <- exp(logFmort)
  fshslx <- exp(logfshslx)
  srvslx <- exp(logsrvslx)
  Faa <- matrix(data = 1, nrow = nyear, ncol = nage) 
  for(y in 1:nyear) Faa[y,] = Fmort[y,] * fshslx 
  
  Z <- Faa+M
  
  jnll <- 0
  
  ssb <- rep(NA, nyear)
  for(y in 1:nyear) ssb[y] <- calc_ssb(exp(logN[y,]),Faa[y,],M[y,],waa,mature,spawnTimes)
  
  predlogR <- rep(0, nyear)
  for(y in 1:nyear){
    thisSSB <- ifelse((y-minAge-1)>(-.5),ssb[y-minAge],ssb[1]) 
    if(srmode==0){ #RW
      if (y == 1){
        predlogR[y] <- logN[y,1] # need to fix this later
      }else{
        predlogR[y] <- logN[y-1,1]
      }
    }
    if(srmode==1){ #Ricker
      predlogR[y] <- rickerpar[1]+log(thisSSB)-exp(rickerpar[2])*thisSSB
    }
    if(srmode==2){ #BH
      predlogR[y] <- bhpar[1]+log(thisSSB)-log(1.0+exp(bhpar[2])*thisSSB)
    }
    if(!(srmode%in%c(0,1,2))){
      stop(paste("srmode", srmode,"not implemented yet"))
    }      
    jnll <- jnll - dnorm(logN[y,1],predlogR[y],sigR,log=TRUE)
  }  
  
  predlogN <- matrix(0, nrow=nyear, ncol=nage)
  for(y in 2:nyear){
    for(a in 2:nage){
      predlogN[y,a] <- logN[y-1,a-1]-Faa[y-1,a-1]-M[y-1,a-1]
      if(a==nage){
        predlogN[y,a] <- log(exp(predlogN[y,a])+exp(logN[y-1,a]-Faa[y-1,a]-M[y-1,a]))
      }
      jnll <- jnll - dnorm(logN[y,a],predlogN[y,a],sigN,log=TRUE)
    }
  }

  # predicted catch ----
  
  # get_pred_logCaa(idx, # lkup vector linking to predObs, aux with appropriate flt opts
  #              Z, logN, Faa) # Z, logN, Faa are also vectors
  # need to modify this to allow for multiple fleets
  logpredcatchatage <- logN-log(Z)+log(1-exp(-Z))+log(Faa)

  # obs_type 0 is aggregate catch in weight (need to figure out units)
  for (i in which(predObs$obs_type == 0)){
    y <- which(year == predObs$year[i])
    predObs$pred[i] <- log(sum(exp(logpredcatchatage[y,]) * waa)/1e6) # waa in g and aggregate catch in t
  }
  
  # predicted survey biomass ----
  logpredindexatage <- logQ + logsrvslx + logN - Z * sampleTimes

  # obs_type 1 is survey biomass in weight (need to figure out units)
  for (i in which(predObs$obs_type == 1)){
    y <- which(year == predObs$year[i])
    predObs$pred[i] <- log(sum(exp(logpredindexatage[y,]) * waa)/1e6) # waa in g and aggregate srv biom in t
  }
  
  # age comps (age error not included) ----
  
  # fishery
  tmp <- exp(logpredcatchatage)
  tmptot <- rowSums(tmp)
  tmp <- tmp/tmptot
  # subsets out 1989 and 2022, because they are not in predObs
  tmp <- tmp[which(year %in% predObs$year[which(predObs$obs_type == 2 & predObs$fleet == 1)]),]
  
  # survey
  tmp2 <- exp(logpredindexatage)
  tmptot2 <- rowSums(tmp2)
  tmp2 <- tmp2/tmptot2
  tmp2 <- tmp2[which(year %in% predObs$year[which(predObs$obs_type == 2 & predObs$fleet == 2)]),]
  
  # combine and vectorize
  tmp3 <- rbind(tmp,tmp2)
  out <- tmp3[1,]
  
  # wow!
  for(i in seq_along(tmp3[,1])[-1]) out <- c(out, tmp3[i,])
  predObs$pred[which(predObs$obs_type == 2)] <- out 
  
  # observation likelihoods
  
  for (i in unique(predObs$likelihood_index[!is.na(predObs$likelihood_index)])){
    
    tmp <- predObs[which(predObs$likelihood_index==i), ] # drop retains matrix structure
    
    unique_nll_type <- unique(tmp$nll_type)
    if(length(unique_nll_type)>1) stop("multiple nll types within tmp")
    
    # dnorm for catches, indices
    if(unique_nll_type==0) {
      browser()
      jnll <- jnll - RTMB::dnorm(tmp$obs, tmp$pred, tmp$obserror, 1)
    }
    # multinomial for comps
    if(unique_nll_type==1) {
      # for()
      # sum(tmp$pred)
      # sum(tmp$obs + 1e-10)
      jnll <- jnll - RTMB::dmultinom(tmp$obserror * tmp$obs, NULL, tmp$pred, 1)
    }
  }
  

  # jnll <- jnll - dnorm(obs, pred, obserror, 1)
  # # age comp
  # tmp = exp(logsrvslx + logN - Z * sampleTimes)
  # pred = ageerror * tmp / sum(tmp) # this is where ageing error would be applied
  # jnll <- jnll - dmultinom(obserror*obs, pred, 1)
  
  #section to calculate population abundance
  
  
  # #section to calculate predicted values
  # dat_aux_pred 
  # 
  # #section to calculate likelihoods
  # for(i in unique(dat$aux$likelihood_index)){
  #   
  # }
  # 
  # 
  # for(y in 2:nrow){
  #   for(a in 2:ncol){
  #     pred <- logN[y-1,a-1]-F[y-1,a-1]-M[y-1,a-1]
  #     if(a==ncol){
  #       pred <- log(exp(pred)+exp(logN[y-1,a]-F[y-1,a]-M[y-1,a]))
  #     }
  #     jnll <- jnll - dnorm(logN[y,a],pred,sdS,log=TRUE)
  #   }
  # }
  # logPred <- numeric(nobs)  
  # for(i in 1:nobs){
  #   y <- aux[i,1]-minYear+1
  #   f <- aux[i,2]
  #   a <- aux[i,3]-minAge+1
  #   Z <- F[y,a]+M[y,a]
  #   if(fleetTypes[f]==0){
  #     logPred[i] <- logN[y,a]-log(Z)+log(1-exp(-Z))+log(F[y,a])
  #   }
  #   if(fleetTypes[f]==2){  
  #     logPred[i] <- logQ[keyQ[f,a]]+logN[y,a]-Z*sampleTimes[f]
  #   }
  #   if(!(fleetTypes[f]%in%c(0,2))){  
  #     stop("This fleet type is has not been implemented yet")
  #   }
  # }
  # 
  # REPORT(logPred)
  # logssb<-log(ssb)
  # ADREPORT(logssb)
  jnll
}    

obj <- MakeADFun(f, par, 
                 #random=c("logN", "logF", "missing"), 
                 #map=list(logsdF=as.factor(rep(0,length(par$logsdF)))), 
                 silent=FALSE)
opt <- nlminb(obj$par, obj$fn, obj$gr, control=list(eval.max=1000, iter.max=1000))
opt$objective

sdr <- sdreport(obj)
plr <- as.list(sdr,report=TRUE, "Est")
plrsd <- as.list(sdr,report=TRUE, "Std")
lines(dat$year, exp(plr$logssb), lwd=3, col="darkred")
lines(dat$year, exp(plr$logssb-2*plrsd$logssb), lwd=3, col="darkred", lty="dotted")
lines(dat$year, exp(plr$logssb+2*plrsd$logssb), lwd=3, col="darkred", lty="dotted")
