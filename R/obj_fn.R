baby_dmultinom <- function(x, N, p, log = FALSE){
  xp1 = x+1
  logres = lgamma(N + 1) - sum(lgamma(x+1)) + sum(x*log(p))
  if(log) return(logres)
  else return(exp(logres))
}

obj_fn<-function(par){ # note dat isn't an argument in the fxn
  getAll(par, dat) # RTMB's attach
  obs <- OBS(obs) # access simulation, OSA residuals
  predObs <- rep(0, nrow(aux))
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

  # Recruitment ----
  ssb <- rep(0, nyear)
  for(y in 1:nyear) ssb[y] <- calc_ssb(exp(logN[y,]),Faa[y,],M[y,],waa,mature,spawnTimes)
  predlogR <- rep(0, nyear)
  for(y in 1:nyear){
    thisSSB <- ifelse((y-minAge-1)>(-.5),ssb[y-minAge],ssb[1]) 
    if(srmode==0){ # RW
      if (y == 1){
        predlogR[y] <- logN[y,1] # need to fix this later
      }else{
        predlogR[y] <- logN[y-1,1]
      }
    }
    if(srmode==1){ # Ricker
      predlogR[y] <- rickerpar[1]+log(thisSSB)-exp(rickerpar[2])*thisSSB
    }
    if(srmode==2){ # BH
      predlogR[y] <- bhpar[1]+log(thisSSB)-log(1.0+exp(bhpar[2])*thisSSB)
    }
    if(!(srmode %in% c(0, 1, 2))){
      stop(paste("srmode", srmode, "not implemented yet"))
    }      
    jnll <- jnll - dnorm(logN[y,1],predlogR[y],sigR,log=TRUE)
  }  
  
  # N matrix
  predlogN <- matrix(0, nrow=nyear, ncol=nage)
  for(y in 2:nyear){
    for(a in 2:nage){
      # full state-space model 
      predlogN[y,a] <- logN[y-1,a-1]-Faa[y-1,a-1]-M[y-1,a-1]
      if(a==nage){
        predlogN[y,a] <- log(exp(predlogN[y,a])+exp(logN[y-1,a]-Faa[y-1,a]-M[y-1,a]))
      }
      if(logN_mode == 0){ # fixed or random effects recruitment with deterministic N matrix
        logN[y,a] <- predlogN[y,a]
      }else{
        jnll <- jnll - dnorm(logN[y,a],predlogN[y,a],sigN,log=TRUE)
      }
    }
  }

  # predicted catch ----
  
  # get_pred_logCaa(idx, # lkup vector linking to aux, aux with appropriate flt opts
  #              Z, logN, Faa) # Z, logN, Faa are also vectors
  # need to modify this to allow for multiple fleets
  logpredcatchatage <- logN-log(Z)+log(1-exp(-Z))+log(Faa)

  # obs_type 0 is aggregate catch in weight (need to figure out how we want to input units)
  for (i in which(aux$obs_type == 0)){
    y <- which(year == aux$year[i])
    predObs[i] <- log(sum(exp(logpredcatchatage[y,]) * waa)/1e6) # waa in g and aggregate catch in t
  }
  
  # predicted survey biomass ----
  logpredindexatage <- logQ + logsrvslx + logN - Z * sampleTimes

  # obs_type 1 is survey biomass in weight (need to figure out units)
  for (i in which(aux$obs_type == 1)){
    y <- which(year == aux$year[i])
    predObs[i] <- log(sum(exp(logpredindexatage[y,]) * waa)/1e6) # waa in g and aggregate srv biom in t
  }
  
  # age comps (age error not included) ----
  
  # fishery
  tmp <- exp(logpredcatchatage)
  tmptot <- rowSums(tmp)
  tmp <- tmp/tmptot

  # survey
  tmp2 <- exp(logpredindexatage)
  tmptot2 <- rowSums(tmp2)
  tmp2 <- tmp2/tmptot2

  # combine and vectorize
  tmp3 <- rbind(tmp,tmp2)
  out <- tmp3[1,]
  
  # wow!
  for(i in seq_along(tmp3[,1])[-1]) out <- c(out, tmp3[i,])
  predObs[which(aux$obs_type == 2)] <- out 
  
  # length comps ----
  
  tmp5 <- exp(logpredcatchatage)
  tmptot5 <- rowSums(tmp5)
  tmp5 <- tmp5/tmptot5

  predcatchatlength <- tmp5[,rep(1,length(sizeage[1,]))]
  
  # This is a matrix version of the elementwise loop below
  #
  # predcatchatlength2 <- tmp5[,rep(1,length(sizeage[1,]))]
  # for(i in seq_along(year)){
  #   predcatchatlength2[i,]<- tmp5[i,] %*% sizeage 
  # }
  
  for(i in seq_along(year)){
    for(j in seq_along(len)){
      predcatchatlength[i,j] <- sum(sizeage[,j]*tmp5[i,])
    }
  }
  
  tmp6 <- exp(logpredcatchatage)
  tmptot6 <- rowSums(tmp6)
  tmp6 <- tmp6/tmptot6

  predcatchatlength2 <- tmp6[,rep(1,length(sizeage[1,]))]
  
  for(i in seq_along(year)){
    for(j in seq_along(len)){
      predcatchatlength2[i,j] <- sum(sizeage[,j]*tmp6[i,])
    }
  }
  
  # combine and vectorize
  tmp7 <- rbind(predcatchatlength,predcatchatlength2)
  out2 <- tmp7[1,]
  
  # wow!
  for(i in seq_along(tmp7[,1])[-1]) out2 <- c(out2, tmp7[i,])
  predObs[which(aux$obs_type == 3)] <- out2 
  
  # observational likelihoods ----
  
  for (i in unique(aux$likelihood_index[!is.na(aux$likelihood_index)])){
    
    tmp <- aux[which(aux$likelihood_index==i), ] 
    tmppred <- predObs[which(aux$likelihood_index==i)]
    
    unique_nll_type <- unique(tmp$nll_type)
    if(length(unique_nll_type)>1) stop("multiple nll types within tmp")
    
    # dnorm for catches, indices
    if(unique_nll_type==0) {
      #browser()
      jnll <- jnll - RTMB::dnorm(tmp$obs, tmppred, tmp$obserror, log=TRUE)
    }
    # multinomial for comps
    if(unique_nll_type==1) {
      # jnll <- jnll - RTMB::dmultinom(x=tmp$obserror * tmp$obs, size=sum(tmp$obserror * tmp$obs), prob=tmppred, log=TRUE)
      jnll <- jnll - baby_dmultinom(tmp$obserror * tmp$obs, sum(tmp$obserror * tmp$obs), tmppred, log=TRUE)
    }
  }

  REPORT(predObs)
  
  ADREPORT(predlogR)
  logssb<-log(ssb)
  ADREPORT(logssb)
  jnll
}    
