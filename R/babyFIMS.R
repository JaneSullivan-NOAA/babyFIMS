library(RTMB)
library(dplyr) 
library(tidyr)
library(ggplot2)
library(ggthemes)
compiler::enableJIT(0)

load("data/am2022.RData")
load("data/sizeage_matrix.RData")
input$sizeage <- sizeage; rm(sizeage)
source("R/helper.R")
source("R/obj_fn.R")

head(input$obsdf, 5) # long format with all observations
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
dat$len <-  input$lens
dat$minAge <- min(dat$age)
dat$sampleTimes <- input$srv_frac
dat$spawnTimes <- input$sp_frac
dat$waa <- input$waa
dat$mature <- input$maturity
dat$sizeage <- input$sizeage
dat$fleetTypes <- unique(input$obsdf$fleet) 
dat$srmode <- 0 #
dat$logN_mode <- 0 # 0 = deterministic SCAA, 1 = sigR estimated, logR estimated as ranef, 2 = full state-space with/ shared sigN for a > 1 (same as n_NAA_sigma in WHAM)

# prediction data frame
dat$aux <- get_pred(dat$aux, input)

# parameter ----
par <- list()
par$logsigR <- log(input$sigr)
par$logsigN <- if(dat$logN_mode==2){log(0.5)}else{numeric(0)}
par$logQ <- 0
# is M a constant in FIMS or by year/age?
par$logM <- matrix(log(input$natmort), nrow=length(dat$year), ncol=length(dat$age))
par$rickerpar <- if(dat$srmode==1){c(1,1)}else{numeric(0)}
par$bhpar <- if(dat$srmode==2){c(1,1)}else{numeric(0)}
par$logN <- matrix(10, nrow=length(dat$year), ncol=length(dat$age)) # Tim suggested initializing at 10 rather than zero
par$logFmort <- matrix(0, nrow=length(dat$year), ncol=1)
par$logfshslx <- log(input$fsh_slx) # need parametric selectivity
par$logsrvslx <- log(input$srv_slx)

# assumes vectors at age are supplied to function for N, F, M, W, and mature
calc_ssb <- function(Naa, Faa, M, waa, mature, spawnTimes){
  sum(Naa*exp((-Faa-M)*spawnTimes)*mature*waa)/1e3
}

# model ----
   

fill_vals <- function(x,vals){rep(as.factor(vals), length(x))}
map <- list()
map$logsigR <- if(dat$logN_mode==0){fill_vals(par$logsigR, NA)}else{factor(1)}
# map$logQ <- fill_vals(par$logQ, NA)
map$logM <- fill_vals(par$logM, NA)
map$logfshslx <- fill_vals(par$logfshslx, NA)
map$logsrvslx <- fill_vals(par$logsrvslx, NA)

nyr <- length(dat$year)
nage <- length(dat$age)
tmp <- matrix(data = NA, ncol = nage, nrow = nyr)
tmp[,1] <- 1:nyr
tmp[1,2:nage] <- (nyr+1):(nyr+nage-1)
map$logN <- as.factor(as.vector(tmp))

obj <- MakeADFun(obj_fn, par, 
                 map=map,
                 random=NULL,
                 silent=FALSE)
opt <- nlminb(obj$par, obj$fn, obj$gr, control=list(eval.max=1000, iter.max=1000))
opt$objective
obj$report()
names(opt$par)
opt$par

sdr <- sdreport(obj)
sdr
plr <- as.list(sdr,report=TRUE, "Est")
plrsd <- as.list(sdr,report=TRUE, "Std")
load('data/orig_am2022.Rdata')

Rec <- as.data.frame(arep$R)
names(Rec) <- c('year', 'R', 'Rec_sd', 'Rec_lci', 'Rec_uci')
Rec <- Rec %>%  mutate(version = 'amak') %>% 
  bind_rows(data.frame(year = dat$year,
               R = exp(plr$predlogR)/1e6,
               Rec_uci = exp(plr$predlogR+2*plrsd$predlogR),
               Rec_lci = exp(plr$predlogR-2*plrsd$predlogR),
               version = 'babyFIMS'))
#View(Rec)
ggplot(Rec %>% filter(year >= 1978), aes(year, (R), col = version)) +
  geom_point() +
  geom_line() + ylim(0,NA) + ggthemes::theme_few()


ssb <- as.data.frame(arep$SSB)
names(ssb) <- c('year', 'ssb', 'ssb_sd', 'ssb_lci', 'ssb_uci')
ssb <- ssb %>% 
  mutate(version = 'amak') %>% 
  bind_rows(data.frame(year = dat$year,
               ssb = exp(plr$logssb),
               ssb_uci = exp(plr$logssb+2*plrsd$logssb),
               ssb_lci = exp(plr$logssb-2*plrsd$logssb),
               version = 'babyFIMS'))

ggplot(ssb %>% filter(year >= 1978), aes(year, (ssb), col = version)) +
  geom_point() +
  geom_line() + ylim(0,NA) + ggthemes::theme_few()

################
#tjm
################

obj$report() #doesn't work for me

source("R/obj_fn.R") #an edited version of f using a hand coded dmultinom, and some adjustment to the nll penalty for recruitment
obj_fn(par)

other_obj <- MakeADFun(obj_fn, par, 
                 map=map,
                 random=NULL,
                 silent=FALSE)

other_obj$report() #works
other_opt <- nlminb(other_obj$par, other_obj$fn, other_obj$gr, control=list(eval.max=1000, iter.max=1000))
other$opt$obj - opt$obj # same
other_sdr <- sdreport(other_obj)
