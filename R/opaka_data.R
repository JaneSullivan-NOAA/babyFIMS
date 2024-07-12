# prepare data for opaka
library(r4ss)
library(dplyr)
library(tidyr)
library(ggplot2)

dat <- r4ss::SS_readdat(file = 'data/opaka/data_updated.ss')
ss3 <- r4ss::SS_read(dir = 'data/opaka')
names(ss3)
ss3$dat$Comments
names(ss3$ctl)
names(ss3$dat)
names(ss3$ctl)[grepl('Grow|grow|LAA', names(ss3$ctl))]
ss3$ctl$GrowthModel

MGpar <- ss3$ctl$MG_parms

# stock info ----
stock_info <- "opakapaka"
years <- ss3$dat$styr:ss3$dat$endyr
ages <- 0:(ss3$dat$Nages-1)
lens <- ss3$dat$lbin_vector
sp_frac <- 0/12 # Q: assume spawning at beg of yr
srv_frac <- 6/12 # Q: assume survey midsummer
# waa <- 999 # Q: No longer a mandatory input but there should a warning of either provide WAA or LAA/LW
natmat <- MGpar$INIT[rownames(MGpar)=="NatM_p_1_Fem_GP_1"]

# maturity -----
mat_l50 <- MGpar$INIT[rownames(MGpar)=="Mat50%_Fem_GP_1"]
matk <- MGpar$INIT[rownames(MGpar)=="Mat_slope_Fem_GP_1"]
mat <- data.frame(length_cm = lens) %>% 
  mutate(mature =  1 / (1 + exp(matk * (lens - c(mat_l50)))))
ggplot(mat, aes(length_cm, mature)) + 
  geom_point() + geom_line()
ggsave('data/opaka/maturity.png')

# length-at-age -----
len_amin <- MGpar$INIT[rownames(MGpar)=="L_at_Amin_Fem_GP_1"]
len_amax <- MGpar$INIT[rownames(MGpar)=="L_at_Amax_Fem_GP_1"]
schnute_kappa <- MGpar$INIT[rownames(MGpar)=="VonBert_K_Fem_GP_1"]
cv <- MGpar$INIT[rownames(MGpar)=="CV_young_Fem_GP_1"]
# cv_amax <- MGpar$INIT[rownames(MGpar)=="CV_old_Fem_GP_1"]
amin <- ss3$ctl$Growth_Age_for_L1
amax <- ss3$ctl$Growth_Age_for_L2
laa <- data.frame(age = ages,
                  laa = len_amin + (len_amax - len_amin) * (1 - exp(-schnute_kappa * (ages - amin))) / 
                    (1 - exp(-schnute_kappa * (amax - amin))))
ggplot(laa, aes(age, laa)) + 
  geom_point() + geom_line()
ggsave('data/opaka/lengthatage.png')

# length-weight -----
lw_alpha <- MGpar$INIT[rownames(MGpar)=="Wtlen_1_Fem_GP_1"]
lw_beta <- MGpar$INIT[rownames(MGpar)=="Wtlen_2_Fem_GP_1"]
lw <- data.frame(length_cm = lens,
                  weight_kg = lw_alpha * lens ^ lw_beta)
ggplot(lw, aes(length_cm, weight_kg)) + 
  geom_point() + geom_line()
ggsave('data/opaka/lengthweight.png')


AtoL <- function(ages, len_amin, len_amax, schnute_kappa, amin, amax){
  
  laa = len_amin + (len_amax - len_amin) * (1 - exp(-schnute_kappa * (ages - amin))) / 
    (1 - exp(-schnute_kappa * (amax - amin)))
  
}

# amax <- 43
# Linf <- 67.5
# schnute_kappa <- 0.242
# a0 <- -.384677
# amin <- 6
# cv <- cv_amax
# ages <- 0:amax
# len_bins <- 0:ceiling((1+3*cv)*len_amax)

mean_length_age <- AtoL(ages, len_amin, len_amax, schnute_kappa, amin, amax)

sizeage <- matrix(NA, nrow=length(lens), ncol=length(ages))

for(i in seq_along(ages)){
  
  #Calculate mean length at age to spread lengths around
  mean_length <- AtoL(ages[i], len_amin, len_amax, schnute_kappa, amin, amax)
  
  #Calculate the cumulative proportion shorter than each composition length
  temp_len_probs <- pnorm(q=lens, mean=mean_length, sd=mean_length*cv)
  
  #Reset the first length proportion to zero so the first bin includes all
  #density smaller than that bin
  temp_len_probs[1] <- 0
  
  #subtract the offset length probabilities to calculate the proportion in each
  #bin. For each length bin the proportion is how many fish are larger than this
  #length but shorter than the next bin length.
  temp_len_probs <- c(temp_len_probs[-1],1)-temp_len_probs
  sizeage[,i] <- temp_len_probs
}
image(t(sizeage))
colSums(sizeage) # should sum to 1 within each age class


# catches ----

# FLAG: assuming single fishery!
catch <- ss3$dat$catch %>% 
  dplyr::group_by(year) %>% 
  dplyr::summarize(obs = log(sum(catch)), 
                   obserror = mean(catch_se)) %>% 
  dplyr::filter(year != -999) %>% 
  dplyr::mutate(obs_type = 0, 
                nll_type = 0,
                fit_data = 1,
                fleet = 1,
                age = NA,
                len = NA)

catchyears <- min(catch$year):max(catch$year)
if (!all(catchyears %in% catch$year)) {
  stop("missing years in landings")
}

# indices of abundance
ss3$dat$CPUE %>% 
  dplyr::group_by(fleet, year) %>% 
  dplyr::summarize(obs = log(sum(catch)), 
                   obserror = mean(catch_se)) %>% 
  dplyr::filter(year != -999) %>% 
  dplyr::mutate(obs_type = 0, 
                nll_type = 0,
                fit_data = 1,
                fleet = 1,
                age = NA,
                len = NA)

