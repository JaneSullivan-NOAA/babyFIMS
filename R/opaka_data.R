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
names(ss3$ctl)[grepl('Gro', names(ss3$ctl))]
ss3$ctl$GrowthModel

MGpar <- ss3$ctl$MG_parms

# stock info ----
stock_info <- "opakapaka"
years <- ss3$dat$styr:ss3$dat$endyr
ages <- 0:(ss3$dat$Nages)
lens <- ss3$dat$lbin_vector
sp_frac <- 0/12 # Q: assume spawning at beg of yr
srv_frac <- 6/12 # Q: assume survey midsummer
# waa <- 999 # Q: No longer a mandatory input but there should a warning of either provide WAA or LAA/LW
natmat <- MGpar$INIT[rownames(MGpar)=="NatM_p_1_Fem_GP_1"]

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

# Schnute parameterization of the von Bertalanffy growth function
AtoL <- function(ages, len_amin, len_amax, schnute_kappa, amin, amax){
  
  laa = len_amin + (len_amax - len_amin) * (1 - exp(-schnute_kappa * (ages - amin))) / 
    (1 - exp(-schnute_kappa * (amax - amin)))
  return(laa)
  
}

# allometric function
LtoW <- function(laa, lw_alpha, lw_beta) {
  waa <- lw_alpha * laa ^ lw_beta
  return(waa)
}

mean_length_age <- AtoL(ages, len_amin, len_amax, schnute_kappa, amin, amax)
mean_weight_age <- LtoW(laa = mean_length_age, lw_alpha, lw_beta)

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

sizeage_df <- as.data.frame(sizeage) %>% 
  mutate(length = lens)
names(sizeage_df) <- c(ages, "length")

sizeage_df <- sizeage_df %>% 
  tidyr::gather("age", "obs", -c("length")) %>% 
  mutate(age = as.numeric(age))

# axis <- tickr(al_key, age, 5)
ggplot(sizeage_df, aes(x = age, y = length, size = obs)) +
  geom_point(shape = 21, fill = "black") +
  scale_size(range = c(0, 3)) +
  guides(size = FALSE) +
  labs(x = "Observed age", y = "Length (cm)", size = NULL) +
  theme_bw()

ggsave("data/opaka/sizeage.png")


# maturity -----
mat_l50 <- MGpar$INIT[rownames(MGpar)=="Mat50%_Fem_GP_1"]
matk <- MGpar$INIT[rownames(MGpar)=="Mat_slope_Fem_GP_1"]
mat <- data.frame(length_cm = lens) %>% 
  mutate(mature =  1 / (1 + exp(matk * (lens - c(mat_l50)))))
ggplot(mat, aes(length_cm, mature)) + 
  geom_point() + geom_line()
ggsave('data/opaka/maturity_at_length.png')

maturity_at_length <- mat$mature
maturity_at_age <- maturity_at_length %*% sizeage %>% as.vector()

ggplot(data.frame(age = ages, mature = maturity_at_age), 
       aes(age, mature)) + 
  geom_point() + geom_line()

# numbers-at-age vector -> numbers-at-length vector (with sizeage transition matrix)
# numbers-at-length vector -> mature numbers-at-length vector (by multiplying by maturity at length vector)
# mature numbers-at-length vector -> mature numbers-at-weight (by multiplying by
# weight at length vector, or fecundity relationship if fecundity not
# proportional to weight)

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
# 1=fishery, 2=camera on research fishing survey, 4 = research fishing survey
indices <- ss3$dat$CPUE %>% 
  dplyr::group_by(index, year) %>% 
  dplyr::summarize(obs = log(sum(obs)), 
                   obserror = mean(se_log)) %>% 
  dplyr::filter(year != -999) %>% 
  dplyr::rename(fleet = index) %>%
  dplyr::mutate(obs_type = 1, 
                nll_type = 0,
                fit_data = ifelse(fleet == 2, 1, 0), # Q: only fit to survey camera survey
                age = NA,
                len = NA) 

# length comps ----
names(ss3$dat)[grepl('lencomp', names(ss3$dat))]
lencomps <- ss3$dat$lencomp %>%
  tidyr::pivot_longer(cols = -c(Yr, Seas, FltSvy, Gender, Part, Nsamp), 
                      names_to = 'len', values_to = 'obs')  %>% 
  dplyr::mutate(len = gsub(pattern = 'l', replacement = '', x = len)) %>% 
  dplyr::select(year = Yr, fleet = FltSvy, len, obs, obserror = Nsamp) %>% 
  dplyr::mutate(fleet = ifelse(fleet < 0, -1*fleet, fleet)) %>% 
  dplyr::group_by(year, fleet) %>% 
  dplyr::mutate(obs = obs / sum(obs)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(obs_type = 3, 
                nll_type = 1,
                fit_data = ifelse(fleet == 2, 1, 0), # Q: only fit to survey camera survey
                age = NA,
                len = as.numeric(len)) 
 
tstcomps <- lencomps %>% 
  dplyr::group_by(year, fleet) %>% 
  dplyr::summarise(tst = round(sum(obs),0)) %>% 
  dplyr::pull(tst)

if (!all(tstcomps %in% 1)) {
  stop("length comps not summing to 1")
}

obsdf <- catch %>% 
  bind_rows(indices) %>% 
  bind_rows(lencomps) %>% 
  select(obs_type, nll_type, fit_data, fleet, year, age, len, obs, obserror)

ss3out <- r4ss::SS_output('data/opaka')
names(ss3out)

# selectivity ----
sel <- ss3out$ageselex
fsh_slx <- sel[sel$Factor=="Asel2" & sel$Fleet == 1,-c(1:7)] %>% colMeans() %>% unname()
srv_slx <- sel[sel$Factor=="Asel2" & sel$Fleet == 2,-c(1:7)] %>% colMeans() %>% unname()
srv_slx_fleet4 <- sel[sel$Factor=="Asel2" & sel$Fleet == 4,-c(1:7)] %>% colMeans() %>% unname()
data.frame(age = ages,
           camera_survey = srv_slx,
           fishing_survey = srv_slx_fleet4,
           fishery = fsh_slx) %>% 
  pivot_longer(cols = -age, names_to = 'fleet', values_to = 'selectivity') %>% 
  ggplot(aes(age, selectivity, col = fleet)) +
  geom_point() + geom_line() +
  scale_x_continuous(limits = c(0,12))
ggsave("data/opaka/slx.png")

sigR <- ss3out$sigma_R_in
ss3waa <- colMeans(unname(ss3out$wtatage[-c(1:6)]))
mean_weight_age
plot(x = ages, y = ss3waa)
points(x = ages, y = mean_weight_age, col = 'blue')
# not sure why there is a discrepancy but we'll use ss3's output for waa
waa <- ss3waa

# catchability ----
pars <- ss3out$parameters
rownames(pars)[grepl('Q', rownames(pars))]
names(ss3out)[grepl('batage', names(ss3out))]
# ss3out$fleet_type
q <- exp(pars$Value[rownames(pars)=="LnQ_base_BFISH(2)"])
fisherycpue_q <- exp(pars$Value[rownames(pars)=="LnQ_base_FRS(1)"])
exp(pars$Value[rownames(pars)=="LnQ_base_BFISH_ResFish(4)"])
ss3out$parameters

input <- list(stock_info = stock_info,
              obsdf = obsdf, 
              years = years, ages = ages, lens = lens,
              sp_frac = sp_frac, srv_frac = srv_frac,
              waa = waa, maturity = maturity_at_age,
              sizeage = t(sizeage), fsh_slx = fsh_slx,
              srv_slx = srv_slx, sigr = sigR, 
              natmort = natmat, q = q)
save(input, file = 'data/opaka/opaka.Rdata')
