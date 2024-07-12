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

names(input)
input$years

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
kappa <- MGpar$INIT[rownames(MGpar)=="VonBert_K_Fem_GP_1"]
cv_amin <- MGpar$INIT[rownames(MGpar)=="CV_young_Fem_GP_1"]
cv_amax <- MGpar$INIT[rownames(MGpar)=="CV_old_Fem_GP_1"]
amin <- ss3$ctl$Growth_Age_for_L1
amax <- ss3$ctl$Growth_Age_for_L2
laa <- data.frame(age = ages,
                  laa = len_amin + (len_amax - len_amin) * (1 - exp(-kappa * (ages - amin))) / 
                    (1 - exp(-kappa * (amax - amin))))
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

# obs_df -----

obs_type nll_type fit_data fleet year age len       obs obserror
1        0        0        1     1 1977  NA  NA  9.987967     0.05
2        0        0        1     1 1978  NA  NA 10.096131     0.05
3        0        0        1     1 1979  NA  NA 10.054662     0.05
4        0        0        1     1 1980  NA  NA  9.927595     0.05
5        0        0        1     1 1981  NA  NA  9.887765     0.05
6        0        0        1     1 1982  NA  NA  9.897168     0.05


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

