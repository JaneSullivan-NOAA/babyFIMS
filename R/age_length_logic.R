#length/age model process test

#Generic growth data to create length at age this could be replaced
#by empirical length at age in practice if available but provides 
#an easy way to parameterize a generic example population
max_age <- 15
Linf <- 100
K <- .3
a0 <- 0
cv <- .2
alpha <- 0.001
beta <- 3
rzero <- 10000
steep <- 1

#Set up length/age composition data ranges
comp_ages <- 0:max_age
comp_lengths <- 0:ceiling((1+3*cv)*Linf)

#Set up number of growth morphes to track
n_growth_morphs <- 7
length_props <- 1+(-3:3)*cv
rec_props <- dnorm(length_props,1,cv)
rec_props <- rec_props/sum(rec_props)

#Vector to specify the ages tracked in the model
#For simplicity start with assumption that ages in 
#comps and model are the same
ages <- comp_ages
n_ages <- length(ages)
#Function to get expected length at age from growth params
AtoL <- function(a_f,Linf_f,K_f,a_0){
  L <- Linf_f*(1-exp(-K_f*(a_f-a_0)))+0.001
}

#Function to get expected weight at length from growth params
LtoW <- function(L_f,alpha_f,beta_f){
  W <- alpha_f*L_f^beta_f
}

logistic_mature <- function(x_f, inflection_point_f, slope_f){
  prop <- 1 / (1 + exp(-1 * slope_f * (x_f - inflection_point_f)))
}

get_recruit <- function(ssb_f,rzero_f,steep_f,phi_0_f){
  recruits <- (0.8*rzero_f*steep_f*ssb_f)/(0.2*phi_0_f*rzero_f*(1.0-steep_f)+ssb_f*(steep_f-0.2))
}
#Vector of model years 
start_year <- 1
end_year <- 10
years <- start_year:end_year
n_years <- length(years)

#Set up indexing vectors for all tracking across years, ages, and lengths

#For simplicity start with the assumption that these are not time varying 
#Assume that an equal number of lengths are tracked
#for each age. This simplifies setting things up but I 
#don't think will be stricktly necessary

year <- sort(rep(years,n_ages*n_growth_morphs))

#Assume that an equal number of lengths are tracked
#for each age. This simplifies setting things up but I 
#don't think will be strickly necessary
ages_x_lengths <- sort(rep(ages,n_growth_morphs))

age <- rep(ages_x_lengths,n_years)

#For this example I'm setting the lengths tracked for each 
#age based on a VB growth curve and 3 standard deviations
#either side of the mean in 1 sd steps. This should mean that
#almost all individuals (99.7%) are inside these bounds 
#assuming a normal distribution.
lengths_x_ages <- AtoL(ages_x_lengths,
                       rep(Linf,length(ages_x_lengths)),
                       rep(K,length(ages_x_lengths)),
                       rep(a0,length(ages_x_lengths))
                       )*rep(length_props,length(ages))

len <- rep(lengths_x_ages,n_years)

#For this example I'm setting the weights as fixed weight at length
weight_x_length <- LtoW(len,
                       rep(alpha,length(length)),
                       rep(beta,length(length)))

weight <- weight_x_length

#For this example set a single natural mortality for all years, ages, and lengths
nat_mort <- rep(0.2,n_years*n_ages*n_growth_morphs)

#For this example set constant 50% female
prop_female <- rep(0.5,n_years*n_ages*n_growth_morphs)

#For this example set proportion mature to logistic maturity by age
age_50_mat <- 6
slope <- 2
mature_x_age<-logistic_mature(ages,age_50_mat,slope)
prop_mature <- rep(sort(rep(mature_x_age,n_growth_morphs)),n_years)

#Data frame to hold abundance data specific to year, age, and length this could 
#probably just be a vector of abundance if the indexes are used.
# Abundance <- data.frame(abundance=rep(0,n_years*n_ages*n_growth_morphs),
#                         year=year,
#                         age=age,
#                         length=length,
#                         weight=weight)

abundance <- rep(0,n_years*n_ages*n_growth_morphs)
unfished_abundance <- rep(0,n_years*n_ages*n_growth_morphs)

init_NAAL <- rep(0,n_ages*n_growth_morphs)
init_unfished_NAAL <- rep(0,n_ages*n_growth_morphs)
init_age <- ages_x_lengths
init_length <- lengths_x_ages
#Initialize ssb
ssb <- rep(NA,n_years)
unfished_ssb <- rep(NA,n_years)

#setup transition pairs
growth_source <- NULL
growth_sink <- NULL
for(i in seq_along(years)){
  for(j in rev(seq_along(ages))){
    for(k in 1:n_growth_morphs){
      if(i==1){
        #In first year no transfer so just leave in same year 
        #an initial population function with fill these.
        growth_source<-c(growth_source,which(year==years[i] & 
                                             age==ages[j])[k])
        
        growth_sink<-c(growth_sink,which(year==years[i] & 
                                         age==ages[j])[k])
      }else if(j==1){
        #The first age class will be filled using a recruitment function.
        growth_source<-c(growth_source,which(year==years[i] & 
                                               age==ages[j])[k])
        
        growth_sink<-c(growth_sink,which(year==years[i] & 
                                           age==ages[j])[k])
      }else{
        #For all other years this will define the proportional transfer
        #between age/length bins each year. Currently simplified to 100% 
        #along a range of growth morphs. Abundance will be adjusted by 
        #natural and fishing mortality impacts.
        growth_source<-c(growth_source,which(year==years[i-1] & 
                                               age==ages[j-1])[k])
        
        growth_sink<-c(growth_sink,which(year==years[i] & 
                                           age==ages[j])[k])
      }
    }
  }
}

#Specify the transfer of abundance between year/age/length bins
transition_proportion <- rep(1,length(growth_source))   

#Setup fleet and survey structure
n_fleets <- 2

if(n_fleets>0){
  fleets <- 1:n_fleets
}else{
  fleets <- NULL
}

n_surveys <- 1

if(n_surveys>0){
  surveys <- 1:n_surveys
}else{
  surveys <- NULL
}

fleet <- sort(rep(fleets,n_years*n_ages*n_growth_morphs))
fleet_year <- rep(year,n_fleets)
fleet_age <- rep(age,n_fleets)
fleet_length <- rep(len,n_fleets)

fish_mort <- rep(0,n_years*n_fleets*n_ages*n_growth_morphs)

total_mort <- rep(0,n_years*n_ages*n_growth_morphs)

catch_abun <- rep(0,n_years*n_fleets*n_ages*n_growth_morphs)

survey <- sort(rep(surveys,n_years*n_ages*n_growth_morphs))
survey_year <- rep(year,n_surveys)
survey_age <- rep(age,n_surveys)
survey_length <- rep(len,n_surveys)

survey_q <- rep(0,n_years*n_surveys*n_ages*n_growth_morphs)

survey_abun <- rep(0,n_years*n_surveys*n_ages*n_growth_morphs)

for(i in seq_along(transition_proportion))
{
  sink <- growth_sink[i]
  source <- growth_source[i]
  trans_prop <- transition_proportion[i]
  fleet_sink <- which(fleet_year==year[sink] &
                      fleet_age==age[sink] &
                      fleet_length==len[sink])
  fleet_source <- which(fleet_year==year[source] &
                        fleet_age==age[source] &
                        fleet_length==len[source])
  survey_sink <- which(survey_year==year[sink] &
                       survey_age==age[sink] &
                       survey_length==len[sink])
  survey_source <- which(survey_year==year[source] &
                         survey_age==age[source] &
                         survey_length==len[source])
  
  if(year[sink]==start_year){
    abundance[sink] <- init_NAAL[which(init_age==age[sink] & init_length==len[sink])]
    unfished_abundance[sink] <- init_unfished_NAAL[which(init_age==age[sink] & init_length==len[sink])]
  }else if(age[sink]==0){
    if(is.na(ssb[which(years==year[sink])])){
      ssb[which(years==year[sink])] <- sum(abundance[which(year==year[sink])]*
                                           weight[which(year==year[sink])]*
                                           prop_mature[which(year==year[sink])]*
                                           prop_female[which(year==year[sink])])
      
      unfished_ssb[which(years==year[sink])] <- sum(unfished_abundance[which(year==year[sink])]*
                                                    weight[which(year==year[sink])]*
                                                    prop_mature[which(year==year[sink])]*
                                                    prop_female[which(year==year[sink])])
    }
    abundance[sink] <- get_recruit(ssb[which(years==year[sink])],rzero,steep,(unfished_ssb[1]/rzero))
    unfished_abundance[sink] <- get_recruit(unfished_ssb[which(years==year[sink])],rzero,steep,unfished_ssb[1]/rzero)
  }else if(age[sink]==max_age){
    
    total_mort[sink] <- sum(nat_mort[sink],fish_mort[fleet_sink])
    
    catch_abun[fleet_sink]<-catch_abun[fleet_sink] + (fish_mort[fleet_sink]/total_mort[sink])*abundance[sink]*(1-exp(-total_mort[sink]))
    
    survey_abun[survey_sink]<-survey_abun[survey_sink] + survey_q[survey_sink]*abundance[sink]*(1-exp(-total_mort[sink])) 
    
    abundance[sink] <- abundance[sink]*exp(-total_mort[sink])
    
    total_mort[source] <- sum(nat_mort[source],fish_mort[fleet_source])
    
    abundance[sink] <- abundance[sink] + trans_prop*abundance[source]*exp(-total_mort[source])

    catch_abun[fleet_sink]<-catch_abun[fleet_sink] + (fish_mort[fleet_source]/total_mort[fleet_source])*trans_prop*abundance[source]*(1-exp(-total_mort[source])) 
    
    survey_abun[survey_sink]<-survey_abun[survey_sink] + survey_q[survey_source]*trans_prop*abundance[source]*(1-exp(-total_mort[source])) 
    
  }else{
    total_mort[source] <- sum(nat_mort[source],fish_mort[fleet_source])
    
    abundance[sink] <- trans_prop*abundance[source]*exp(-total_mort[source])
    
    catch_abun[fleet_sink]<-catch_abun[fleet_sink] + (fish_mort[fleet_sink]/total_mort[source])*trans_prop*abundance[source]*(1-exp(-total_mort[source]))    
    
    survey_abun[survey_sink]<-survey_abun[survey_sink] + survey_q[survey_source]*trans_prop*abundance[source]*(1-exp(-total_mort[source])) 
  }
}

#Here we calculate interpolated relative abundance at length 
#param len #is the length that relative abundance is being calculated for
#param year #is the time at which this is being calculated (fixed at whole years for now)
#param pop_val_id #is a vector of values identifying the row id of model estimates
#                  to include in the interpolation

comp_lengths
ages

#Need to add fleets to this. Should set up so that there can be fleet specific
#comp bins. Just start with single population comps for now though.

#First set up vectors to store indexing values of year, age, composition length,
#the model length for lower interpolation bound, model length for upper interpolation
#bound, and the abundance for the composition component.
length_comp_year <- sort(rep(years,(length(ages)*length(comp_lengths))))
length_comp_age <- rep(sort(rep(ages,length(comp_lengths))),length(years))
length_comp_length <- (rep(comp_lengths,length(ages)*length(years)))
length_comp_lower_interp <- rep(NA,length(length_comp_length))
length_comp_upper_interp <- rep(NA,length(length_comp_length))
length_comp_abund <- rep(NA,length(length_comp_length))

#Setup temporary storage vectors to calculate interpolation bounds and abundance
temp_abun_per_rec<-NULL
temp_comp_lower_interp<-rep(NA,length(comp_lengths))
temp_comp_upper_interp<-rep(NA,length(comp_lengths))
comp_lower_interp <- NULL
comp_upper_interp <- NULL

#Sequence along all the model ages to calculate the mean length at that age
#and a normal distribution of lengths around the mean as well as the closest
#model lengths at age to identify interpolation bounds.
for(i in seq_along(ages)){
  
  mean_length <- AtoL(ages[i],Linf,K,a0)
  
  temp_len_probs<-pnorm(q=comp_lengths,mean=mean_length,sd=mean_length*cv)
  temp_len_probs[1]<-0
  temp_len_probs<- c(temp_len_probs[-1],1)-temp_len_probs
  
  temp_abun_per_rec <- c(temp_abun_per_rec,temp_len_probs)
  
  pop_lengths <- mean_length*length_props
  
  temp_comp_lower_interp[comp_lengths<pop_lengths[1]] <- 0
  temp_comp_upper_interp[comp_lengths<pop_lengths[1]] <- 1
  for(j in seq_along(pop_lengths)[-1]){
    temp_comp_lower_interp[comp_lengths<pop_lengths[j] &
                             comp_lengths>=pop_lengths[j-1] ] <- j-1
    
    temp_comp_upper_interp[comp_lengths<pop_lengths[j] &
                             comp_lengths>=pop_lengths[j-1] ] <- j
  }
  temp_comp_lower_interp[comp_lengths>=pop_lengths[length(pop_lengths)]] <- length(pop_lengths)
  temp_comp_upper_interp[comp_lengths>=pop_lengths[length(pop_lengths)]] <- length(pop_lengths)+1

  comp_lower_interp <- c(comp_lower_interp,temp_comp_lower_interp)
  comp_upper_interp <- c(comp_upper_interp,temp_comp_upper_interp)
}
#replicate abundance and interpolation values across years 
length_comp_no_mort_abund_per_rec <- rep(temp_abun_per_rec,length(years))
length_comp_lower_interp <- rep(comp_lower_interp,length(years))
length_comp_upper_interp <- rep(comp_upper_interp,length(years))
length_comp_prop_upper <- (length_comp_length-length_comp_lower_interp)/
                          (length_comp_upper_interp-length_comp_lower_interp)

#Loop over all length components to calculate abundance
for(i in seq_along(length_comp_year)){
  length_comp_abund[i] <- recruitment(year,age)*
                          length_comp_no_mort_abund_per_rec[i]*(
                            
                          (1-length_comp_prop_upper[i])*abundance[year==length_comp_year[i] &
                                     age==length_comp_age[i] &
                                     len==length_comp_lower_interp[i]]
                          +
                          length_comp_prop_upper[i]*abundance[year==length_comp_year[i] &
                                       age==length_comp_age[i] &
                                       len==length_comp_upper_interp[i]]
                          )
}

#Add code to sum length components across ages 





