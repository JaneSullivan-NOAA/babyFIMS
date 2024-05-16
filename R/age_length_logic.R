#length/age model process

#Generic growth data to create length at age 
max_age <- 15
Linf <- 100
K <- .3
a0 <- 0
cv <- .2
alpha <- 0.001
beta <- 3
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
AtoL <- function(a,Linf,K){
  L <- Linf*(1-exp(-K*a))+0.001
}

#Function to get expected weight at length from growth params
LtoW <- function(L,alpha,beta){
  W <- alpha*L^beta
}

logistic_mature <- function(x, inflection_point, slope){
  prop <- 1 / (1 + exp(-1 * slope * (x - inflection_point)))
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

year_index <- sort(rep(years,n_ages*n_growth_morphs))

#Assume that an equal number of lengths are tracked
#for each age. This simplifies setting things up but I 
#don't think will be strickly necessary
ages_x_lengths <- sort(rep(ages,n_growth_morphs))

age_index <- rep(ages_x_lengths,n_years)

#For this example I'm setting the lengths tracked for each 
#age based on a VB growth curve and 3 standard deviations
#either side of the mean in 1 sd steps. This should mean that
#almost all individuals (99.7%) are inside these bounds 
#assuming a normal distribution.
lengths_x_ages <- AtoL(ages_x_lengths,
                       rep(Linf,length(ages_x_lengths)),
                       rep(K,length(ages_x_lengths))
                       )*rep(length_props,length(ages))

length_index <- rep(lengths_x_ages,n_years)

#For this example I'm setting the weights as fixed weight at length
weight_x_length <- LtoW(length_index,
                       rep(alpha,length(length_index)),
                       rep(beta,length(length_index)))

weight_index <- weight_x_length

#For this example set a single natural mortality for all years, ages, and lengths
nat_mort_index <- rep(0.2,n_years*n_ages*n_growth_morphs)

#For this example set constant 50% female
prop_female_index <- rep(0.5,n_years*n_ages*n_growth_morphs)

#For this example set proportion mature to logistic maturity by age
age_50_mat <- 6
slope <- 2
mature_x_age<-logistic_mature(ages,age_50_mat,slope)
prop_Mature_index <- rep(sort(rep(mature_x_age,n_growth_morphs)),n_years)

#Data frame to hold abundance data specific to year, age, and length this could 
#probably just be a vector of abundance if the indexes are used.
#QUESTION: Should we add weight to this?
# Abundance <- data.frame(abundance=rep(0,n_years*n_ages*n_growth_morphs),
#                         year=year_index,
#                         age=age_index,
#                         length=length_index,
#                         weight=weight_index)

abundance <- rep(0,n_years*n_ages*n_growth_morphs)
unfished_abundance <- rep(0,n_years*n_ages*n_growth_morphs)


#setup transition pairs
growth_source <- NULL
growth_sink <- NULL
for(i in seq_along(years)){
  for(j in rev(seq_along(ages))){
    for(k in 1:n_growth_morphs){
      if(i==1){
        #In first year no transfer so just leave in same year 
        #an initial population function with fill these.
        growth_source<-c(growth_source,which(year_index==years[i] & 
                                             age_index==ages[j])[k])
        
        growth_sink<-c(growth_sink,which(year_index==years[i] & 
                                         age_index==ages[j])[k])
      }else if(j==1){
        #The first age class will be filled using a recruitment function.
        growth_source<-c(growth_source,which(year_index==years[i] & 
                                               age_index==ages[j])[k])
        
        growth_sink<-c(growth_sink,which(year_index==years[i] & 
                                           age_index==ages[j])[k])
      }else{
        #For all other years this will define the proportional transfer
        #between age/length bins each year. Currently simplified to 100% 
        #along a range of growth morphs. Abundance will be adjusted by 
        #natural and fishing mortality impacts.
        growth_source<-c(growth_source,which(year_index==years[i-1] & 
                                               age_index==ages[j-1])[k])
        
        growth_sink<-c(growth_sink,which(year_index==years[i] & 
                                           age_index==ages[j])[k])
      }
    }
  }
}

#Data frame to specify the transfer of abundance between year/age/length bins
#after jugling this in my head I think just specifying the row of the abundance
#data frame to move from and to is the easiest way to do this??
Growth_transition <- data.frame(proportion=rep(1,length(growth_source)),
                                growth_source=growth_source,
                                growth_sink=growth_sink)    



#Population steps

#Fill the initial abundance in year 1 values 

#Fill recuitment at age 0 values each year

#Calculate the transition of abundance accounting for catch and mortality

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

ssb <- rep(NA,n_years)
unfished_ssb <- rep(NA,n_years)

catch_year <- rep(year_index,n_fleets)
catch_age <- rep(age_index,n_fleets)
catch_length <- rep(length_index,n_fleets)

fish_mort <- rep(0,n_years*n_fleets*n_ages*n_growth_morphs)

total_mort <- rep(0,n_years*n_ages*n_growth_morphs)

catch <- rep(0,n_years*n_fleets*n_ages*n_growth_morphs)

survey_year <- rep(year_index,n_surveys)
survey_age <- rep(age_index,n_surveys)
survey_length <- rep(length_index,n_surveys)

survey_obs_frac <- rep(0,n_years*n_surveys*n_ages*n_growth_morphs)

survey_intercept <- rep(0,n_years*n_surveys*n_ages*n_growth_morphs)

for(i in seq_along(Growth_transition[,1]))
{
  sink_index <- Growth_transition$growth_sink[i]
  source_index <- Growth_transition$growth_source[i]
  trans_prop <- Growth_transition$proportion[i]
  if(year_index[sink_index]==start_year){
    abundance[sink_index] <- get_init_NAA[age_index[sink_index]]
    unfished_abundance[sink_index] <- get_init_unfished_NAA[age_index[sink_index]]
  }else if(age_index[sink_index]==0){
    if(is.na(ssb[which(years==year_index[sink_index])])){
      ssb[which(years==year_index[sink_index])] <- sum(abundance[which(year_index==year_index[sink_index])]*
                                                       weight_index[which(year_index==year_index[sink_index])]*
                                                       prop_mature_index[which(year_index==year_index[sink_index])]*
                                                       prop_female_index[which(year_index==year_index[sink_index])])
      
      unfished_ssb[which(years==year_index[sink_index])] <- sum(unfished_abundance[which(year_index==year_index[sink_index])]*
                                                                weight_index[which(year_index==year_index[sink_index])]*
                                                                prop_mature_index[which(year_index==year_index[sink_index])]*
                                                                prop_female_index[which(year_index==year_index[sink_index])])
    }
    
    abundance[sink_index] <- get_recruit[ssb,steep,r_zero,phi_zero]
    unfished_abundance[sink_index] <- get_recruit[unfished_ssb,steep,r_zero,phi_zero]
  }else if(age_index[sink_index]==max_age){
    
    total_mort[sink_index] <- sum(nat_mort_index[sink_index],fish_mort[which(catch_year==year_index[sink_index] &
                                                                                   catch_age==age_index[sink_index] &
                                                                                   catch_length==length_index[sink_index])])
    
    catch[which(catch_year==year_index[sink_index] &
                  catch_age==age_index[sink_index] &
                  catch_length==length_index[sink_index])]<-catch[which(catch_year==year_index[sink_index] &
                                                                          catch_age==age_index[sink_index] &
                                                                          catch_length==length_index[sink_index])] + (fish_mort[which(catch_year==year_index[sink_index] &
                                                                                                                                        catch_age==age_index[sink_index] &
                                                                                                                                        catch_length==length_index[sink_index])]/total_mort[sink_index])*abundance[sink_index]*(1-exp(-total_mort[sink_index]))
    
    survey_intercept[which(survey_year==year_index[sink_index] &
                  survey_age==age_index[sink_index] &
                  survey_length==length_index[sink_index])]<-survey_intercept[which(survey_year==year_index[sink_index] &
                                                                          survey_age==age_index[sink_index] &
                                                                          survey_length==length_index[sink_index])] + abundance[sink_index]*(1-exp(-survey_obs_frac[which(catch_year==year_index[sink_index] &
                                                                                                                                                                            catch_age==age_index[sink_index] &
                                                                                                                                                                            catch_length==length_index[sink_index])])) 
    
    abundance[sink_index] <- abundance[sink_index]*exp(-total_mort[sink_index])
    
    total_mort[source_index] <- sum(nat_mort_index[source_index],fish_mort[which(catch_year==year_index[source_index] &
                                                                                   catch_age==age_index[source_index] &
                                                                                   catch_length==length_index[source_index])])
    
    abundance[sink_index] <- abundance[sink_index] + trans_prop*abundance[source_index]*exp(-total_mort[source_index])

    catch[which(catch_year==year_index[sink_index] &
                  catch_age==age_index[sink_index] &
                  catch_length==length_index[sink_index])]<-catch[which(catch_year==year_index[sink_index] &
                                                                          catch_age==age_index[sink_index] &
                                                                          catch_length==length_index[sink_index])] + (fish_mort[which(catch_year==year_index[source_index] &
                                                                                                                                        catch_age==age_index[source_index] &
                                                                                                                                        catch_length==length_index[source_index])]/total_mort[source_index])*trans_prop*abundance[source_index]*(1-exp(-total_mort[source_index])) 
    
  }else{
    total_mort[source_index] <- sum(nat_mort_index[source_index],fish_mort[which(catch_year==year_index[source_index] &
                                                                     catch_age==age_index[source_index] &
                                                                     catch_length==length_index[source_index])])
    
    abundance[sink_index] <- trans_prop*abundance[source_index]*exp(-total_mort[source_index])
    
    catch[which(catch_year==year_index[sink_index] &
                  catch_age==age_index[sink_index] &
                  catch_length==length_index[sink_index])]<-catch[which(catch_year==year_index[sink_index] &
                                                                          catch_age==age_index[sink_index] &
                                                                          catch_length==length_index[sink_index])] + (fish_mort[which(catch_year==year_index[source_index] &
                                                                                                                                      catch_age==age_index[source_index] &
                                                                                                                                      catch_length==length_index[source_index])]/total_mort[source_index])*trans_prop*abundance[source_index]*(1-exp(-total_mort[source_index]))    
  }
}



#TODO: Still need to write up logic for interpolating between abundance values
#to calculate length comps.
#
for(i in fleets){
  
}

#
#


