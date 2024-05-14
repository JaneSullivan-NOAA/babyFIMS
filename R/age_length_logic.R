#length/age model process

#Generic growth data to create length at age 
max_age <- 15
Linf <- 100
K <- .3
a0 <- 0
cv <- .2

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
  L <- Linf*(1-exp(-K*a))
}

#Vector of model years 
years <- 1:10
n_years <- length(years)

#Set up indexing vectors for all tracking across years, ages, and lengths


#For simplicity start with the assumption that these are not time varying 

#Assume that an equal number of lengths are tracked
#for each age. This simplifies setting things up but I 
#don't think will be strickly necessary

year_index <- sort(rep(years,n_ages*n_growth_morphs))

#Data frame to hold the length/age combinations
#that will be tracked in the model. For simplicity
#start with the assumption that these are not time varying 

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

#Data frame to hold abunance data specific to year, age, and length this could 
#probably just be a vector of abundance if the indexes are used.
#QUESTION: Should we add weight to this?
Abundance <- data.frame(abundance=rep(0,n_years*n_ages*n_growth_morphs),
                        year=year_index,
                        age=age_index,
                        length=length_index)


#setup transition pairs
abundance_source <- NULL
abundance_sink <- NULL
for(i in seq_along(years[-1])){
  for(j in seq_along(ages[-1])){
    for(k in 1:n_growth_morphs){
      abundance_source<-c(abundance_source,which(Abundance[,"year"]==years[i] & 
                                        Abundance[,"age"]==ages[j])[k])
      
      abundance_sink<-c(abundance_sink,which(Abundance[,"year"]==years[i+1] & 
                                Abundance[,"age"]==ages[j+1])[k])
    }
  }
}

#Data frame to specify the transfer of abundance between year/age/length bins
#after jugling this in my head I think just specifying the row of the abundance
#data frame to move from and to is the easiest way to do this??
Growth_transition <- data.frame(proportion=rep(1,length(abundance_source)),
                                abundance_source=abundance_source,
                                abundance_sink=abundance_sink)    


#Data frame to hold catch data specific to year, fleet, age, and length
Catch <- data.frame(Catch=c(),year=c(),fleet=c(),age=c(),length=c())

#Population steps

#Fill the initial abundance in year 1 values 

#Fill recuitment at age 0 values each year

#Calculate the transition of abundance accounting for catch and mortality

#TODO: I think I need to modify Growth transition to include 
#the year 1 and age 0 points and have secondary functions to 
#calculate those values.
for(i in seq_along(Growth_transition[,1]))
{
  
}

#TODO: Still need to write up logic for interpolating between abundance values
#to calculate length comps.
#
#
#


