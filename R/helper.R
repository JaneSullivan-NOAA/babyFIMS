
#Add a likelihood calculation index to the obs data frame
get_likelihood_index <- function(myaux) {
  # myaux <- input$obsdf
  myaux$likelihood_index <- NA
  like_index <- 1
  # Loop over all observations to identify multinomial
  # observations with linked likelihoods
  for(i in seq_along(myaux$obs_type)){
    if(i==1){
      #Set the first observation at index 1
      myaux$likelihood_index[i] <- like_index 
    }else if(myaux$obs_type[i]!=myaux$obs_type[i-1]){
      #If observation type changes then increment
      like_index <- like_index + 1
      myaux$likelihood_index[i] <- like_index 
    }else if(myaux$obs_type[i]<2){
      #All catches and indices are independent values so always increment
      like_index <- like_index + 1
      myaux$likelihood_index[i] <- like_index 
    }else if(myaux$year[i]!=myaux$year[i-1]){
      #Increment if the year of composition changes
      like_index <- like_index + 1
      myaux$likelihood_index[i] <- like_index 
    }else if(myaux$fleet[i]!=myaux$fleet[i-1]){
      #Increment if the fleet of composition changes
      like_index <- like_index + 1
      myaux$likelihood_index[i] <- like_index 
    }else{
      #If year and fleet don't change then comps will be combined for a 
      #single multinomial likelihood calculation
      myaux$likelihood_index[i] <- like_index 
    }
  }
  return(myaux)
}