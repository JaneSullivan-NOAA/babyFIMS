# allows data inputs to the dmultinom to be proportions (this is the same
# function as TMB's dmultinom)
baby_dmultinom <- function(x, N, p, dolog = FALSE){
  xp1 = x+1
  logres = lgamma(N + 1) - sum(lgamma(x+1)) + sum(x*log(p))
  if(dolog) return(logres)
  else return(exp(logres))
}
