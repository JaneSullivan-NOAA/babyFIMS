// Schnute parameterization of the von Bertalanffy growth model

// Inputs: vectors of specimen-level observed lengths and ages
// Outputs: length-at-age predictions and a size-age transition matrix

// Reference:
// Schnute and Fournier 1980: Growth structure in length-frequency analysis


#include <TMB.hpp>
#include <numeric>

#define see(object) std::cout << #object ":\n" << object << "\n";

template <class Type>
Type square(Type x){return x * x;}

template<class Type>
Type objective_function<Type>::operator() ()
{

  DATA_VECTOR(obs_age);  // observations: ages
  DATA_VECTOR(obs_len);  // observations: length
  DATA_SCALAR(amin);     // reference age 1
  DATA_SCALAR(amax);     // reference age 2
  DATA_VECTOR(age);      // user-defined vector of ages for predictions (e.g., for plotting fitted curve)
  DATA_VECTOR(mod_abin); // model age bins (for size-age transition matrix)
  DATA_VECTOR(mod_lbin); // model length bins (for size-age transition matrix)

  PARAMETER(log_len_amin);   // log-space length at reference age 1
  PARAMETER(log_len_amax);   // log-space length at reference age 2
  PARAMETER(log_kappa);      // growth rate
  PARAMETER(log_sigma_amin); // SD of length at reference age 1
  // PARAMETER(log_sigma_amax); // SD of length at reference age 2
  PARAMETER(log_sigma_slope); // log-space slope between sigma_amin and sigma_amax

  Type len_amin = exp(log_len_amin);
  Type len_amax = exp(log_len_amax);
  Type kappa = exp(log_kappa);
  Type sigma_amin = exp(log_sigma_amin);
  // Type sigma_amax = exp(log_sigma_amax);
  Type sigma_slope = exp(log_sigma_slope);

  int nobs = obs_age.size();
  int nage = age.size();
  int nabin = mod_abin.size();
  int nlbin = mod_lbin.size();
  
  // predicted values and residuals
  vector<Type> pred_len(nobs);
  vector<Type> mean_log_pred_len(nobs);
  vector<Type> residuals(nobs);
  vector<Type> pred_sigma(nobs);
  // for plotting
  vector<Type> laa(nage);
  // size-age transition matrix
  matrix<Type> sizeage_mat(nabin, nlbin);
  vector<Type> sizeage_vec(nabin); 
  vector<Type> sizeage_sigma(nabin); 
  pred_len.setZero();
  mean_log_pred_len.setZero();
  residuals.setZero();
  pred_sigma.setZero();
  laa.setZero();
  sizeage_mat.setZero();
  sizeage_vec.setZero();
  sizeage_sigma.setZero();
  
  Type obj_fun = 0;
  
  // assume CV of growth increases linearly with age
  // Type sigma_slope = (sigma_amax - sigma_amin) / (amax - amin);
  Type sigma_amax = sigma_slope * (amax - amin) + sigma_amin;
  
  // iterate through observed values
  for(int i = 0; i < nobs; i++) {
    // expected values
    pred_len(i) = len_amin + (len_amax - len_amin) * (1 - exp(-kappa * (obs_age(i) - amin))) / (1 - exp(-kappa * (amax - amin)));
    
    // age-specific sigmas
    if(obs_age(i) <= amin) { // SD amin
      pred_sigma(i) = sigma_amin;
    } else {
      if(obs_age(i) >= amax) { // SD amax
        pred_sigma(i) = sigma_amax; 
      } else { // linear interpolation
        pred_sigma(i) = sigma_amin + sigma_slope * obs_age(i);
      }
    }
    // negative log-likelihood
    mean_log_pred_len(i) = log(pred_len(i)) - Type(0.5) * pred_sigma(i) * pred_sigma(i);
    obj_fun -= dnorm(log(obs_len(i)), mean_log_pred_len(i), pred_sigma(i), 1);
    
    // residuals
    residuals(i) = (obs_len(i) - pred_len(i)) / pred_sigma(i);
  }
  
  // predicted length-at-age for plotting length-at-age
  for(int i = 0; i < nage; i++) {
    laa(i) = len_amin + (len_amax - len_amin) * (1 - exp(-kappa * (age(i) - amin))) / (1 - exp(-kappa * (amax - amin)));
  }
  
  // repeat steps to construct size-age transition matrix
  // temporary variables 
  Type Lminp = min(mod_lbin); // where lengths are uniform length bins
  Type Lmaxp = max(mod_lbin);
  Type lbin_width = mod_lbin(1) - mod_lbin(0); // length bin interval (e.g. 2.0)
  Type Ll1p = 0.0;
  Type Llp = 0.0;
  Type fac1 = 0.0; 
  Type fac2 = 0.0;
  
  for(int i = 0; i < nabin; i++) {
    
    // predicted length-at-age
    sizeage_vec(i) = len_amin + (len_amax - len_amin) * (1 - exp(-kappa * (mod_abin(i) - amin))) / (1 - exp(-kappa * (amax - amin)));
    
    // age-specific sigmas
    if(mod_abin(i) <= amin) { // SD amin
      sizeage_sigma(i) = sigma_amin;
    } else {
      if(mod_abin(i) >= amax) { // SD amax
        sizeage_sigma(i) = sigma_amax; 
      } else { // linear interpolation
        sizeage_sigma(i) = sigma_amin + sigma_slope * mod_abin(i);
      }
    }
    
    // convert log se to arithmetic scale sd
    sizeage_sigma(i) = sqrt(exp(sizeage_sigma(i)*sizeage_sigma(i))-1)*sizeage_vec(i);
    
    // construct size-age matrix
    for(int j = 0; j < nlbin; j++) {
      // upper limit smallest len bin, important colsums = 0
      if(j == 0) {
        fac1 =  (Lminp + lbin_width - sizeage_vec(i))/sizeage_sigma(i);
        sizeage_mat(i,j) = pnorm(fac1);
      } else {
        if(j == (nlbin-1)) {
          fac1 = (Lmaxp - sizeage_vec(i))/sizeage_sigma(i);
          sizeage_mat(i,j) = Type(1.0) - pnorm(fac1);
        } else {
          Ll1p = mod_lbin(j+1);
          Llp = mod_lbin(j);
          fac1 = (Ll1p - sizeage_vec(i))/sizeage_sigma(i);
          fac2 = (Llp - sizeage_vec(i))/sizeage_sigma(i);
          sizeage_mat(i,j) = pnorm(fac1) - pnorm(fac2);
        }
      }
    }
  }
  
  ADREPORT(laa);
  ADREPORT(len_amin);
  ADREPORT(len_amax);
  ADREPORT(kappa);
  ADREPORT(sigma_amin);
  ADREPORT(sigma_amax);
  
  REPORT(residuals); 
  REPORT(pred_len);
  REPORT(pred_sigma);
  REPORT(sigma_slope);
  REPORT(sizeage_mat);
  REPORT(sizeage_vec);
  REPORT(sizeage_sigma);
  
  return obj_fun;

}
