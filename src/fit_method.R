 
source("src/utils.R")


fit_all_methods <- function(method){
  ODM_env = new.env()
  for (s in set_power(as.numeric(1:(ncol(DB)-2)))){
    if (set_is_empty(s) != T){
      ODM_env[[glue("method{set_names(s)}")]] = fit(method, s)
    }
  }
  ODM_env <<-ODM_env
}

fit <- function(method, S){
  ##Pending description##
  
  
  #Methods:
  
  if (method == "mahalanobis"){
    S #Needed for fixing a bug regarding R lazy evaluation in function closures
    result = function(x){distmah(S,x) > critval(S,verb=F)}
  }
  
  return(result)  
}