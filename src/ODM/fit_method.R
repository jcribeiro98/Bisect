library(dbscan)
source("src/registry_extras/utils.R")


fit_all_methods <- function(method,...){
  glue("Fitting method: {cyan$underline(method)} to all combination of subspaces: \n")

  ODM_env = new.env()
  for (s in set_power(as.numeric(1:(ncol(DB)-2)))){
    if (set_is_empty(s) != T){
      
      glue("Fitting in the feature space : {set_names(s, sep = ' & ')}")
      ODM_env[[glue("method{set_names(s)}")]] = fit(method, s,...)
    }
  }
  ODM_env <<-ODM_env
}

fit <- function(method, S,...){
  ##Pending description##
  
  
  #Methods:
  
  if (method == "mahalanobis"){
    
    S #Needed for fixing a bug regarding R lazy evaluation in function closures
    result = function(x){distmah(S,x) > critval(S,verb=F)}
  }
  
  if (method == "LOF"){
    sS = set_subspace_grab(S)
    DB_new = DB

    result = function(x){
      DB_new[nrow(DB) + 1,sS] = x
      scores = lof(DB_new[sS],...)
      crit_val = quantile(scores, .95)
      return(scores[nrow(DB) + 1] > crit_val)
    }
  }
  
  return(result)  
}
