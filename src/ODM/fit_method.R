library(dbscan)
library(reticulate)
source("src/registry_extras/utils.R")

get_subspaces <- function(){
  if (ncol(DB)-2 <= 11){
    print(glue("Calculating all possible subspaces..."))
    supS = set_power(as.numeric(1:(ncol(DB)-2)))
    print(glue("Done! Number of total subspaces: {length(supS)}"))
  }else{
    print(glue("Selecting random subspaces..."))
    supS = random_subspace()  
    print(glue("Done! Number of subspaces selected: {length(supS)}"))
  }
  return(supS)
}

random_subspace <- function(){
  supS = set()
  
  for (i in 1:(2^11-2)){ #With replacement (feature bagging)
    print(glue("{i/(2^11 - 2)}%"))
    d = sample(1:(ncol(DB)-3),1)
    s = as.numeric(sample(1:(ncol(DB)-3), d))
    s = as.set(s)
    supS = set_union(supS, set(s))
  }
  index = array()
  j = 1
  for (i in 1:(ncol(DB)-2)){
    for(s in supS){
      for (e in s){
        if(e == i){ #Other checking methods doesn't work
          index[j] = i
          j = j + 1
          break
          }
        }
    }
  }
  j = 1
  elements_not_included = array()
  for (i in 1:(ncol(DB)-2)){
    if(!(i %in% index)){
      elements_not_included[j] = i
      j = j + 1
    }
  }
  if(!elements_not_included %is% NA){
    supS = set_union(supS, set(as.set(as.numeric(elements_not_included))))
  }
  supS = set_union(supS, set(as.set(as.numeric(1:(ncol(DB)- 2)))))
  return(supS)
}


fit_all_methods <- function(method,...){
  #' @title Fit all the methods
  #' 
  #' @description Fits the selected ODM in all possible susbpace of the DB
  #' 
  #' Arguments:
  #' @param method: ODM
  #' @param ...: Params that are passed to the fit function.
  print(glue(
  "Fitting method: {cyan$underline(method)} to all combination of subspaces: \n"
  ))
  if (method %in% c("DeepSVDD", "fast_ABOD", "ECOD","pyod_LOF")){ #Initialize python connection 
    init_python()
  }
  ODM_env = new.env()
  supS = get_subspaces()
  ODM_env$supS = supS
  
  for (s in supS){
    if (set_is_empty(s) != T){
      
      print(glue("Fitting in the feature space : {set_names(s, sep = ' & ')}"))
      ODM_env[[glue("method{set_names(s)}")]] = fit(method, s,...)
    }
  }
  ODM_env <<-ODM_env
} 
  
fit <- function(method, S,...){
  #' @title Fit function
  #' 
  #' @description Fits the selected method in the desired subspace
  #' 
  #' Arguments
  #' @param method: ODM
  #' @param S: Subspace to which fit the ODM into. 
  #' @param ...: Params passed to the fitting methods.
  
  
  #Methods:
  
  if (method == "mahalanobis"){
    S #Needed for fixing a bug regarding R lazy evaluation in function closures
    result = function(x){distmah(S,x) > critval(S,verb=F)}
  }
  
  if (method == "LOF"){
    sS = set_subspace_grab(S)
    DB_new = DB
  
    result = function(x,...){
      DB_new[nrow(DB) + 1,sS] = x[sS]
      scores = lof(DB_new[sS],...)
      crit_val = quantile(scores, .95)
      return(scores[nrow(DB) + 1] > crit_val)
    }
  }
  if (method == "pyod_LOF"){
    model <- import("pyod.models.lof")
    lof <- model$LOF(...)
    
    sS = set_subspace_grab(S)
    lof$fit(DB[sS])
    
    result = function(x){
      prediction = lof$predict(matrix(x[sS], ncol = length(sS)))
      return(prediction == 1)
    }
  }
  if (method == "DeepSVDD"){
    model <- import("pyod.models.deep_svdd")
    dsvdd <- model$DeepSVDD(...)
    
    sS = set_subspace_grab(S)
    dsvdd$fit(DB[sS])
    result <- function(x){
      prediction = dsvdd$predict(matrix(x[sS], ncol = length(sS))) 
      return(prediction == 1 )
    }
  }
  
  if (method == "fast_ABOD"){
    abod <- import("pyod.models.abod")
    fast_abod <- abod$ABOD(...)
    
    sS = set_subspace_grab(S)
    fast_abod$fit(DB[sS])
    result <- function(x){
      prediction = fast_abod$predict(matrix(x[sS], ncol = length(sS))) 
      return(prediction == 1 )
    }
  }
  
  if (method == "ECOD"){
    ecod <- import("pyod.models.ecod")
    ECOD = ecod$ECOD(...)
    
    sS = set_subspace_grab(S)
    ECOD$fit(DB[sS])
    result <- function(x){
      prediction = ECOD$predict(matrix(x[sS], ncol = length(sS)))
      return(prediction == 1 )
    }
  }
    
  return(result)  
}
