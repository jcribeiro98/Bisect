library(sets)
library(glue)
library(dplyr)
source("src/utils.R")
source("src/fit_method.R")


inference <- function(x,S, method = "mahalanobis"){
  ##Pending description##
  
  
  if(exists("ODM_env", envir = globalenv()) != T){
    for (s in set_power(as.numeric(1:(ncol(DB)-2)))){
      if (set_is_empty(s) != T){
        ODM_env[[glue("method{set_names(s)}")]] = fit(method, S)
        ODM_env <<-ODM_env
        rm(ODM_env,envir = enviroment())
      }
    }
  }
  
  #Methods: Given how each method is defined in the method environment, call them
  #to perform inference 
  
  if(method == "mahalanobis"){
    result = ODM_env[[glue("method{set_names(S)}")]](x)
  }
  
  return(result)
}