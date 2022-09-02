library(sets)
library(glue)
library(dplyr)
source("src/registry_extras/utils.R")
source("src/ODM/fit_method.R")


inference <- function(x,S, method = "mahalanobis"){
  ##Pending description##
  
  
  if(exists("ODM_env", envir = globalenv()) != T){
    fit_all_methods(method)
  }
  
  #Methods: Given how each method is defined in the method environment, call them
  #to perform inference 
  
  if(method == "mahalanobis"){
    result = ODM_env[[glue("method{set_names(S)}")]](x)
  }
  
  return(result)
}