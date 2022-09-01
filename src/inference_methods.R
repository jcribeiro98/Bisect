library(sets)
library(glue)
library(dplyr)
source("src/utils.R")


inference <- function(x,S, method = "mahalanobis"){
  
  if(exists("ODM_env", envir = globalenv()) != T){
    for (s in set_power(as.numeric(1:(ncol(DB)-2)))){
      if (set_is_empty(s) != T){
        
        ODM_env[[glue("method{set_names(s)}")]] = 2#' Here we are going to add 
                                                   #' the call of a new function
                                                   #' that fits the desired 
                                                   #' method to the subspace.
        ODM_env <<-ODM_env
        rm(ODM_env,envir = enviroment())
      }
    }
  }
  
  #Methods: Given how each method is defined in the method environment, call them
  #to perform inference 
  
  if(method == "mahalanobis"){
    result = ODM_env[[glue("method{set_names(S)}")]]
  }
  
  return(result)
}