library(sets)
library(glue)
library(dplyr)
source("src/registry_extras/utils.R")
source("src/ODM/fit_method.R")


inference <- function(x,S, method = "mahalanobis",...){
  #' @title Inference function
  #' 
  #' @description Manages the inference for each model.
  #' Each model is fitted within an enviroment saved in the Glob env.
  #' This function will manage how each method is going to be called 
  #' in the enviroment to perform inference. If no method is fitted, 
  #' the function will call the fit_all_methods function to fit all 
  #' the method in all of the possible subspaces. 
  #' 
  #' Arguments:
  #' @param x: Value to infer
  #' @param S: Subspace to where perform the inferece
  #' @param method: ODM
  #' @param ...: Params passed to some ODMs.
  
  
  if(exists("ODM_env", envir = globalenv()) != T){
    fit_all_methods(method,...)
  }
  
  #Methods:Given how each method is defined in the method environment, call them
  #to perform inference 
  
  if(method %in% c("mahalanobis", "DeepSVDD", "fast_ABOD", "ECOD")){
    result = ODM_env[[glue("method{set_names(S)}")]](x)
  }
  
  if(method == "LOF"){
    result = ODM_env[[glue("method{set_names(S)}")]](x,...)
  }
  
  return(result)
}