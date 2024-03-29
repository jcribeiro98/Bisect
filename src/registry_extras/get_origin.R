library(sets)
library(glue)
library(reticulate)
source("src/HOGen/outlier_check.R")
source("src/ODM/inference_methods.R")


get_origin <- function(type){
  #' @title Get the origin for Bisect
  #' @description Auxiliary function to get the origin for BISECT. Multiple 
  #' different origin selection methods are included. The one presented in 
  #' *Efficient Generation of Hidden Outliers for Improved Outlier Detection*
  #' corresponds to *type == 'weighted'*
  #' 
  #' Arguments:
  #' @param type: Controls the origin type. Currently, we have implemented 
  #'              -centroid: õ is the centroid of the data.
  #'              -least_outliers: Given a fitted LOF, õ is selected as the 
  #'              point that obtained the lowest outlier score.
  #'              -random: Full random selection.
  #'              -weighted: Random selection with a weighted distribution using 
  #'              LOF scores. 

  if (type == "centroid"){
    print(glue("Origin method: {type}"))
    origin = colMeans(DB[2:(ncol(DB) - 1)])
  }
  if (type == "least_outlier"){
    print(glue("Origin method: {type}"))
    model <- import("pyod.models.lof")
    lof <- model$LOF()
    lof$fit(DB[2:(ncol(DB) - 1)])
    index = which.max(lof$predict_proba(DB[2:(ncol(DB) - 1)])[,1])
    
    origin = DB[index, 2:(ncol(DB) - 1)]
    names = names(origin)
    origin = as.matrix(origin)
    dim(origin) <- NULL
    
    names(origin) = names
  }
  if (type == "random"){
    print(glue("Origin method: {type}"))
    out_DB = DB %>% filter(Out == 0)
    
    index <- sample(1:nrow(out_DB), 1)
    origin <- out_DB[index, 2:(ncol(out_DB) - 1)]
    names = names(origin)
    origin = as.matrix(origin)
    dim(origin) <- NULL
    
    names(origin) = names
  }
  if (type == "weighted"){
    if (exists("proba_vector") != TRUE){
      print(glue("Origin method: {type}"))
      print(glue("Calculating probability vector..."))
      model <- import("pyod.models.lof")
      lof <- model$LOF()
      lof$fit(DB[2:(ncol(DB) - 1)])
      proba_vector = lof$predict_proba(DB[2:(ncol(DB) - 1)])[,1]
      proba_vector <<- proba_vector/sum(proba_vector)
      #DB["Out"] <<- lof$predict(DB[2:(ncol(DB)-1)])
      print(glue("Done!"))
      }
    
    DB[ncol(DB) + 1] = proba_vector
    names(DB)[ncol(DB)] = "proba_vector"
    out_DB = DB %>% filter(Out == 0)
    index <- sample(1:nrow(out_DB), 1, prob = out_DB$proba_vector)
    origin <- out_DB[index, 2:(ncol(out_DB) - 2)]
    names = names(origin)
    origin = as.matrix(origin)
    dim(origin) <- NULL
    
    names(origin) = names
  }
  
  return(origin)
}