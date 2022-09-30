library(sets)
library(glue)
source("src/ODM/inference_methods.R")
source("src/registry_extras/classes.R")
source("src/HOGen/outlier_check.R")

hidden_sample = function(gen_points = 100, 
                         eps, 
                         l = min(DB[2:(ncol(DB)-1)]), 
                         u = max(DB[2:(ncol(DB)-1)])){
  #' @title Sampling method given by Georg's paper
  #' 
  #' @describeIn Sampling method described by Georg in his paper: Hiding 
  #' outliers in high dimensional data. 
  #' 
  #' Arguments:
  #' 
  #' @param B : Number of samplings to perform (defaults at 1000).
  #' @param eps : *eps* value as given by Georg in his paper.
  #' @param l : Lower bounding of the data space (defaults to the minimum value 
  #'            of the data base)
  #' @param u : Upper bounding of the data space (defaults to the maximum value
  #'            of the data base)
  
  
  Y = DB[sample(nrow(DB), gen_points, replace = T),2:(ncol(DB)-1)]
  x = matrix(nrow = nrow(Y), ncol = ncol(Y))
  name= matrix(0,nrow = 1, ncol = ncol(x))
  for (i in 1:ncol(x)){
    name[i] = glue('y{i}')
  }
  colnames(x) = name
  
  for (j in 1:nrow(Y)){
    y = Y[j,]
    
    k = 1
    for (i in y){
      lower_Change = runif(1, min = 0, max = i - l)
      upper_Change = runif(1, min = 0, max = u - i)
      
      if (runif(1, min = 0, max = 1) < 0.5){ x[j,k] = i - eps*lower_Change}
      else{x[j,k] = i + eps*upper_Change}
      
      k = k + 1
    }
  }
  return(x)
}


main_hidden <-function(gen_points = 100, 
                       eps = eps, 
                       l = min(DB[2:(ncol(DB)-1)]), 
                       u = max(DB[2:(ncol(DB)-1)]),
                       method = "mahalanobis",
                       check_version = "fast",
                       dev_opt = FALSE){
  #' @title Main function for the hidden algorithm
  #' 
  #' @description Performs the hidden sampling methods as described in Georg's 
  #' work and then checks if the sampling points are hidden outliers or not. 
  #' The function then returns a matrix containing all the values found to be 
  #' hidden outliers. 
  #' 
  #' Arguments:
  #'         -Identical to *hidden_sample*-
  
  
  x_list = hidden_sample(gen_points = gen_points, 
                         eps = eps, 
                         l = min(DB[2:(ncol(DB)-1)]), 
                         u = max(DB[2:(ncol(DB)-1)]))
  
  hidden_x_list = matrix(0, nrow = nrow(x_list), ncol = ncol(x_list))
  hidden_x_type = matrix(0, nrow = nrow(x_list), ncol = 1)
  tic()
  for (i in 1:nrow(x_list)){
    if (check_version == "fast"){
      check_if_outlier = outlier_check_fast(x_list[i,])
    }else{
      check_if_outlier = outlier_check(x_list[i,])
    }
    if(check_if_outlier[[1]] == 0){ 
      print(glue('x in {check_if_outlier[[2]]}'))
      hidden_x_list[i,] = x_list[i,] 
      hidden_x_type[i,] = check_if_outlier[[2]]}
  }
  exec_time = toc()
  exec_time = exec_time$callback_msg
  
  name= matrix(0,nrow = 1, ncol = ncol(hidden_x_list)) 
  for (i in 1:ncol(hidden_x_list)){
    name[i] = glue('y{i}')
  }
  colnames(hidden_x_list) = name#We implement the same 
                                #naming structure just 
                                #because it's easier then 
                                #to manually debug some
                                #results.
  
  hidden_x_list = hidden_x_list[rowSums(hidden_x_list) != 0,]
  gen_result = hog_method(DB, gen_points, method, "Hidden", 
                          ODM_env, hidden_x_list, hidden_x_type, exec_time)
  
  if(dev_opt == F){
    rm(ODM_env, envir = globalenv())}
  return(gen_result)
}
