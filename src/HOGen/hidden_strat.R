library(sets)
library(glue)
source("src/ODM/inference_methods.R")
source("src/registry_extras/classes.R")
source("src/HOGen/outlier_check.R")

hidden_parall_routine <- function(x_list,check_version){
  hidden_results <- foreach (i = 1:nrow(x_list), .combine = rbind) %dopar% {
    if (check_version == "fast"){
      check_if_outlier = outlier_check_fast(x_list[i,])
    }else{
      check_if_outlier = outlier_check(x_list[i,])
    }
    if(check_if_outlier[[1]] == 0){ 
      print(glue('x in {check_if_outlier[[2]]}'))
      result_point = x_list[i,] 
      result_point[length(result_point) + 1] = check_if_outlier[[2]]
    }else{
      result_point = matrix(0,1,ncol(x_list) + 1)}
    result_point  
  }
  return(hidden_results)
} 

hidden_sc_routine <- function(x_list, check_version){
  hidden_x_list = matrix(0, nrow = nrow(x_list), ncol = ncol(x_list))
  hidden_x_type = matrix(0, nrow = nrow(x_list), ncol = 1) 
  for (i in 1:nrow(x_list)){
  if (check_version == "fast"){
    check_if_outlier = outlier_check_fast(x_list[i,])
  }else{
    check_if_outlier = outlier_check(x_list[i,])
    }
  if(check_if_outlier[[1]] == 0){ 
    print(glue('x in {check_if_outlier[[2]]}'))
    hidden_x_list[i,] = x_list[i,] 
    hidden_x_type[i,] = check_if_outlier[[2]]
    }
  }
  return(cbind(hidden_x_list,hidden_x_type))
}


hidden_sample = function(gen_points = 100, 
                         eps, 
                         l = min(DB[2:(ncol(DB)-1)]), 
                         u = max(DB[2:(ncol(DB)-1)])){
  #' @title Sampling method given by Steinbuss,G. and Boehm,K
  #' 
  #' @describeIn Sampling method described by Steinbuss, G. & Boehm, K. in: 
  #' Hiding outliers in high dimensional data. 
  #' 
  #' Arguments:
  #' 
  #' @param gen_points : Number of samplings to perform (defaults at 1000).
  #' @param eps : *eps* value.
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
                       num_workers = detectCores()/2,
                       until_gen = FALSE,
                       dev_opt = FALSE,...){
  #' @title Main function for the HIDDEN algorithm
  #' 
  #' @description Performs the hidden sampling methods as described in HIDDEN's
  #' paper and then checks if the sampling points are hidden outliers or not. 
  #' The function then returns a matrix containing all the values found to be 
  #' hidden outliers. 
  #' 
  #' Arguments:
  #' 
  #' @param gen_points: Number of points to check
  #' @param eps: Epsilon value to employ 
  #' @param l: Left part of each hypercube's side
  #' @param u: Right part of each hypercube's side
  #' @param method: Adversary to employ for Hidden
  #' @param check_version: *Do not touch* Parameter that controls the subspace 
  #' checking version (old one proposed by Steinbuss and Boehm vs faster 
  #' implementation) 
  #' @param num_workers: Number of workers for parallelization 
  #' @param until_gen: Boolean controlling if HIDDEN should keep generating until
  #' it finds "gen_points" hidden outliers.
  #' @param dev_opt: Activate developers options. Whenever activated all the 
  #' fitted ODM will be stored in the global environment (inside an extra 
  #' environment called *ODM_env*). Any hog_model class object will utilize these
  #' very same fitted detectors (i.e., if left in memory, after HIDDEN generates
  #' hidden outliers using the fitted detectors in*ODM_env*, BISECT will use
  #' the same detectors if *ODM_env* is still in the global environment).
  #' @param ...: Passes extra parameters to the Adversary selected. 
  
  if(exists("ODM_env", envir = globalenv()) != T){
    fit_all_methods(method,...)
  }
  tic()
  if (!until_gen){
  x_list = hidden_sample(gen_points = gen_points, 
                         eps = eps, 
                         l = l, 
                         u = u)
  
  registerDoParallel(num_workers)
  
  
  hidden_results <- foreach (i = 1:nrow(x_list), .combine = rbind) %dopar% {
    if (check_version == "fast"){
      check_if_outlier = outlier_check_fast(x_list[i,])
    }else{
      check_if_outlier = outlier_check(x_list[i,])
    }
    if(check_if_outlier[[1]] == 0){ 
      print(glue('x in {check_if_outlier[[2]]}'))
      result_point = x_list[i,] 
      result_point[length(result_point) + 1] = check_if_outlier[[2]]
    }else{
        result_point = matrix(0,1,ncol(x_list) + 1)}
      result_point  
  }
  exec_time = toc()
  exec_time = exec_time$callback_msg
  stopImplicitCluster()
  }else{
    registerDoParallel(num_workers)
    hidden_count = 0
    tic()
    while(hidden_count <= gen_points){
      print(glue("Finding Outliers... {hidden_count/gen_points * 100}%"))
      x_list = hidden_sample(gen_points = 4*num_workers, 
                             eps = eps, 
                             l = l, 
                             u = u)
    if (num_workers == 1){
      dummy_hidden_results <- hidden_sc_routine(x_list,check_version)
    }else{
      dummy_hidden_results <- hidden_parall_routine(x_list,check_version)
    }
    if(hidden_count == 0){
      hidden_results <- dummy_hidden_results
    }else{  
    hidden_results  <- rbind(dummy_hidden_results, hidden_results)}
    hidden_count <- nrow(matrix(as.numeric(hidden_results[,1:ncol(x_list)]), 
                                nrow(hidden_results))[
                                  rowSums(matrix(as.numeric(hidden_results[
                                    ,1:ncol(x_list)]), nrow(hidden_results))) 
                                  != 0, ,drop = FALSE])
    
    }
    stopImplicitCluster()
    exec_time = toc()
    exec_time = exec_time$callback_msg
    print(glue("Done!"))
  }
  
  hidden_x_list <- matrix(as.numeric(hidden_results[,1:ncol(x_list)]), 
                          nrow(hidden_results)) #' We added a char into a 
                                                   #' numeric array, so 
                                                   #' everything is a char now
  hidden_x_type <- as.matrix(hidden_results[,(ncol(x_list)+1)])
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





