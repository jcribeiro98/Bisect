library(sets)
library(crayon)
library(glue)
library(car)
library(dplyr)
library(uniformly)
library(tictoc)
source('src/registry_extras/utils.R')
source('src/ODM/inference_methods.R')
source('src/registry_extras/classes.R')


bisect <- function(L, x, iternum=1000, verb = T){
  #' @title Bisection method implementation for function f
  #' 
  #' @description Performs the bisection algorithm over the function 
  #'                        f*: [0,L]-->R^2--->R
  #'                             c      cx   f(cx),
  #'              given the function f and direction x.   
  #' 
  #' Arguments:
  #' @param f : (function) Function to which apply the algorithm 
  #'            (defaults to a function in the global environment named f). 
  #'            For clarity we will use f(cx) instead of f*(c).
  #'            
  #' @param L : (real) Maximum value of the interval where f* is defined. 
  #' @param x : (array) Vector in R^2 that defines a direction to which perform
  #'            the bisection algorithm.
  #' @param iternum : (integer) Number of maximal iterations for the bisection 
  #'                  algorithm.
  #' @param verb : (logical) Controls if the function has to be verbose or not.
  
  
  a = 0; b = L
  for (i in 1:iternum){
    c = (b+a)/2
    
    check_if_outlier = f(c*x + colMeans(DB[2:(ncol(DB) - 1)]))
    outlier_indicator = check_if_outlier[[1]]
    outlier_type = check_if_outlier[[2]]
    
    if (outlier_indicator > 0){b = c} #Direction vector 
                                                   #(x*c) + an origin (colMeans) 
    else if (outlier_indicator < 0){a = c}
    else{ print(glue("x in {outlier_type}")); break}
  }
  return(list(c, outlier_type))
}


f <- function(x,verb = F, h_index = F, method = "mahalanobis"){
  #' @title Bisection main function 
  #' 
  #' @description Defines the function that is going to be used in the bisection
  #' rule during the main execution of the algorithm. The function itself is 
  #' really simple, it checks if a point x is an outlier in all of the subspaces
  #' of the total space. f is defined as
  #' 
  #'        1, if x is outside of bounds,
  #' f(x)= -1, if x is in the join acceptance area,
  #'        0, if x is in H1 U H2.
  #' 
  #' @note The function only works currently with mahanalobis distances as an 
  #' outlier detection method (will probably remain like that in this file).
  #'         
  #' Arguments:
  #' @param x :(array) Point to check
  #' @param verb : (Logical) Value controlling if f should be verbose or not (
  #'               defaults to FALSE)
  #' @param h_index: (logical) Controls if the function should keep track if
  #'                 the outlier was in H1 or H2 by adding an extra value to
  #'                 the result ('h1' or 'h2'). *WARNING*: This will change the
  #'                 result's class from an integer to a mixed array.
  
  
  h2 = matrix(0, ncol = 2^(ncol(DB)-2)-1, nrow = 1)
  h2[2^(ncol(DB)-2)-1] = 1
  
  supS = set_power(as.numeric(1:(ncol(DB)-2)))
  index = matrix(0,ncol = 2^(ncol(DB)-2)-1, nrow=1)
  
  j = 0
  for (S in supS){
    if(set_is_empty(S) != T){if (inference(x, S, method)){ index[j] = 1 }}
    j = j + 1}
  
  if(sum(index[1:(2^(ncol(DB)-2)-2)]) > 0 && index[2^(ncol(DB)-2)-1] == 0){
    result = list(0,"H1")}
  else if(isTRUE(all.equal(index, h2))){ 
    result = list(0, 'H2')}
  else if(sum(index[1:(2^(ncol(DB)-2)-2)]) > 0 && index[2^(ncol(DB)-2)-1] == 1){
    if(verb==T){
      print(glue('x Outside of bounds'))}
    result = list(1, 'OB')}
  else{
    if(verb==T){print(glue('x in the total acceptance area'))} 
    result = list(-1, 'IL')}
  
  return(result)
}


main <- function(B=100, method = "mahalanobis", seed = F,
                 verb = T, dev_opt = F){
  #' @title Main function for the bisect experiment
  #' 
  #' @description Given the number of data that you want to generate (in this 
  #' case, B=100 is default), the function tries to find a hidden outlier in the 
  #' direction of those B generated data, by the bisection rule defined with 
  #' function f.
  #' 
  #' Arguments:
  #' @param B : (integer) number of data points whose directions to check                 
  #' @param seed : (logical/integer) If seed=F, the entry will be ignored, 
  #'              if seed=T, a random seed will be generated, and then 
  #'              saved in a file within the set working directory.
  #'              If class(seed) == numerical, then it will execute the given
  #'              seed.
  #' @param verb : (logical) Control if the function have to be verbose or not.
  #' @param dev_opt : (logical) Control some developer/debug options
  
  
  if(class(seed) == 'numeric'){
    set.seed(seed)
  }
  if(seed == T){
    seed = sample(1:10^9,1)
    print("A seed has been generated")
    dir.create("seed", showWarnings = F)
    write.table(seed, file = "seed/seed.txt", sep = " ")
    set.seed(seed)
  }
  if(dev_opt == T){warning(glue("The developer option has been activated. 
                                This means that the ODM environment will not 
                                be deleted, and therefore some collusion might
                                occur with following experiments. Poceed with
                                caution
                                "))}
  
  x_list = runif_on_sphere(n = B, d = ncol(DB) - 2, r = 1) #sample the 
                                                           #directional vectors
  L = max(sqrt(rowSums(DB[2:(ncol(DB)-1)]^2)))
  hidden_x_list = matrix(0, nrow = nrow(x_list), ncol = ncol(x_list))
  hidden_x_type = matrix(0, nrow = nrow(x_list), ncol = 1)
  
  tic()
  for (i in 1:nrow(x_list)){
    bisection_results = bisect(L=L, x = x_list[i,], verb = verb)
    hidden_c = bisection_results[[1]]
    outlier_type = bisection_results[[2]]
    
    if(outlier_type %in% c("H1", "H2")){ 
      hidden_x_list[i,] = hidden_c*(x_list[i,]) + colMeans(DB[2:(ncol(DB)-1)])
      hidden_x_type[i,] = outlier_type
      }
  }
  exec_time = toc()
  exec_time = exec_time$callback_msg
  
  hidden_x_list = hidden_x_list[rowSums(hidden_x_list) != 0,]
  gen_result = hog_method(DB, B, method, "bisection", 
                          ODM_env, hidden_x_list, hidden_x_type, exec_time)
  
  if(dev_opt == F){
  rm(ODM_env, envir = globalenv())}
  return(gen_result)
}

