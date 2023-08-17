library(sets)
library(glue)
library(uniformly)
library(tictoc)
library(foreach)
library(doParallel)
source("src/registry_extras/utils.R")
source("src/ODM/inference_methods.R")
source("src/registry_extras/classes.R")
source("src/HOGen/outlier_check.R")


bisect <- function(l, x, iternum = 100, method, 
                   verb = TRUE, check_version = "fast",
                   l_val_option = "fixed",...) {
  #' @title Bisection method implementation for function f (DEPRICATED)
  #'
  #' @note THIS FUNCTION PERFORMS AN OLDER VERSION OF BISECT WITHOUT THE 
  #' CUT TRICK. 
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

  if (l_val_option != "fixed"){
    l = l + runif(1, min = -l, max = l)
  }
  a <- 0
  b <- l
  for (i in 1:iternum) {
    c <- (b + a) / 2
    
    if (check_version == "fast"){
      check_if_outlier <- outlier_check_fast(c * x + 
                                               colMeans(DB[2:(ncol(DB) - 1)]),
                                        method = method, verb = verb)
    }else {
      check_if_outlier <- outlier_check(c * x + 
                                               colMeans(DB[2:(ncol(DB) - 1)]),
                                        method = method, verb = verb)
    }
    outlier_indicator <- check_if_outlier[[1]]
    outlier_type <- check_if_outlier[[2]]

    if (outlier_indicator > 0) {
      b <- c
    } # Direction vector
    # (x*c) + an origin (colMeans)
    else if (outlier_indicator < 0) {
      a <- c
    } else {
      print(glue("x in {outlier_type} in {i} iterations"))
      break
    }
  }
  return(list(c, outlier_type))
}

grab_bisect_results <- function(i){
  #' @note THIS FUNCTION PERFORMS AN OLDER VERSION OF BISECT WITHOUT THE 
  #' CUT TRICK. 
  bisection_results <- bisect( l = l, x = x_list[i, ], method = method, 
                               verb = verb, check_version = check_version
  )
  hidden_c <- bisection_results[[1]]
  outlier_type <- bisection_results[[2]]
  
  if (outlier_type %in% c("H1", "H2")) {
    hidden_x_list[i, ] <<- hidden_c * (x_list[i, ]) +
      colMeans(DB[2:(ncol(DB) - 1)])
    hidden_x_type[i, ] <<- outlier_type
  }
}

main <- function(gen_points = 100, method = "mahalanobis", seed = FALSE,
                 verb = FALSE, dev_opt = FALSE, check_version = "fast",
                 num_workers = detectCores()/2, l_val_option = "fixed", ...) {
  #' @title Main function for the bisect experiment
  #' 
  #' @note THIS FUNCTION PERFORMS AN OLDER VERSION OF BISECT WITHOUT THE 
  #' CUT TRICK. 
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


  if (class(seed) == "numeric") {
    set.seed(seed)
  }
  if (seed == T) {
    seed <- sample(1:10^9, 1)
    print("A seed has been generated")
    dir.create("seed", showWarnings = F)
    write.table(seed, file = "seed/seed.txt", sep = " ")
    set.seed(seed)
  }
  if (dev_opt == T) {
    warning(glue("The developer option has been activated.
                                This means that the ODM environment will not
                                be deleted, and therefore some collusion might
                                occur with following experiments. Poceed with
                                caution
                                "))
  }
  if(exists("ODM_env", envir = globalenv()) != T){
    fit_all_methods(method,...)
  }

  x_list <- runif_on_sphere(n = gen_points, d = ncol(DB) - 2, r = 1) 
  l <- max(sqrt(rowSums(DB[2:(ncol(DB) - 1)]^2)))
  
  registerDoParallel(num_workers)
  
  tic()
  bisection_results <- foreach (i =  1:nrow(x_list), .combine = rbind) %dopar% {
    bisection_results <- bisect( l = l, x = x_list[i, ], method = method, 
                                 verb = verb, check_version = check_version, 
                                 l_val_option = l_val_option)
    hidden_c <- bisection_results[[1]]
    outlier_type <- bisection_results[[2]]

    if (outlier_type %in% c("H1", "H2")) {
      result_point <- hidden_c * (x_list[i, ]) +
        colMeans(DB[2:(ncol(DB) - 1)])
      result_point[length(result_point) + 1] = outlier_type #'doPar has a weird way of dealing with 
                                     #'the outcome of the loops, so we need 
                                     #'to handle the results this way
    }else{result_point = matrix(0,1,ncol(DB)-1)}
    result_point
  }
  exec_time <- toc()
  exec_time <- exec_time$callback_msg
  stopImplicitCluster()
  
  hidden_x_list <- matrix(as.numeric(bisection_results[,1:ncol(x_list)]), 
                          nrow(bisection_results)) #' We added a char into a 
                                                   #' numeric array, so 
                                                   #' everything is a char now
  hidden_x_type <- as.matrix(bisection_results[,(ncol(x_list)+1)])
  hidden_x_list <- hidden_x_list[rowSums(hidden_x_list) != 0, ]
  gen_result <- hog_method(
    DB, gen_points, method, "bisection",
    ODM_env, hidden_x_list, hidden_x_type, exec_time, directions = x_list
  )

  if (dev_opt == F) {
    rm(ODM_env, envir = globalenv())
  }
  return(gen_result)
}
