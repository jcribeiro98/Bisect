

outlier_check_fast <- function(x, verb = FALSE, method = "mahalanobis", ...) {
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
  #' @note This function is the fast version of the other one included in this
  #' file. We do not recommend the use of the slowest version, as the only 
  #' difference between the two is the speed (this version does not approximate
  #' anything, the results are the same, just faster with this version). 
  #'
  #' Arguments:
  #' @param x :(array) Point to check
  #' @param verb : (Logical) Value controlling if it should be verbose or not (
  #'               defaults to FALSE)
  #' @param method: ODM to check with. If no ODM_env is in memory, one would be
  #' created for the selected method. If ODM_env is in memory, this argument 
  #' will be ignored.
  

  # Preface -----------------------------------------------------------------
  
  h2 <- matrix(0, ncol = length(ODM_env$supS), nrow = 1)
  h2[length(ODM_env$supS)] <- 1
  supS <- ODM_env$supS
  index <- matrix(0, ncol = length(ODM_env$supS), nrow = 1)
  if (inference(x, 1:(ncol(DB) - 2), method, ...)) {
    index[length(ODM_env$supS)] = 1
  }
  
  # Checking loops  ---------------------------------------------------------
  
  j <- 0
    for (S in supS) {
      if (set_is_empty(S) != TRUE) {
        if (inference(x, S, method, ...)) {
          index[j] <- 1
          break 
        }
      }
      j <- j + 1}
 
  # Logical assortment ------------------------------------------------------
  
  if (sum(index[1:(length(ODM_env$supS) - 1)]) > 0 &&
      index[length(ODM_env$supS)] == 0) {
    result <- list(0, "H1")
  } else if (isTRUE(all.equal(index, h2))) {
    result <- list(0, "H2")
  } else if (sum(index[1:(length(ODM_env$supS) - 1)]) > 0 &&
             index[length(ODM_env$supS)] == 1) {
    if (verb == T) {
      print(glue("x Outside of bounds"))
    }
    result <- list(1, "OB")
  } else {
    if (verb == T) {
      print(glue("x in the total acceptance area"))
    }
    result <- list(-1, "IL")
  }
  
  return(result)
}


outlier_check <- function(x, verb = F, method = "mahalanobis", ...) {
  #' @title Bisection main function
  #'
  #' @description Defines the function that is going to be used in the bisection
  #' rule during the main execution of the algorithm. The function itself is
  #' really simple, it checks if a point x is an outlier in all of the subspaces
  #' of the total space. f is defined as
  #'
  #'        1, if x is outside of bounds,
  #' f(x)= -1, if x is in the joint acceptance area,
  #'        0, if x is in H1 U H2.
  #'
  #' Arguments:
  #' @param x :(array) Point to check
  #' @param verb : (Logical) Value controlling if it should be verbose or not (
  #'               defaults to FALSE)
  #' @param method: ODM to check with. If no ODM_env is in memory, one would be
  #' created for the selected method. If ODM_env is in memory, this argument 
  #' will be ignored.
  
  
  # Preface -----------------------------------------------------------------
  
  h2 <- matrix(0, ncol = length(ODM_env$supS), nrow = 1)
  h2[length(ODM_env$supS)] <- 1
  supS <- ODM_env$supS
  index <- matrix(0, ncol = length(ODM_env$supS), nrow = 1)
  
  # Checking loops  ---------------------------------------------------------
  
  j <- 0
  for (S in supS) {
    if (set_is_empty(S) != T) {
      if (inference(x, S, method, ...)) {
        index[j] <- 1
      }
    }
    j <- j + 1
  }
  
  # Logical assortment ------------------------------------------------------
  
  if (sum(index[1:(length(ODM_env$supS) - 1)]) > 0 &&
      index[length(ODM_env$supS)] == 0) {
    result <- list(0, "H1")
  } else if (isTRUE(all.equal(index, h2))) {
    result <- list(0, "H2")
  } else if (sum(index[1:(length(ODM_env$supS) - 1)]) > 0 &&
             index[length(ODM_env$supS)] == 1) {
    if (verb == T) {
      print(glue("x Outside of bounds"))
    }
    result <- list(1, "OB")
  } else {
    if (verb == T) {
      print(glue("x in the total acceptance area"))
    }
    result <- list(-1, "IL")
  }
  
  return(result)
}
