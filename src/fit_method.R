 
source("src/utils.R")

fit <- function(method, S){
  ##Pending description##
  
  
  #Methods:
  
  if (method == "mahalanobis"){
    result = function(x){distmah(S,x) > critval(S,verb=F)}
  }
  
  return(result)  
}