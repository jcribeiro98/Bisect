library(sets)
library(glue)
library(dplyr)


DB_gen <- function(db){
#' @title Formating of a database for the bisection method
#' 
#' @description The methods developed in this project utilizes an specific
#' data.frame arquitecture to ease a lot of the calls needed during the 
#' *bisect* method. We included a function to reshape any given matrix/df 
#' into the desired arquitecture. The final df will consist of 
#'            ~id |  y1 |  ... | yn | ~dM
#' Where id is a vector of ids, y1...yn are a renaming of the numerical columns
#' presented on the prior dataframe and dM is the Mahanalobis distance of all 
#' the data.
#' 
#' @note ~ -> Indicates that the df feature might change in further releases.
#' 
#' Arguments
#' @param db :(matrix/data.frame) Database to reshape.

  
  n = ncol(db)
  DB = matrix(0, nrow=nrow(db),ncol= ncol(db) + 1)
  DB = as.data.frame(DB)
  
  name = matrix(0,ncol=ncol(db) + 1,nrow=1)
  for (i in 1:ncol(db)){
    name[i+1] = glue('y{i}')
  }
  name[1]='id'
  colnames(DB) = name
  
  for (i in 1:ncol(db)){
    DB[glue('y{i}')] = db[,i]
  }
  DB['id'] = 1:nrow(DB)
  
  dM = mahalanobis(as.matrix(DB[2:(ncol(db) + 1)],header=T),colMeans(as.matrix(DB[2:(ncol(db) + 1)],header=T)),cov(as.matrix(DB[2:(ncol(db) + 1)],header=T)))
  DB[ncol(db)+2]=dM
  colnames(DB)[ncol(db)+2] = 'dM'
  return(DB)
}

critval <-function(S,verb = F){
  #' @title Critical value given by a normal DB while using the MD. 
  #' 
  #' Arguments
  #' @param S : (set) Set of indices of the  features conforming the DB.
  #' @param verb : (logical) Indicates if the function has to be verbose or not.
  
  
  if (verb == T){
    print(glue('cutoff value on dM = {qchisq(0.95, length(S))} \n'))}
  return(qchisq(0.95, length(S)))}


distmah <- function(S,x){
  #' @title Calculation of the Mah. distance.
  #' 
  #' @description Calculates the Mah. distance given a value x, and a set of 
  #'  indices of the features conforming the space.
  #'  
  #'  Arguments
  #'  @param x : (array) Value to which calculate its MD wrt the mean.
  #'  @param S : (set) Set of indices which conforms the subspace.
  
  
  sS = c(0)
  j = 1
  for (i in S){
    sS[j] = glue('y{i}')
    j = j + 1
  }
  x = as.data.frame(t(x))

  
  return(mahalanobis(x %>% select(colnames(x[sS])), 
                     colMeans(DB %>% select(colnames(DB[sS]))), 
                     cov(DB %>% select(colnames(DB[sS])))))
}