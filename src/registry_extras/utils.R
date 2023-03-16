library(sets)
library(glue)
library(dplyr)


DB_gen <- function(db, true_inliers = FALSE, method = "pyod_LOF",...){
#' @title Formating of a database for the bisection method
#' 
#' @description The methods developed in this project utilizes an specific
#' data.frame architecture to ease a lot of the calls needed during the 
#' *bisect* method. We included a function to reshape any given matrix/df 
#' into the desired architecture. The final df will consist of 
#'            ~id |  y1 |  ... | yn | ~dM
#' Where id is a vector of ids, y1...yn are a renaming of the numerical columns
#' presented on the prior dataframe and dM is the Mahalanobis distance of all 
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
  
  out = array(0, dim = nrow(DB))
  if(true_inliers == T){
    DB <<- DB
    fit_all_methods(method,...)
  }
  if(exists("ODM_env", envir = globalenv())){
    model = ODM_env[[glue("method{set_names(1:(ncol(DB)-2))}")]]
    for ( i in 1:nrow(DB)){
      x = DB[i,]
      names = names(x)
      x = as.matrix(x)
      dim(x) <- NULL
      
      names(x) = names
      if(model(x)){
          out[i] = 1
      }
    }
  }
    
  DB[ncol(db)+2]=out
  colnames(DB)[ncol(db)+2] = 'Out'
  DB <<- DB
  return(DB)
}

critval <- function(S,verb = F){
  #' @title Critical value given by a normal DB while using the MD. 
  #' 
  #' Arguments
  #' @param S : (set) Set of indices of the  features conforming the DB.
  #' @param verb : (logical) Indicates if the function has to be verbose or not.
  
  
  if (verb == T){
    print(glue('cutoff value on dM = {qchisq(0.95, length(S))} \n'))}
  return(qchisq(0.95, length(S)))}


set_names <- function(s, sep = "_"){  
  #' @title Names of a set
  #' 
  #' @description Given an array or list, the function obtains a character 
  #' containing the name (as a character) of each element, separated by sep.
  #' 
  #' Arguments:
  #' @param s : (array or list) Array whose names want to be fetched
  #' @param sep : (character) Separator used in between the names (defaults to 
  #'              "_")
  
  
  index = as.vector(matrix(0,1,length(s)))
  j = 1
  for (i in s){index[j] = i; j = j + 1}
  result = paste(as.character(index), collapse = sep)
  return(result)
}

set_subspace_grab <- function(S){
  sS = c(0)
  j = 1
  for (i in S){
    sS[j] = glue('y{i}')
    j = j + 1
  }
  return(sS)
}


distmah <- function(S,x){
  #' @title Calculation of the Mah. distance.
  #' 
  #' @description Calculates the Mah. distance given a value x, and a set of 
  #'  indices of the features conforming the space.
  #'  
  #'  Arguments
  #'  @param x : (array) Value to which calculate its MD wrt the mean.
  #'  @param S : (set) Set of indices which conforms the subspace.
  
  
  sS = set_subspace_grab(S)
  x = as.data.frame(t(x))

  
  return(mahalanobis(x %>% select(colnames(x[sS])), 
                     colMeans(DB %>% select(colnames(DB[sS]))), 
                     cov(DB %>% select(colnames(DB[sS])))))
}

init_python <- function(){
  if (!("hidden_out" %in% conda_list()[,1])){
    stop("No virtual enviroment by the name *hidden_out* has been found.
  Please, read the docs to learn how to install the python methods.
  If you have succesfully managed the installation, please, restart the R session 
  and try again.")
  }
  use_condaenv("hidden_out")
  np <- import("numpy")
  
}


