library(reticulate)
library(progress)

hyperbox <- function(data,nsample){
  #' @title Hyperbox Generator
  #' @description Main function of the Hyperbox generator, as described in 
  #' the related article.
  #' Arguments:
  #' @param data: Data to generate outliers for
  #' @param nsample: Number of outliers to generate
  
  LOF = import("pyod.models.lof")
  LOF = LOF$LOF()
  LOF$fit(data)
  gen_data_matrix = matrix(ncol = ncol(data))
  gen_data_matrix = gen_data_matrix[!rowSums(is.na(gen_data_matrix)),] #First row are NAs
  
  max_values_in_feature = matrix(0,1,ncol = ncol(data))
  min_values_in_feature = max_values_in_feature
  for (j in 1:ncol(data)){
    max_values_in_feature[j] = max(data[,j])
    min_values_in_feature[j] = min(data[,j])
  }
  pb <- progress_bar$new(total = nsample)   
  while(nrow(gen_data_matrix) < nsample){
    candidate = array(0, ncol(data))
    for(j in 1:ncol(data)){
      candidate[j] = runif(1,min_values_in_feature[j],max_values_in_feature[j])
    }
    if(LOF$predict(t(candidate)) == 1){
      gen_data_matrix = rbind(gen_data_matrix,candidate)
      pb$tick()
    }
  }
  row.names(gen_data_matrix) <- NULL
  colnames(gen_data_matrix) <- colnames(data)
  return(gen_data_matrix)
}
