
# Balancing with cWGAN ----------------------------------------------------
source("experiments/new_balancing_experiments.R")
use_condaenv("cGAN")
cgan_module = import_from_path("wgan.imblearn", path = "experiments/src_extra_detectors/GANbalanced/")


## Functions ---------------------------------------------------------------

cwgan_exp <- function(dataset, tsplit=.2,seed_vector = c(47385123, 12345, 54321, 
                                                         010101, 121212, 19, 33),
                      out_ratio = "default",
                      meths = c("rf")){
  #'@title: Experiment launcher for the cWGAN oversampling method
  #'@description  Launch an experiment on the given dataset with multiple
  #'repetition using the given vector of seeds. 
  #'
  #' Arguments:
  #' @param dataset: Chose a dataset given as a loader function
  #' @param tsplit: Split for the train-test
  #' @param seed_vector: Vector of seeds for multiple reproductions of the 
  #' experiment
  #' @param out_ratio: Outlier ratio of the downloaded datasets
  #' @param meths: Two-class classification method to oversample. Options are:
  #'               -svmLinear2: SVM with linear kernel
  #'               -Any classifier on the "caret" package. For example, 'rf'
  #'               which is random forest.
  #'

  loader <- match.fun(glue("load_{dataset}"))
  auc_matrix_smote = matrix(0, nrow = length(meths), ncol = length(seed_vector))
  rownames(auc_matrix_smote) = meths
  i = 1
  for (meth in meths){
    j = 1
    for (seed in seed_vector){
      if (out_ratio == "default"){
        data = loader(tsplit = tsplit, seed = seed)}else{
          data = loader(tsplit = tsplit, out_ratio = out_ratio, seed = seed)
        }
      train_data = data[[1]]
      train_labels = data[[2]]
      test_data = data[[3]]
      test_labels = data[[4]]
      set.seed(seed)
      
      cGAN = cgan_module$GANbalancer(idx_cont = np_array(0:(ncol(train_data)-1)), 
                                     categorical = NULL,auxiliary=TRUE,
                                     gan_architecture="fisher",
                                     generator_input= ncol(train_data), 
                                     generator_layers=c(50L,50L), 
                                     critic_layers=c(50L,50L), 
                                     layer_norm=TRUE)
      resampled_dataset = cGAN$fit_resample(train_data,np_array(train_labels,dtype = "int32"))
      
      
      final_labels = resampled_dataset[[2]]
      final_data = resampled_dataset[[1]]
      if(meth == "svmLinear2"){
        rf_classifier <-svm(final_data,final_labels, type = "C-classification", 
                            kernel = "radial", probability = TRUE)
        predictions <- as.numeric(attr(predict(rf_classifier, test_data,  probability = TRUE), "probabilities")[,1])
      }else{
        rf_classifier <- train(final_data,as.factor(final_labels),
                               method = meth)
        predictions <- predict(rf_classifier, test_data, type = "prob")[,1]
      }
      auc_matrix_smote[i,j] = c(auc(test_labels, predictions)) 
      j = j + 1
    }
    i = i + 1
  }
  return(auc_matrix_smote)
}


run_all_cwgan <- function(datasets = c("HEARTDISEASE", "PIMA",
                                          "STAMPS", "WPBC","PARKINSON",
                                          "IONOSPHERE", "SPAMBASE", "CARDIOTOCOGRAPHY", 
                                          "ANNTHYROID", "PAGEBLOCKS", "ARRYTHMIA",
                                          "WILT"), seed_vector = c(47385123L, 12345L, 54321L, 
                                                                   010101L, 121212L, 19L, 33L)
){
  #' @title Runs a cWGAN with default param. over all the datasets
  #' @description Runs a cWGAN with default param. over all the datasets for the
  #' experiments. The vector of seeds control how many repetitions will be performed
  #' @param datasets: Array of datasets
  #' @param seed_vector: Vector of seeds for multiple repetitions
  
  for (dataset in datasets){ 
    assign(glue("{tolower(dataset)}_cGAN_results"), 
           cwgan_exp(glue("{toupper(dataset)}"),
                     seed_vector = seed_vector), 
           envir = globalenv())
    save.image(glue("experiments/new_balancing_experiments/RData/cGAN_exp_{Sys.time()}.RData"))
    write.csv(get(glue("{tolower(dataset)}_cGAN_results")),
              glue("experiments/new_balancing_experiments/CSVs/{tolower(dataset)}_cGAN_results_{Sys.time()}.csv"))
  }
}

get_cgan_unsup = function(datasets, method = "rf"){
  #' @details Auxiliary function to fetch the results as a matrix from the 
  #' enviroment.

  res = c()
  for (dataset in datasets){
    res = rbind(res, median(get(glue("{tolower(dataset)}_cGAN_results"))[method,]))
  }
  rownames(res) = datasets
  return(res)
} 


## Experiments ------------------------------------------------------------

run_all_cwgan()
datasets = c("WILT", "PIMA", "STAMPS", "PAGEBLOCKS", "HEARTDISEASE", "ANNTHYROID",
             "CARDIOTOCOGRAPHY", "PARKINSON", "IONOSPHERE", "WPBC", "SPAMBASE",
             "ARRYTHMIA")
get_cgan_unsup(datasets)
