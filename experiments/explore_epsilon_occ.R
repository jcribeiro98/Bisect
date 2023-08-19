source("experiments/one_class_classification.R")
library(R.utils)

run_time <- function(dataset, epsilons){
  run_time_vector <- array(0, length(epsilons))
  i = 0
  for (eps in epsilons){
    loader <- match.fun(glue("load_{dataset}"))
    data = loader(tsplit = .8, seed = 47385123)
    train_data = data[[1]]; train_labels = data[[2]] 
    test_data = data[[3]]; test_labels = data[[4]]
    set.seed(47385123)
    hog_model = main_hidden(gen_points = 100,
                            method = "pyod_LOF",
                            num_workers = 1, until_gen = TRUE,
                            eps = eps, dev_opt = TRUE)
    i = i + 1
    run_time_vector[i] = as.numeric(gsub("[secelapsed ]","", hog_model$exec_time))
  }
  rm(ODM_env, envir = globalenv())
  return(run_time_vector)
}

datasets = c("WILT","PIMA","STAMPS","PARKINSON","HEARTDISEASE","ANNTHYROID",
             "CARDIOTOCOGRAPHY","IONOSPHERE","WPBC","SPAMBASE","ARRYTHMIA")

run_time_all <- function(datasets,epsilons = c(0.1,0.5,0.75)){
  for (dataset in datasets){ 
    exp_res = run_time(glue("{toupper(dataset)}"), epsilons = epsilons) 
    write.csv(exp_res,glue("experiments/explore_epsilon/OCC/{tolower(dataset)}_exp_{Sys.time()}.csv"))
  }
}

run_time_all(datasets)