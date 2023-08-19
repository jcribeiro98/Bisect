# Balancing task --------------------------------------------------
library(miceadds)
source.all("src/HOGen/")
source("experiments/lof_based_generator.R")
library(caret)
library(pROC)
library(ROCR)
library(DMwR)
library(farff)
library(e1071)
library(smotefamily)
use_condaenv("hidden_out")

## Auxiliary Functions ---------------------------------------------------------------
generate_balancing <- function(seed_set = 47385123,
                               train_data, train_labels, test_data, test_labels, synth_data,synth_label,
                               meth = "rf"){
  #' @title Trains a given method with a balanced dataset
  #' @description Function to handle the balancing experiment. Given a two-class
  #' method, train_data, synth_data and test_data; it trains the balanced classifier
  #' and then tests on the test data, outputting the auc of the training. 
  #' @param seed_set: Seed for the single experiment
  #' @param train_data: Train data to balance
  #' @param train_labels: Labels for the train data
  #' @param test_data: Test data to test with
  #' @param test_labels: Labels for the test data
  #' @param synth_data: Synthetic data generated with an oversampling method
  #' @param synth_labels: Labels for the synthetic data
  
  i = 0
    tryCatch({
      i = i + 1
      set.seed(seed_set)
      index = sample(1:nrow(synth_data), length(train_labels[train_labels == 0]) - length(train_labels[train_labels == 1]))
      synth_data_dummy = synth_data[index,]
      train_data_dummy = rbind(train_data, synth_data_dummy)
      rownames(train_data_dummy) = NULL #Renumbers the dataframe
      train_labels_dummy = c(train_labels,synth_label[index])
      rm(index)
      
      set.seed(seed_set)
      if(meth == "svmLinear2"){
        rf_classifier <-svm(train_data_dummy,train_labels_dummy, type = "C-classification", 
                            kernel = "radial", probability = TRUE)
        predictions <- as.numeric(attr(predict(rf_classifier, test_data,  probability = TRUE), "probabilities")[,1])
        }else{
      rf_classifier <- train(train_data_dummy, 
                             factor(train_labels_dummy), method = meth)
      
      predictions <- predict(rf_classifier, test_data, type = "prob")[,1]
      }
      auc_result <- auc(test_labels, predictions)
      }
    ,
    error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
  return(auc_result)
}

## Detection methods  ----------------------------------------------------
### MBisect -----------------------------------------------------------------

bisect_exp <- function(dataset,tsplit = .2, out_ratio = "default", methods, 
                       seed_vector = c(47385123L, 12345L, 54321L, 
                                       010101L, 121212L, 19L, 33L),
                       meths = c("rf","svmLinear2","knn")){
  #'@title: Experiment launcher for the Bisect oversampling method
  #'@description  Launch an experiment on the given dataset with multiple
  #'repetition using the given vector of seeds. 
  #'
  #' Arguments:
  #' @param dataset: Chose a dataset given as a loader function
  #' @param tsplit: Split for the train-test
  #' @param methods: Origin selection method to test
  #' @param seed_vector: Vector of seeds for multiple reproductions of the 
  #' experiment
  #' @param out_ratio: Outlier ratio of the downloaded datasets
  #' @param meths: Two-class classification method to oversample. Options are:
  #'               -svmLinear2: SVM with linear kernel
  #'               -Any classifier on the "caret" package. For example, 'rf'
  #'               which is random forest.
  #'
  
  loader <- match.fun(glue("load_{dataset}"))
  auc_matrix_bisect = matrix(0,nrow = length(meths)*length(methods),
                             ncol = length(seed_vector))
  auc_matrix_bisect_names = array(0, length(meths)*length(methods))
  
  j = 1
  for (seed in seed_vector){
    i = 1
    if (out_ratio == "default"){
      data = loader(tsplit = tsplit, seed = seed)}else{
        data = loader(tsplit = tsplit, out_ratio = out_ratio, seed = seed)
      }
    train_data = data[[1]]
    train_labels = data[[2]]
    test_data = data[[3]]
    test_labels = data[[4]]
    for(method in methods){
      hog_model = main_multibisect(gen_points = 2*nrow(train_data),
                                   method = "pyod_LOF", seed = seed, 
                                   num_workers = 1, l_val_option = "not_fixed",
                                   type = method, dev_opt = TRUE)
      
      synth_data = as.data.frame(hog_model$ho_list)
      synth_label = rep(1, nrow(synth_data))
      names(synth_data) = names(train_data) #Just so it can bind in the next line
      for(meth in meths){ #Same hidden outliers have been put in the training set
        auc_matrix_bisect[i,j] = generate_balancing(seed,train_data,train_labels,
                                                    test_data,test_labels,
                                                    synth_data,synth_label,meth)
        auc_matrix_bisect_names[i] = glue("{method}_{meth}")
        i = i + 1
      }
    }
    j = j + 1 #For each seed, the same subspaces have been used in the definition of hidden outlier
    rm(ODM_env, envir = globalenv())
    rm(proba_vector, envir = globalenv())
  }
  rownames(auc_matrix_bisect) = auc_matrix_bisect_names
  return(auc_matrix_bisect)
}

### Hidden -----------------------------------------------------------------

hidden_exp <- function(dataset, tsplit = .2, epsilons, out_ratio = "default",
                       seed_vector = c(47385123L, 12345L, 54321L, 
                                       010101L, 121212L, 19L, 33L), 
                       meths = c("rf","svmLinear2","knn")
                       ){
  #'@title: Experiment launcher for the Hidden oversampling method
  #'@description  Launch an experiment on the given dataset with multiple
  #'repetition using the given vector of seeds. 
  #'
  #' Arguments:
  #' @param dataset: Chose a dataset given as a loader function
  #' @param tsplit: Split for the train-test
  #' @param epsilons: Epsilons to test
  #' @param seed_vector: Vector of seeds for multiple reproductions of the 
  #' experiment
  #' @param out_ratio: Outlier ratio of the downloaded datasets
  #' @param meths: Two-class classification method to oversample. Options are:
  #'               -svmLinear2: SVM with linear kernel
  #'               -Any classifier on the "caret" package. For example, 'rf'
  #'               which is random forest.
  #'
  
  loader <- match.fun(glue("load_{dataset}"))
  auc_matrix_hidden = matrix(0,nrow = length(meths)*length(epsilons),
                             ncol = length(seed_vector))
  auc_matrix_hidden_names = array(0, length(meths)*length(epsilons))
  
  j = 1
  for(seed in seed_vector){
    i = 1
    if (out_ratio == "default"){
      data = loader(tsplit = tsplit, seed = seed)}else{
        data = loader(tsplit = tsplit, out_ratio = out_ratio, seed = seed)
      }
    train_data = data[[1]]
    train_labels = data[[2]]
    test_data = data[[3]]
    test_labels = data[[4]]
    for(epsilon in epsilons){
      set.seed(seed)
      hog_model = main_hidden(gen_points = 2*nrow(train_data),
                              method = "pyod_LOF", num_workers = 1,
                              eps = epsilon, dev_opt = TRUE, until_gen = TRUE)
      
      synth_data = as.data.frame(hog_model$ho_list)
      synth_label = rep(1, nrow(synth_data))
      names(synth_data) = names(train_data) #Just so it can bind in the next line
      for (meth in meths){
        auc_matrix_hidden[i,j] = generate_balancing(seed,train_data,train_labels,
                                                    test_data,test_labels,synth_data,
                                                    synth_label,meth)
        auc_matrix_hidden_names[i] = glue("{meth}_eps{epsilon}")
        i = i + 1
      }
    }
    j = j + 1
    rm(ODM_env, envir = globalenv())
  }
  rownames(auc_matrix_hidden) = auc_matrix_hidden_names
  return(auc_matrix_hidden)
}

### Nothing --------------------------------------------------------------

nothing_exp <- function(dataset, tsplit = .2, out_ratio = "default",
                        seed_vector, methods = c("rf","svmLinear2", "knn")){
  #' @title: Experiment launcher for no oversampling method
  #' @description  Launch an experiment on the given dataset with multiple
  #' repetition using the given vector of seeds. 
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
  auc_matrix_nothing = matrix(0,nrow = length(methods), ncol = length(seed_vector))
  rownames(auc_matrix_nothing) = methods
  i = 1
  for (meth in methods){
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
      if(meth == "svmLinear2"){
        rf_classifier <-svm(train_data,train_labels, type = "C-classification", 
                            kernel = "radial", probability = TRUE)
        predictions <- as.numeric(attr(predict(rf_classifier, test_data,  probability = TRUE), "probabilities")[,1])
      }else{
        rf_classifier <- train(factor(train_labels)~.,data = cbind(train_data,train_labels), method = meth)
        predictions <- predict(rf_classifier, test_data, type = "prob")[,1]
      }
      print(glue("working with method: {meth}"))
      auc_matrix_nothing[i,j] <- c(auc(test_labels, predictions))  
      j = j + 1
    }
    i = i + 1
  }
  return(auc_matrix_nothing)
}



### SMOTE --------------------------------------------------------------

smote_exp <- function(dataset, tsplit=.2,seed_vector = c(47385123, 12345, 54321, 
                                                         010101, 121212, 19, 33),
                      out_ratio = "default",
                      meths = c("rf","svmLinear2", "knn"), K = 3){
  #' @title: Experiment launcher for the SMOTE oversampling method
  #' @description  Launch an experiment on the given dataset with multiple
  #' repetition using the given vector of seeds. 
  #'
  #' Arguments:
  #' @param dataset: Chose a dataset given as a loader function
  #' @param tsplit: Split for the train-test
  #' @param K: Number of NN
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
      synth_data = SMOTE(train_data, train_labels, K = K,
                         dup_size = length(train_labels[train_labels == 0])/
                           length(train_labels[train_labels == 1]))$syn_data
      synth_label = as.vector(synth_data["class"]$class)
      synth_data = synth_data[,-c(ncol(synth_data))]
      names(synth_data) = names(train_data) #Just so it can bind in the next line
      if(meth == "svmLinear2"){
        rf_classifier <-svm(rbind(train_data,synth_data),c(train_labels,synth_label), type = "C-classification", 
                            kernel = "radial", probability = TRUE)
        predictions <- as.numeric(attr(predict(rf_classifier, test_data,  probability = TRUE), "probabilities")[,1])
      }else{
        rf_classifier <- train(rbind(train_data,synth_data), 
                               factor(c(train_labels,synth_label)),
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

dbsmote_exp <- function(dataset, tsplit=.2,seed_vector = c(47385123, 12345, 54321, 
                                                           010101, 121212, 19, 33),
                        out_ratio = "default",
                        meths = c("rf","svmLinear2", "knn")){
  #' @title: Experiment launcher for the DBSMOTE oversampling method
  #' @description  Launch an experiment on the given dataset with multiple
  #' repetition using the given vector of seeds. 
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
  auc_matrix_dbsmote = matrix(0, nrow = length(meths), ncol = length(seed_vector))
  rownames(auc_matrix_dbsmote) = meths
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
      synth_data = DBSMOTE(train_data, train_labels,
                           dupSize = length(train_labels[train_labels == 0])/
                             length(train_labels[train_labels == 1]))$syn_data
      synth_label = as.vector(synth_data["class"]$class)
      synth_data = synth_data[,-c(ncol(synth_data))]
      names(synth_data) = names(train_data) #Just so it can bind in the next line
      if(meth == "svmLinear2"){
        rf_classifier <-svm(rbind(train_data,synth_data),c(train_labels,synth_label), type = "C-classification", 
                            kernel = "radial", probability = TRUE)
        predictions <- as.numeric(attr(predict(rf_classifier, test_data,  probability = TRUE), "probabilities")[,1])
      }else{
        rf_classifier <- train(rbind(train_data,synth_data), 
                               factor(c(train_labels,synth_label)),
                               method = meth)
        predictions <- predict(rf_classifier, test_data, type = "prob")[,1]
      }
      auc_matrix_dbsmote[i,j] = c(auc(test_labels, predictions)) 
      j = j + 1
    }
    i = i + 1
  }
  return(auc_matrix_dbsmote)
}

adasyn_exp <- function(dataset, tsplit=.2,seed_vector = c(47385123, 12345, 54321, 
                                                          010101, 121212, 19, 33),
                       out_ratio = "default",
                       meths = c("rf","svmLinear2", "knn"), K = 3){
  #' @title: Experiment launcher for the Adasyn oversampling method
  #' @description  Launch an experiment on the given dataset with multiple
  #' repetition using the given vector of seeds. 
  #'
  #' Arguments:
  #' @param dataset: Chose a dataset given as a loader function
  #' @param tsplit: Split for the train-test
  #' @param K: Number of NN
  #' @param seed_vector: Vector of seeds for multiple reproductions of the 
  #' experiment
  #' @param out_ratio: Outlier ratio of the downloaded datasets
  #' @param meths: Two-class classification method to oversample. Options are:
  #'               -svmLinear2: SVM with linear kernel
  #'               -Any classifier on the "caret" package. For example, 'rf'
  #'               which is random forest.
  #'
  
  loader <- match.fun(glue("load_{dataset}"))
  auc_matrix_adasyn = matrix(0, nrow = length(meths), ncol = length(seed_vector))
  rownames(auc_matrix_adasyn) = meths
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
      synth_data = ADAS(train_data, train_labels, K = K)$syn_data
      synth_label = as.vector(synth_data["class"]$class)
      synth_data = synth_data[,-c(ncol(synth_data))]
      names(synth_data) = names(train_data) #Just so it can bind in the next line
      if(meth == "svmLinear2"){
        rf_classifier <-svm(rbind(train_data,synth_data),c(train_labels,synth_label), type = "C-classification", 
                            kernel = "radial", probability = TRUE)
        predictions <- as.numeric(attr(predict(rf_classifier, test_data,  probability = TRUE), "probabilities")[,1])
      }else{
        rf_classifier <- train(rbind(train_data,synth_data), 
                               factor(c(train_labels,synth_label)),
                               method = meth)
        predictions <- predict(rf_classifier, test_data, type = "prob")[,1]
      }
      auc_matrix_adasyn[i,j] = c(auc(test_labels, predictions)) 
      j = j + 1
    }
    i = i + 1
  }
  return(auc_matrix_adasyn)
}

### HYPERBOX -----------------------------------------------------------

run_hyperbox <- function(dataset,
                         seed_vector = c(47385123, 12345, 54321, 010101, 
                                         121212, 19, 33),tsplit = .2, m_vector = c(1)){
  #' @title: Experiment launcher for the Hyperbox oversampling method
  #' @description  Launch an experiment on the given dataset with multiple
  #' repetition using the given vector of seeds. 
  #'
  #' Arguments:
  #' @param dataset: Chose a dataset given as a loader function
  #' @param tsplit: Split for the train-test
  #' @param m_vector: Number of points to generate (in multiples of the train 
  #' data number of points). If an array is given, multiple experiments are
  #' performed
  #' @param seed_vector: Vector of seeds for multiple reproductions of the 
  #' experiment
  #' @param out_ratio: Outlier ratio of the downloaded datasets
  #' @param meths: Two-class classification method to oversample. Options are:
  #'               -svmLinear2: SVM with linear kernel
  #'               -Any classifier on the "caret" package. For example, 'rf'
  #'               which is random forest.
  #'
  
  loader <- match.fun(glue("load_{dataset}"))
  
  for (seed in seed_vector){
    for(m in m_vector){
      data = loader(tsplit = tsplit, seed = seed)
      train_data = data[[1]]; train_labels = data[[2]] 
      test_data = data[[3]]; test_labels = data[[4]]
      synth_data <- as.data.frame(hyperbox(train_data,floor(m*nrow(train_data))))
      names(synth_data) <- names(train_data)  
      assign(glue("{m}_{tolower(dataset)}_hyperbox_results_{seed}"), 
             generate_casting(seed = seed,
                              train_data, train_labels, test_data, test_labels,
                              synth_data))
    }
    rm(ODM_env, envir = globalenv())
    rm(proba_vector, envir = globalenv())
  }
  
  for(m in m_vector){
    result = c()
    for(seed in seed_vector){
      result = cbind(result, get(glue("{m}_{tolower(dataset)}_hyperbox_results_{seed}")))
    }
    rownames(result) <- c("ROC", "Prec-Rec") 
    assign(glue("{m}_{tolower(dataset)}_hyperbox_results"),result, envir = globalenv())
  }
}


## Experiment runners  ---------------------------------------------------

### Bisect -----------------------------------------------------------

run_all_bisect <- function(methods = c("weighted"), 
                           out_ratio = "default",
                           datasets = c("GLASS", "HEARTDISEASE", "HEPATITITS",
                                        "WBC", "PIMA", "SHUTTLE", "STAMPS", 
                                        "WPBC", "ARRYTHMIA"), 
                           seed_vector = c(47385123L, 12345L, 54321L, 010101L, 121212L, 19L, 33L)){
  #' @details Runs all experiments using Bisect
  
  for (dataset in datasets){
    assign(glue("{tolower(dataset)}_bisect_results"), 
           bisect_exp(glue("{toupper(dataset)}"), methods = methods, 
                      out_ratio = out_ratio,seed_vector = seed_vector), 
           envir = globalenv())
    save.image(glue("experiments/new_balancing_experiments/RData/bisect_exp_{Sys.time()}.RData")) 
  }
}


### Hidden -----------------------------------------------------------

run_all_hidden <- function(eps = c(.1),seed_vector = c(47385123L, 12345L, 54321L, 
                                                       010101L, 121212L, 19L, 
                                                       33L),
                           out_ratio = "default",datasets){
  #' @details Runs all experiments using Hidden
  
  for (dataset in datasets){
    assign(glue("{tolower(dataset)}_hidden_results"), 
           hidden_exp(glue("{toupper(dataset)}"),  epsilons = eps, 
                      seed_vector = seed_vector, out_ratio = out_ratio), 
           envir = globalenv())
    save.image(glue("experiments/new_balancing_experiments/RData/{eps}_hidden_exp_{Sys.time()}.RData")) 
  }
}

### Synthetic Oversampling -----------------------------------------------------------

run_all_synthetic_oversamp <- function(K = 5,seed_vector = c(47385123L, 12345L, 
                                                             54321L, 010101L, 
                                                             121212L, 19L, 33L),
                                       out_ratio = "default",datasets){
  #' @details Runs all experiments using all synth. minority oversampling methods.
  
  for (dataset in datasets){
    assign(glue("{tolower(dataset)}_smote_results"), 
           smote_exp(glue("{toupper(dataset)}"),
                     seed_vector = seed_vector, out_ratio = out_ratio, K = K), 
           envir = globalenv())
    assign(glue("{tolower(dataset)}_dbsmote_results"), 
           dbsmote_exp(glue("{toupper(dataset)}"),
                       seed_vector = seed_vector, out_ratio = out_ratio), 
           envir = globalenv())
    assign(glue("{tolower(dataset)}_adasyn_results"), 
           adasyn_exp(glue("{toupper(dataset)}"),
                      seed_vector = seed_vector, out_ratio = out_ratio, K = K), 
           envir = globalenv())
    save.image(glue("experiments/new_balancing_experiments/RData/smote_exp_{Sys.time()}.RData"))
  }
}


### Hyperbox -------------------------------------------------------------

run_all_hyperbox <- function(datasets = c("HEARTDISEASE", "PIMA",
                                          "STAMPS", "WPBC",
                                          "IONOSPHERE", "SPAMBASE", 
                                          "CARDIOTOCOGRAPHY", 
                                          "ANNTHYROID", "PAGEBLOCKS", 
                                          "ARRYTHMIA", "WILT"), 
                             seed_vector = c(47385123L, 12345L, 54321L, 010101L,
                                             121212L, 19L, 33L)){
  #' @details Runs all experiments using HB
  
  for (dataset in datasets){ 
    run_hyperbox(glue("{toupper(dataset)}"), 
                 seed_vector = seed_vector) 
    save.image(glue("experiments/new_balancing_experiments/RData/{m}_hyperbox_exp_{Sys.time()}.RData"))
  }
}



## Dataset loaders ---------------------------------------------------------

### HEART DISEASE -----------------------------------------------------------
load_HEARTDISEASE <- function(tsplit = .2, out_ratio = "02", hog_exp = TRUE, seed = 47385123){
  #' @title Loads the named dataset from the dataset directories. To download the
  #' datasets follow the instructions of the repository. It outputs a list
  #' containing the train set, the test set, the train labels and the test labels
  #' Arguments:
  #' @param tsplit: Split of train-test
  #' @param out_ratio: Ratio of outliers for the train test
  #' @param hog_exp: Selector to check if the experiment is using hidden outliers.
  #' It allows extra computations (generates the DB object)
  #' @param seed: Seed for the split
   
  data = readARFF(glue("databases/HeartDisease_withoutdupl_norm_{out_ratio}_v10.arff"))
  true_size_data = readARFF("databases/true_size/HeartDisease_withoutdupl_norm_44.arff")
  
  target = as.vector(data[,15])
  target[target == "yes"] = rep(1, length(target[target == "yes"]))
  target[target == "no"] = rep(0, length(target[target == "no"]))
  set.seed(seed)
  index = as.vector(createDataPartition(target, p = tsplit, list = FALSE))
  train_data <- data[index,-c(14,15)]
  train_labels <- target[index]
  test_data <- rbind(data[-index,-c(14,15)],filter(true_size_data, !id %in% data$id)[,-c(1,15)])
  test_labels <- c(target[-index],rep(1,length(filter(true_size_data, !id %in% data$id)[,15])))
  
  if(hog_exp == TRUE){
    DB_gen(train_data)
    rm(ODM_env, envir = globalenv())
  }
  
  return(list(train_data,train_labels,test_data,test_labels))
}


### STAMPS ---------------------------------------------------------------
load_STAMPS <- function(tsplit = .2, out_ratio = "02", hog_exp = TRUE, seed = 47385123){
  #' @title Loads the named dataset from the dataset directories. To download the
  #' datasets follow the instructions of the repository. It outputs a list
  #' containing the train set, the test set, the train labels and the test labels
  #' Arguments:
  #' @param tsplit: Split of train-test
  #' @param out_ratio: Ratio of outliers for the train test
  #' @param hog_exp: Selector to check if the experiment is using hidden outliers.
  #' It allows extra computations (generates the DB object)
  #' @param seed: Seed for the split
  
  data = readARFF(glue("databases/Stamps_withoutdupl_norm_{out_ratio}_v10.arff"))
  true_size_data = readARFF("databases/true_size/Stamps_withoutdupl_norm_09.arff")
  
  target = as.vector(data[,11])
  target[target == "yes"] = rep(1, length(target[target == "yes"]))
  target[target == "no"] = rep(0, length(target[target == "no"]))
  set.seed(seed)
  index = as.vector(createDataPartition(target, p = tsplit, list = FALSE))
  train_data <- data[index,-c(10,11)]
  train_labels <- target[index]
  test_data <- rbind(data[-index,-c(10,11)],filter(true_size_data, !id %in% data$id)[,-c(10,11)])
  test_labels <- c(target[-index],rep(1,length(filter(true_size_data, !id %in% data$id)[,11])))
  
  if(hog_exp == TRUE){
    DB_gen(train_data)
    rm(ODM_env, envir = globalenv())
  }
  
  return(list(train_data,train_labels,test_data,test_labels))
}

### PIMA -----------------------------------------------------------------
load_PIMA <- function(tsplit = .2, out_ratio = "02", hog_exp = TRUE, seed = 47385123){
  #' @title Loads the named dataset from the dataset directories. To download the
  #' datasets follow the instructions of the repository. It outputs a list
  #' containing the train set, the test set, the train labels and the test labels
  #' Arguments:
  #' @param tsplit: Split of train-test
  #' @param out_ratio: Ratio of outliers for the train test
  #' @param hog_exp: Selector to check if the experiment is using hidden outliers.
  #' It allows extra computations (generates the DB object)
  #' @param seed: Seed for the split
  
  data = readARFF(glue("databases/Pima_withoutdupl_norm_{out_ratio}_v10.arff"))
  true_size_data = readARFF("databases/true_size/Pima_withoutdupl_norm_35.arff")
  
  target = as.vector(data[,10])
  target[target == "yes"] = rep(1, length(target[target == "yes"]))
  target[target == "no"] = rep(0, length(target[target == "no"]))
  set.seed(seed)
  index = as.vector(createDataPartition(target, p = tsplit, list = FALSE))
  train_data <- data[index,-c(9,10)]
  train_labels <- target[index]
  test_data <- rbind(data[-index,-c(9,10)],filter(true_size_data, !id %in% data$id)[,-c(1,10)])
  test_labels <- c(target[-index],rep(1,length(filter(true_size_data, !id %in% data$id)[,10])))
  
  if(hog_exp == TRUE){
    DB_gen(train_data)
    rm(ODM_env, envir = globalenv())
  }
  
  return(list(train_data,train_labels,test_data,test_labels))
}

### WPBC -----------------------------------------------------------------
load_WPBC <- function(tsplit = .2, out_ratio = "02", hog_exp = TRUE, seed = 47385123){
  #' @title Loads the named dataset from the dataset directories. To download the
  #' datasets follow the instructions of the repository. It outputs a list
  #' containing the train set, the test set, the train labels and the test labels
  #' Arguments:
  #' @param tsplit: Split of train-test
  #' @param out_ratio: Ratio of outliers for the train test
  #' @param hog_exp: Selector to check if the experiment is using hidden outliers.
  #' It allows extra computations (generates the DB object)
  #' @param seed: Seed for the split

  out_ratio = as.integer(out_ratio)/100
  true_size_data = readARFF("databases/WPBC_withoutdupl_norm.arff")
  set.seed(47385123) #So we always get the same downsampling
  index = sample(1:47, 47 - round(nrow((true_size_data[true_size_data["outlier"] == "no",]))*out_ratio/(1-out_ratio))) 
  out_idx = (true_size_data[true_size_data["outlier"] == "yes",])[index,"id"]
  data = true_size_data[-out_idx,]
  
  
  target = as.vector(data[,35])
  target[target == "yes"] = rep(1, length(target[target == "yes"]))
  target[target == "no"] = rep(0, length(target[target == "no"]))
  set.seed(seed)
  index = as.vector(createDataPartition(target, p = tsplit, list = FALSE))
  train_data <- data[index,-c(34,35)]
  train_labels <- target[index]
  test_data <- rbind(data[-index,-c(34,35)],filter(true_size_data, !id %in% data$id)[,-c(34,35)])
  test_labels <- c(target[-index],rep(1,length(filter(true_size_data, !id %in% data$id)[,35])))
  
  if(hog_exp == TRUE){
    DB_gen(train_data)
    rm(ODM_env, envir = globalenv())
  }
  
  return(list(train_data,train_labels,test_data,test_labels))
}


### ARRYTHMIA ------------------------------------------------------------
load_ARRYTHMIA <- function(tsplit = .2, out_ratio = "02", hog_exp = TRUE, seed = 47385123){
  #' @title Loads the named dataset from the dataset directories. To download the
  #' datasets follow the instructions of the repository. It outputs a list
  #' containing the train set, the test set, the train labels and the test labels
  #' Arguments:
  #' @param tsplit: Split of train-test
  #' @param out_ratio: Ratio of outliers for the train test
  #' @param hog_exp: Selector to check if the experiment is using hidden outliers.
  #' It allows extra computations (generates the DB object)
  #' @param seed: Seed for the split
  
  data = readARFF(glue("databases/Arrhythmia_withoutdupl_norm_{out_ratio}_v10.arff"))
  true_size_data = readARFF("databases/true_size/Arrhythmia_withoutdupl_norm_46.arff")
  
  target = as.vector(data[,261])
  target[target == "yes"] = rep(1, length(target[target == "yes"]))
  target[target == "no"] = rep(0, length(target[target == "no"]))
  set.seed(seed)
  index = as.vector(createDataPartition(target, p = tsplit, list = FALSE))
  train_data <- data[index,-c(260,261)]
  train_labels <- target[index]
  test_data <-  rbind(data[-index,-c(260,261)],filter(true_size_data, !id %in% data$id)[,-c(260,261)])
  test_labels <- c(target[-index],rep(1,length(filter(true_size_data, !id %in% data$id)[,261])))
  
  if(hog_exp == TRUE){
    DB_gen(train_data, true_inliers = FALSE)
    rm(ODM_env, envir = globalenv())
  }
  
  return(list(train_data,train_labels,test_data,test_labels))
}

### PARKINSON ------------------------------------------------------------
load_PARKINSON <- function(tsplit = .2, out_ratio = "05", hog_exp = TRUE, seed = 47385123){
  #' @title Loads the named dataset from the dataset directories. To download the
  #' datasets follow the instructions of the repository. It outputs a list
  #' containing the train set, the test set, the train labels and the test labels
  #' Arguments:
  #' @param tsplit: Split of train-test
  #' @param out_ratio: Ratio of outliers for the train test
  #' @param hog_exp: Selector to check if the experiment is using hidden outliers.
  #' It allows extra computations (generates the DB object)
  #' @param seed: Seed for the split
  
  data = readARFF(glue("databases/Parkinson_withoutdupl_norm_{out_ratio}_v10.arff"))
  true_size_data = readARFF("databases/true_size/Parkinson_withoutdupl_norm_75.arff")
  
  target = as.vector(data[,24])
  target[target == "yes"] = rep(1, length(target[target == "yes"]))
  target[target == "no"] = rep(0, length(target[target == "no"]))
  set.seed(seed)
  index = as.vector(createDataPartition(target, p = tsplit, list = FALSE))
  train_data <- data[index,-c(23,24)]
  train_labels <- target[index]
  test_data <-  rbind(data[-index,-c(23,24)],setNames(filter(true_size_data, !id %in% data$id)[,-c(1,24)],names(data[,-c(23,24)])))
  test_labels <- c(target[-index],rep(1,length(filter(true_size_data, !id %in% data$id)[,24])))
  
  if(hog_exp == TRUE){
    DB_gen(train_data, true_inliers = FALSE)
    rm(ODM_env, envir = globalenv())
  }
  return(list(train_data,train_labels,test_data,test_labels))
}

### CARDIOTOCOGRAPHY ------------------------------------------------------------
load_CARDIOTOCOGRAPHY <- function(tsplit = .2, out_ratio = "02", hog_exp = TRUE, seed = 47385123){
  #' @title Loads the named dataset from the dataset directories. To download the
  #' datasets follow the instructions of the repository. It outputs a list
  #' containing the train set, the test set, the train labels and the test labels
  #' Arguments:
  #' @param tsplit: Split of train-test
  #' @param out_ratio: Ratio of outliers for the train test
  #' @param hog_exp: Selector to check if the experiment is using hidden outliers.
  #' It allows extra computations (generates the DB object)
  #' @param seed: Seed for the split
  
  data = readARFF(glue("databases/Cardiotocography_withoutdupl_norm_{out_ratio}_v10.arff"))
  true_size_data = readARFF("databases/true_size/Cardiotocography_withoutdupl_norm_22.arff")
  
  target = as.vector(data[,23])
  target[target == "yes"] = rep(1, length(target[target == "yes"]))
  target[target == "no"] = rep(0, length(target[target == "no"]))
  set.seed(seed)
  index = as.vector(createDataPartition(target, p = tsplit, list = FALSE))
  train_data <- data[index,-c(22,23)]
  train_labels <- target[index]
  test_data <-  rbind(data[-index,-c(22,23)],setNames(filter(true_size_data, !id %in% data$id)[,-c(22,23)],names(data[,-c(22,23)])))
  test_labels <- c(target[-index],rep(1,length(filter(true_size_data, !id %in% data$id)[,23])))
  
  if(hog_exp == TRUE){
    DB_gen(train_data, true_inliers = FALSE)
    rm(ODM_env, envir = globalenv())
  }
  return(list(train_data,train_labels,test_data,test_labels))
}

### SPAMBASE ------------------------------------------------------------
load_SPAMBASE <- function(tsplit = .2, out_ratio = "02", hog_exp = TRUE, seed = 47385123){
  #' @title Loads the named dataset from the dataset directories. To download the
  #' datasets follow the instructions of the repository. It outputs a list
  #' containing the train set, the test set, the train labels and the test labels
  #' Arguments:
  #' @param tsplit: Split of train-test
  #' @param out_ratio: Ratio of outliers for the train test
  #' @param hog_exp: Selector to check if the experiment is using hidden outliers.
  #' It allows extra computations (generates the DB object)
  #' @param seed: Seed for the split
  
  data = readARFF(glue("databases/SpamBase_withoutdupl_norm_{out_ratio}_v10.arff"))
  true_size_data = readARFF("databases/true_size/SpamBase_withoutdupl_norm_40.arff")
  
  target = as.vector(data[,59])
  target[target == "yes"] = rep(1, length(target[target == "yes"]))
  target[target == "no"] = rep(0, length(target[target == "no"]))
  set.seed(seed)
  index = as.vector(createDataPartition(target, p = tsplit, list = FALSE))
  train_data <- data[index,-c(58,59)]
  train_labels <- target[index]
  test_data <-  rbind(data[-index,-c(58,59)],setNames(filter(true_size_data, !id %in% data$id)[,-c(58,59)],names(data[,-c(58,59)])))
  test_labels <- c(target[-index],rep(1,length(filter(true_size_data, !id %in% data$id)[,59])))
  
  if(hog_exp == TRUE){
    DB_gen(train_data, true_inliers = FALSE)
    rm(ODM_env, envir = globalenv())
  }
  return(list(train_data,train_labels,test_data,test_labels))
}

### IONOSPHERE -----------------------------------------------------------------
load_IONOSPHERE <- function(tsplit = .2, out_ratio = "02", hog_exp = TRUE, seed = 47385123){
  #' @title Loads the named dataset from the dataset directories. To download the
  #' datasets follow the instructions of the repository. It outputs a list
  #' containing the train set, the test set, the train labels and the test labels
  #' Arguments:
  #' @param tsplit: Split of train-test
  #' @param out_ratio: Ratio of outliers for the train test
  #' @param hog_exp: Selector to check if the experiment is using hidden outliers.
  #' It allows extra computations (generates the DB object)
  #' @param seed: Seed for the split
  
  out_ratio = as.integer(out_ratio)/100
  true_size_data = readARFF("databases/true_size/Ionosphere_withoutdupl_norm.arff")
  set.seed(47385123) #So we always get the same downsampling
  index = sample(1:126, 126 - round(nrow((true_size_data[true_size_data["outlier"] == "no",]))*out_ratio/(1-out_ratio)))
  out_idx = (true_size_data[true_size_data["outlier"] == "yes",])[index,"id"]
  data = true_size_data[-out_idx,]
  
  
  target = as.vector(data[,34])
  target[target == "yes"] = rep(1, length(target[target == "yes"]))
  target[target == "no"] = rep(0, length(target[target == "no"]))
  set.seed(seed)
  index = as.vector(createDataPartition(target, p = tsplit, list = FALSE))
  train_data <- data[index,-c(33,34)]
  train_labels <- target[index]
  test_data <- rbind(data[-index,-c(33,34)],filter(true_size_data, !id %in% data$id)[,-c(33,34)])
  test_labels <- c(target[-index],rep(1,length(filter(true_size_data, !id %in% data$id)[,34])))
  
  if(hog_exp == TRUE){
    DB_gen(train_data)
    rm(ODM_env, envir = globalenv())
  }
  
  return(list(train_data,train_labels,test_data,test_labels))
}

### ANNTHYROID -----------------------------------------------------------------

load_ANNTHYROID <- function(tsplit = .2, out_ratio = "02", hog_exp = TRUE, seed = 47385123){
  #' @title Loads the named dataset from the dataset directories. To download the
  #' datasets follow the instructions of the repository. It outputs a list
  #' containing the train set, the test set, the train labels and the test labels
  #' Arguments:
  #' @param tsplit: Split of train-test
  #' @param out_ratio: Ratio of outliers for the train test
  #' @param hog_exp: Selector to check if the experiment is using hidden outliers.
  #' It allows extra computations (generates the DB object)
  #' @param seed: Seed for the split
  
  data = readARFF(glue("databases/Annthyroid_withoutdupl_norm_{out_ratio}_v10.arff"))
  true_size_data = readARFF("databases/true_size/Annthyroid_withoutdupl_norm_07.arff")
  
  target = as.vector(data[,23])
  target[target == "yes"] = rep(1, length(target[target == "yes"]))
  target[target == "no"] = rep(0, length(target[target == "no"]))
  set.seed(seed)
  index = as.vector(createDataPartition(target, p = tsplit, list = FALSE))
  train_data <- data[index,-c(22,23)]
  train_labels <- target[index]
  test_data <-  rbind(data[-index,-c(22,23)],setNames(filter(true_size_data, !id %in% data$id)[,-c(22,23)],names(data[,-c(22,23)])))
  test_labels <- c(target[-index],rep(1,length(filter(true_size_data, !id %in% data$id)[,23])))
  
  if(hog_exp == TRUE){
    DB_gen(train_data, true_inliers = FALSE)
    rm(ODM_env, envir = globalenv())
  }
  return(list(train_data,train_labels,test_data,test_labels))
}

### PAGEBLOCKS -----------------------------------------------------------------
load_PAGEBLOCKS <- function(tsplit = .2, out_ratio = "02", hog_exp = TRUE, seed = 47385123){
  #' @title Loads the named dataset from the dataset directories. To download the
  #' datasets follow the instructions of the repository. It outputs a list
  #' containing the train set, the test set, the train labels and the test labels
  #' Arguments:
  #' @param tsplit: Split of train-test
  #' @param out_ratio: Ratio of outliers for the train test
  #' @param hog_exp: Selector to check if the experiment is using hidden outliers.
  #' It allows extra computations (generates the DB object)
  #' @param seed: Seed for the split
  
  data = readARFF(glue("databases/PageBlocks_withoutdupl_norm_{out_ratio}_v10.arff"))
  true_size_data = readARFF("databases/true_size/PageBlocks_withoutdupl_09.arff")
  
  target = as.vector(data[,12])
  target[target == "yes"] = rep(1, length(target[target == "yes"]))
  target[target == "no"] = rep(0, length(target[target == "no"]))
  set.seed(seed)
  index = as.vector(createDataPartition(target, p = tsplit, list = FALSE))
  train_data <- data[index,-c(11,12)]
  train_labels <- target[index]
  test_data <-  rbind(data[-index,-c(11,12)],setNames(filter(true_size_data, !id %in% data$id)[,-c(11,12)],names(data[,-c(11,12)])))
  test_labels <- c(target[-index],rep(1,length(filter(true_size_data, !id %in% data$id)[,12])))
  
  if(hog_exp == TRUE){
    DB_gen(train_data, true_inliers = FALSE)
    rm(ODM_env, envir = globalenv())
  }
  return(list(train_data,train_labels,test_data,test_labels))
}

### WILT -----------------------------------------------------------------
load_WILT <- function(tsplit = .2, out_ratio = "02", hog_exp = TRUE, seed = 47385123){
  #' @title Loads the named dataset from the dataset directories. To download the
  #' datasets follow the instructions of the repository. It outputs a list
  #' containing the train set, the test set, the train labels and the test labels
  #' Arguments:
  #' @param tsplit: Split of train-test
  #' @param out_ratio: Ratio of outliers for the train test
  #' @param hog_exp: Selector to check if the experiment is using hidden outliers.
  #' It allows extra computations (generates the DB object)
  #' @param seed: Seed for the split
  
  data = readARFF(glue("databases/Wilt_norm_{out_ratio}_v10.arff"))
  true_size_data = readARFF("databases/true_size/Wilt_norm_05.arff")
  
  target = as.vector(data[,7])
  target[target == "yes"] = rep(1, length(target[target == "yes"]))
  target[target == "no"] = rep(0, length(target[target == "no"]))
  set.seed(seed)
  index = as.vector(createDataPartition(target, p = tsplit, list = FALSE))
  train_data <- data[index,-c(6,7)]
  train_labels <- target[index]
  test_data <-  rbind(data[-index,-c(6,7)],setNames(filter(true_size_data, !id %in% data$id)[,-c(6,7)],names(data[,-c(6,7)])))
  test_labels <- c(target[-index],rep(1,length(filter(true_size_data, !id %in% data$id)[,7])))
  
  if(hog_exp == TRUE){
    DB_gen(train_data, true_inliers = FALSE)
    rm(ODM_env, envir = globalenv())
  }
  return(list(train_data,train_labels,test_data,test_labels))
}


