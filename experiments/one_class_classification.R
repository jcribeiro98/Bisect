# One class classification --------------------------------------------
library(miceadds)
source.all("src/HOGen/")
source("experiments/lof_based_generator.R")
library(caret)
library(pROC)
library(ROCR)
library(DMwR)
library(farff)
library(e1071)


## Function preamble --------------------------------------------------------
generate_casting <- function(seed,
                             train_data, train_labels, test_data, test_labels,
                             synth_data){
  #' @title Trains a given two-class classifier with self-supervision
  #' @description Function to handle the one-class classification experiments. 
  #' Given a two-class method, train_data, synth_data and test_data; it trains 
  #' the balanced classifier and then tests on the test data, outputting the 
  #' auc of the training. 
  #' @param seed_set: Seed for the single experiment
  #' @param train_data: Train data to balance
  #' @param train_labels: Labels for the train data
  #' @param test_data: Test data to test with
  #' @param test_labels: Labels for the test data
  #' @param synth_data: Synthetic data generated with an oversampling method
  #' @param synth_labels: Labels for the synthetic data
  
  i = 0
  set.seed(seed)
  train_data_dummy = rbind(train_data, synth_data)
  rownames(train_data_dummy) = NULL #Renumbers the dataframe
  train_labels_dummy = c(rep(0, nrow(train_data)),rep(1, nrow(synth_data)))
      
  set.seed(seed)
  cl <- makePSOCKcluster(5)
  registerDoParallel(cl)
  rf_classifier <- train(train_data_dummy, 
                         factor(train_labels_dummy), method = "rf")
  stopCluster(cl)
  predictions <- predict(rf_classifier, test_data, type = "prob")[,1]
  rocauc_vector <- auc(test_labels, predictions)
  prauc_vector <- MLmetrics::PRAUC(predictions,test_labels)
  
  return(rbind(rocauc_vector,prauc_vector))
}



## Detection Methods -----------------------------------------------------------------

### Bisect ------------------------------------------------------------------

run_bisect <- function(dataset, m_vector = 1, origin = "weighted",
                          seed_vector = c(47385123, 12345, 54321, 010101, 
                                          121212, 19, 33),tsplit = .8, method = "pyod_LOF"){
  #' @title: Experiment launcher for the RF+Bisect one-class classifier
  #' @description  Launch an experiment on the given dataset with multiple
  #' repetition using the given vector of seeds. 
  #'
  #' Arguments:
  #' @param dataset: Chose a dataset given as a loader function
  #' @param tsplit: Split for the train-test
  #' @param method: Adversary to use
  #' @param origin: Origin selection method
  #' @param seed_vector: Vector of seeds for multiple reproductions of the 
  #' experiment
  #' @param m_vector: Number of artificial instances to generate. Defaults to 1
  #' @param meths: Two-class classification method to oversample. Options are:
  #'               -svmLinear2: SVM with linear kernel
  #'               -Any classifier on the "caret" package. For example, 'rf'
  #'               which is random forest.
  
  loader <- match.fun(glue("load_{dataset}"))
  
  for (seed in seed_vector){
    for(m in m_vector){
      data = loader(tsplit = tsplit, seed = seed)
      train_data = data[[1]]; train_labels = data[[2]] 
      test_data = data[[3]]; test_labels = data[[4]]
      hog_model = main_multibisect(gen_points = floor(m*nrow(train_data)),
                                   method = method, seed = seed, 
                                   num_workers = 1, l_val_option = "not_fixed",
                                   type = origin, verb = FALSE, dev_opt = TRUE)
      
      synth_data <- as.data.frame(hog_model$ho_list)
      names(synth_data) <- names(train_data)  
      assign(glue("{m}_{tolower(dataset)}_bisect_results_{seed}"), 
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
     result = cbind(result, get(glue("{m}_{tolower(dataset)}_bisect_results_{seed}")))
   }
  rownames(result) <- c("ROC", "Prec-Rec") 
  assign(glue("{m}_{tolower(dataset)}_bisect_results"),result, envir = globalenv())
  }
}

### Hidden ------------------------------------------------------------------

run_hidden <- function(dataset, m_vector, origin = "weighted",
                       seed_vector = c(47385123, 12345, 54321, 010101, 
                                       121212, 19, 33),tsplit = .8,eps = .1){
  #' @title: Experiment launcher for the RF+Hidden one-class classifier
  #' @description  Launch an experiment on the given dataset with multiple
  #' repetition using the given vector of seeds. 
  #'
  #' Arguments:
  #' @param dataset: Chose a dataset given as a loader function
  #' @param tsplit: Split for the train-test
  #' @param method: Adversary to use
  #' @param eps: Value for epsilon.
  #' @param seed_vector: Vector of seeds for multiple reproductions of the 
  #' experiment
  #' @param m_vector: Number of artificial instances to generate. Defaults to 1
  #' @param meths: Two-class classification method to oversample. Options are:
  #'               -svmLinear2: SVM with linear kernel
  #'               -Any classifier on the "caret" package. For example, 'rf'
  #'               which is random forest.
  
  loader <- match.fun(glue("load_{dataset}"))
  
  for (seed in seed_vector){
    for(m in m_vector){
      data = loader(tsplit = tsplit, seed = seed)
      train_data = data[[1]]; train_labels = data[[2]] 
      test_data = data[[3]]; test_labels = data[[4]]
      set.seed(seed)
      hog_model = main_hidden(gen_points = floor(m*nrow(train_data)),
                                   method = "pyod_LOF",
                                   num_workers = 1, until_gen = TRUE,
                                   eps = eps, dev_opt = TRUE)
      
      synth_data <- as.data.frame(hog_model$ho_list)
      names(synth_data) <- names(train_data)  
      assign(glue("{m}_{tolower(dataset)}_hidden_results_{seed}"), 
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
      result = cbind(result, get(glue("{m}_{tolower(dataset)}_hidden_results_{seed}")))
    }
    rownames(result) <- c("ROC", "Prec-Rec") 
    assign(glue("{m}_{tolower(dataset)}_hidden_results"),result, envir = globalenv())
  }
}

### OCC Baselines ------------------------------------------------------------
load_odm <- function(odm){
  #' @description Auxiliary function to load an specific odm for the baselines.
  
  if(odm == "LOF"){
    odm <- import("pyod.models.lof")
    odm = odm$LOF()
  }
  if(odm == "DeepSVDD"){
    odm <- import("pyod.models.deep_svdd")
    odm <- odm$DeepSVDD()
  }
  if(odm == "GAAL"){
    odm <- import("pyod.models.mo_gaal")
    odm <- odm$MO_GAAL()
  }
  if(odm == "KNN"){
    odm <- import("pyod.models.knn")
    odm <- odm$KNN()
  }
  if(odm == "AnoGAN"){
    odm <- import("pyod.models.anogan")
    odm <- odm$AnoGAN()
  }
  if(odm == "GMM"){
    odm <- import("pyod.models.gmm")
    odm <- odm$GMM()
  }
  if(odm == "OCSVM"){
    odm <- import("pyod.models.ocsvm")
    odm <- odm$OCSVM()
  }
  return(odm)
}

run_unsup <- function(dataset,seed_vector = c(47385123L, 12345L, 54321L, 
                                              010101L, 121212L, 19L, 
                                              33L),
                      odm_name = "LOF",tsplit = .8){
  #' @title: Experiment launcher for the selected one-class classifier among the 
  #' baselines. It loads from pyod.
  #' @description  Launch an experiment on the given dataset with multiple
  #' repetition using the given vector of seeds. 
  #'
  #' Arguments:
  #' @param dataset: Chose a dataset given as a loader function
  #' @param tsplit: Split for the train-test
  #' @param odm_name: Name of the ODM to use
  #' @param seed_vector: Vector of seeds for multiple reproductions of the 
  #' experiment
  #' 
  
  loader <- match.fun(glue("load_{dataset}"))
    for(seed in seed_vector){
      data = loader(tsplit = tsplit,seed = seed)
      train_data = data[[1]]; train_labels = data[[2]]
      test_data = data[[3]]; test_labels = data[[4]]
      odm <- load_odm(odm_name)
      odm$fit(train_data)
      predictions = odm$decision_function(test_data)
      rocauc_vector <- auc(test_labels, predictions)
      prauc_vector <- MLmetrics::PRAUC(as.numeric(predictions),test_labels)
      assign(glue("{tolower(dataset)}_{tolower(odm_name)}_results_{seed}"),
             rbind(rocauc_vector, prauc_vector)) 
      
    }
  result = c()
  for(seed in seed_vector){
    result = cbind(result, get(glue("{tolower(dataset)}_{tolower(odm_name)}_results_{seed}")))
  }
  rownames(result) <- c("ROC", "Prec-Rec") 
  assign(glue("{tolower(dataset)}_{tolower(odm_name)}_results"),result, envir = globalenv())
}

### Hyperbox ------------------------------------------------------------------

run_hyperbox <- function(dataset, m_vector,
                       seed_vector = c(47385123, 12345, 54321, 010101, 
                                       121212, 19, 33),tsplit = .8){
  #' @title: Experiment launcher for the RF+HB one-class classifier
  #' @description  Launch an experiment on the given dataset with multiple
  #' repetition using the given vector of seeds. 
  #'
  #' Arguments:
  #' @param dataset: Chose a dataset given as a loader function
  #' @param tsplit: Split for the train-test
  #' @param seed_vector: Vector of seeds for multiple reproductions of the 
  #' experiment
  #' @param m_vector: Number of artificial instances to generate. Defaults to 1
  
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

### Bisect ------------------------------------------------------------------

run_all_bisect <- function(origin = "weighted", m_vector = c(.25,.5,1),
                           datasets = c("WILT", "PIMA", "STAMPS", "PAGEBLOCKS", 
                                        "HEARTDISEASE", "ANNTHYROID",
                                        "CARDIOTOCOGRAPHY", "PARKINSON",
                                        "IONOSPHERE", "WPBC", "SPAMBASE",
                                        "ARRYTHMIA"),
                           seed_vector = c(47385123L, 12345L, 54321L, 
                                           010101L, 121212L, 19L, 33L), method = "pyod_LOF"){
  for (dataset in datasets){ 
    run_bisect(glue("{toupper(dataset)}"), origin , 
               seed_vector = seed_vector, m_vector = m_vector, method = method) 
    save.image(glue("experiments/one_class_classification/RData/bisect_exp_{Sys.time()}.RData"))
  }
}

### Hidden ------------------------------------------------------------------

run_all_hidden <- function(eps = .1, m_vector = c(.25,.5,1),
                           datasets = c("HEARTDISEASE", 
                                        "HEPATITIS",
                                        "PIMA",
                                        "PARKINSON",
                                        "STAMPS", 
                                        "WPBC",
                                        "IONOSPHERE",
                                        "SPAMBASE",
                                        "CARDIOTOCOGRAPHY", 
                                        "ANNTHYROID", 
                                        "PAGEBLOCKS",
                                        "ARRYTHMIA",
                                        "WILT"), 
                           seed_vector = c(47385123L, 12345L, 54321L, 
                                           010101L, 121212L, 19L, 33L)){
  for (dataset in datasets){ 
    run_hidden(glue("{toupper(dataset)}"), eps = eps, 
               seed_vector = seed_vector, m_vector = m_vector) 
    save.image(glue("experiments/one_class_classification/RData/{eps}_hidden_exp_{Sys.time()}.RData"))
  }
}

### OCC Baselines -----------------------------------------------------------

run_all_unsup <- function(datasets = c("HEARTDISEASE", "HEPATITIS", "PIMA",
                                       "PARKINSON", "STAMPS", "WPBC",
                                       "IONOSPHERE", "SPAMBASE", "CARDIOTOCOGRAPHY", 
                                       "ANNTHYROID", "PAGEBLOCKS", "ARRYTHMIA",
                                       "WILT"), seed_vector = c(47385123L, 12345L, 54321L, 
                                                                010101L, 121212L, 19L, 33L),
                          odm = "LOF"){
  for (dataset in datasets){ 
    run_unsup(glue("{toupper(dataset)}"), 
              seed_vector = seed_vector, odm_name = odm) 
    save.image(glue("experiments/one_class_classification/RData/{tolower(odm)}_exp_{Sys.time()}.RData"))
  }
}

### Hyperbox ------------------------------------------------------------------

run_all_hyperbox <- function(datasets = c("HEARTDISEASE", "HEPATITIS", "PIMA",
                                          "PARKINSON", "STAMPS", "WPBC",
                                          "IONOSPHERE", "SPAMBASE", "CARDIOTOCOGRAPHY", 
                                          "ANNTHYROID", "PAGEBLOCKS", "ARRYTHMIA",
                                          "WILT"), seed_vector = c(47385123L, 12345L, 54321L, 
                                                                   010101L, 121212L, 19L, 33L),
                             m_vector = c(.25,.5,1)){
  for (dataset in datasets){ 
    run_hyperbox(glue("{toupper(dataset)}"), 
                 seed_vector = seed_vector, m_vector = m_vector) 
    save.image(glue("experiments/one_class_classification/RData/hyperbox_exp_{Sys.time()}.RData"))
  }
}


## Datasets loaders -----------------------------------------------------------

load_ARRYTHMIA <- function(tsplit = .2, 
                           hog_exp = TRUE, seed = 47385123){
  #' @title Loads the named dataset from the dataset directories. To download the
  #' datasets follow the instructions of the repository. It outputs a list
  #' containing the train set, the test set, the train labels and the test labels
  #' Arguments:
  #' @param tsplit: Split of train-test
  #' @param hog_exp: Selector to check if the experiment is using hidden outliers.
  #' It allows extra computations (generates the DB object)
  #' @param seed: Seed for the split

  true_size_data = readARFF("databases/true_size/Arrhythmia_withoutdupl_norm_46.arff")
  inliers = true_size_data %>% filter(outlier == "no")
  outliers = true_size_data %>% filter(outlier == "yes")
  set.seed(seed)
  
  index = as.vector(createDataPartition(inliers[,"outlier"], p = tsplit, list = FALSE))
  train_data <- inliers[index,-c(260,261)]
  train_labels <- rep(0,length(index))
  test_data <-  rbind(inliers[-index,-c(260,261)],outliers[,-c(260,261)])
  test_labels <- c(rep(0,nrow(inliers[-index,])), rep(1,nrow(outliers)))
  
  if(hog_exp == TRUE){
    DB_gen(train_data, true_inliers = FALSE)
    
  }
  return(list(train_data,train_labels,test_data,test_labels))
}

load_HEARTDISEASE <- function(tsplit = .2, 
                           hog_exp = TRUE, seed = 47385123){
  #' @title Loads the named dataset from the dataset directories. To download the
  #' datasets follow the instructions of the repository. It outputs a list
  #' containing the train set, the test set, the train labels and the test labels
  #' Arguments:
  #' @param tsplit: Split of train-test
  #' @param hog_exp: Selector to check if the experiment is using hidden outliers.
  #' It allows extra computations (generates the DB object)
  #' @param seed: Seed for the split
  
  true_size_data = readARFF("databases/true_size/HeartDisease_withoutdupl_norm_44.arff")
  inliers = (true_size_data %>% filter(outlier == "no"))
  outliers = (true_size_data %>% filter(outlier == "yes"))
  set.seed(seed)
  
  index = as.vector(createDataPartition(inliers[,"outlier"], p = tsplit, list = FALSE))
  train_data <- inliers[index,-c(1,15)]
  train_labels <- rep(0,length(index))
  test_data <-  rbind(inliers[-index,-c(1,15)],outliers[,-c(1,15)])
  test_labels <- c(rep(0,nrow(inliers[-index,])), rep(1,nrow(outliers)))
  
  if(hog_exp == TRUE){
    DB_gen(train_data, true_inliers = FALSE)
    
  }
  return(list(train_data,train_labels,test_data,test_labels))
}

load_HEPATITIS <- function(tsplit = .2, 
                              hog_exp = TRUE, seed = 47385123){
  #' @title Loads the named dataset from the dataset directories. To download the
  #' datasets follow the instructions of the repository. It outputs a list
  #' containing the train set, the test set, the train labels and the test labels
  #' Arguments:
  #' @param tsplit: Split of train-test
  #' @param hog_exp: Selector to check if the experiment is using hidden outliers.
  #' It allows extra computations (generates the DB object)
  #' @param seed: Seed for the split
  
  true_size_data = readARFF("databases/true_size/Hepatitis_withoutdupl_norm_16.arff")
  inliers = (true_size_data %>% filter(outlier == "no"))
  outliers = (true_size_data %>% filter(outlier == "yes"))
  set.seed(seed)
  
  index = as.vector(createDataPartition(inliers[,"outlier"], p = tsplit, list = FALSE))
  train_data <- inliers[index,-c(20,21)]
  train_labels <- rep(0,length(index))
  test_data <-  rbind(inliers[-index,-c(20,21)],outliers[,-c(20,21)])
  test_labels <- c(rep(0,nrow(inliers[-index,])), rep(1,nrow(outliers)))
  
  if(hog_exp == TRUE){
    DB_gen(train_data, true_inliers = FALSE)
    
  }
  return(list(train_data,train_labels,test_data,test_labels))
}

load_PIMA <- function(tsplit = .2, 
                      hog_exp = TRUE, seed = 47385123){
  #' @title Loads the named dataset from the dataset directories. To download the
  #' datasets follow the instructions of the repository. It outputs a list
  #' containing the train set, the test set, the train labels and the test labels
  #' Arguments:
  #' @param tsplit: Split of train-test
  #' @param hog_exp: Selector to check if the experiment is using hidden outliers.
  #' It allows extra computations (generates the DB object)
  #' @param seed: Seed for the split
  
  true_size_data = readARFF("databases/true_size/Pima_withoutdupl_norm_35.arff")
  inliers = (true_size_data %>% filter(outlier == "no"))
  outliers = (true_size_data %>% filter(outlier == "yes"))
  set.seed(seed)
  
  index = as.vector(createDataPartition(inliers[,"outlier"], p = tsplit, list = FALSE))
  train_data <- inliers[index,-c(1,10)]
  train_labels <- rep(0,length(index))
  test_data <-  rbind(inliers[-index,-c(1,10)],outliers[,-c(1,10)])
  test_labels <- c(rep(0,nrow(inliers[-index,])), rep(1,nrow(outliers)))
  
  if(hog_exp == TRUE){
    DB_gen(train_data, true_inliers = FALSE)
    
  }
  return(list(train_data,train_labels,test_data,test_labels))
}

load_PARKINSON <- function(tsplit = .2, 
                      hog_exp = TRUE, seed = 47385123){
  #' @title Loads the named dataset from the dataset directories. To download the
  #' datasets follow the instructions of the repository. It outputs a list
  #' containing the train set, the test set, the train labels and the test labels
  #' Arguments:
  #' @param tsplit: Split of train-test
  #' @param hog_exp: Selector to check if the experiment is using hidden outliers.
  #' It allows extra computations (generates the DB object)
  #' @param seed: Seed for the split
  
  true_size_data = readARFF("databases/true_size/Parkinson_withoutdupl_norm_75.arff")
  inliers = (true_size_data %>% filter(outlier == "no"))
  outliers = (true_size_data %>% filter(outlier == "yes"))
  set.seed(seed)
  
  index = as.vector(createDataPartition(inliers[,"outlier"], p = tsplit, list = FALSE))
  train_data <- inliers[index,-c(1,24)]
  train_labels <- rep(0,length(index))
  test_data <-  rbind(inliers[-index,-c(1,24)],outliers[,-c(1,24)])
  test_labels <- c(rep(0,nrow(inliers[-index,])), rep(1,nrow(outliers)))
  
  if(hog_exp == TRUE){
    DB_gen(train_data, true_inliers = FALSE)
    
  }
  return(list(train_data,train_labels,test_data,test_labels))
}

load_STAMPS <- function(tsplit = .2, 
                           hog_exp = TRUE, seed = 47385123){
  #' @title Loads the named dataset from the dataset directories. To download the
  #' datasets follow the instructions of the repository. It outputs a list
  #' containing the train set, the test set, the train labels and the test labels
  #' Arguments:
  #' @param tsplit: Split of train-test
  #' @param hog_exp: Selector to check if the experiment is using hidden outliers.
  #' It allows extra computations (generates the DB object)
  #' @param seed: Seed for the split
  
  true_size_data = readARFF("databases/true_size/Stamps_withoutdupl_norm_09.arff")
  inliers = (true_size_data %>% filter(outlier == "no"))
  outliers = (true_size_data %>% filter(outlier == "yes"))
  set.seed(seed)
  
  index = as.vector(createDataPartition(inliers[,"outlier"], p = tsplit, list = FALSE))
  train_data <- inliers[index,-c(10,11)]
  train_labels <- rep(0,length(index))
  test_data <-  rbind(inliers[-index,-c(10,11)],outliers[,-c(10,11)])
  test_labels <- c(rep(0,nrow(inliers[-index,])), rep(1,nrow(outliers)))
  
  if(hog_exp == TRUE){
    DB_gen(train_data, true_inliers = FALSE)
    
  }
  return(list(train_data,train_labels,test_data,test_labels))
}

load_WPBC <- function(tsplit = .2, 
                        hog_exp = TRUE, seed = 47385123){
  #' @title Loads the named dataset from the dataset directories. To download the
  #' datasets follow the instructions of the repository. It outputs a list
  #' containing the train set, the test set, the train labels and the test labels
  #' Arguments:
  #' @param tsplit: Split of train-test
  #' @param hog_exp: Selector to check if the experiment is using hidden outliers.
  #' It allows extra computations (generates the DB object)
  #' @param seed: Seed for the split
  
  true_size_data = readARFF("databases/WPBC_withoutdupl_norm.arff")
  inliers = (true_size_data %>% filter(outlier == "no"))
  outliers = (true_size_data %>% filter(outlier == "yes"))
  set.seed(seed)
  
  index = as.vector(createDataPartition(inliers[,"outlier"], p = tsplit, list = FALSE))
  train_data <- inliers[index,-c(34,35)]
  train_labels <- rep(0,length(index))
  test_data <-  rbind(inliers[-index,-c(34,35)],outliers[,-c(34,35)])
  test_labels <- c(rep(0,nrow(inliers[-index,])), rep(1,nrow(outliers)))
  
  if(hog_exp == TRUE){
    DB_gen(train_data, true_inliers = FALSE)
    
  }
  return(list(train_data,train_labels,test_data,test_labels))
}

load_IONOSPHERE <- function(tsplit = .2, 
                        hog_exp = TRUE, seed = 47385123){
  #' @title Loads the named dataset from the dataset directories. To download the
  #' datasets follow the instructions of the repository. It outputs a list
  #' containing the train set, the test set, the train labels and the test labels
  #' Arguments:
  #' @param tsplit: Split of train-test
  #' @param hog_exp: Selector to check if the experiment is using hidden outliers.
  #' It allows extra computations (generates the DB object)
  #' @param seed: Seed for the split
  
  true_size_data = readARFF("databases/true_size/Ionosphere_withoutdupl_norm.arff")
  inliers = (true_size_data %>% filter(outlier == "no"))
  outliers = (true_size_data %>% filter(outlier == "yes"))
  set.seed(seed)
  
  index = as.vector(createDataPartition(inliers[,"outlier"], p = tsplit, list = FALSE))
  train_data <- inliers[index,-c(33,34)]
  train_labels <- rep(0,length(index))
  test_data <-  rbind(inliers[-index,-c(33,34)],outliers[,-c(33,34)])
  test_labels <- c(rep(0,nrow(inliers[-index,])), rep(1,nrow(outliers)))
  
  if(hog_exp == TRUE){
    DB_gen(train_data, true_inliers = FALSE)
    
  }
  return(list(train_data,train_labels,test_data,test_labels))
}

load_SPAMBASE <- function(tsplit = .2, 
                        hog_exp = TRUE, seed = 47385123){
  #' @title Loads the named dataset from the dataset directories. To download the
  #' datasets follow the instructions of the repository. It outputs a list
  #' containing the train set, the test set, the train labels and the test labels
  #' Arguments:
  #' @param tsplit: Split of train-test
  #' @param hog_exp: Selector to check if the experiment is using hidden outliers.
  #' It allows extra computations (generates the DB object)
  #' @param seed: Seed for the split
  
  true_size_data = readARFF("databases/true_size/SpamBase_withoutdupl_norm_40.arff")
  inliers = (true_size_data %>% filter(outlier == "no"))
  outliers = (true_size_data %>% filter(outlier == "yes"))
  set.seed(seed)
  
  index = as.vector(createDataPartition(inliers[,"outlier"], p = tsplit, list = FALSE))
  train_data <- inliers[index,-c(58,59)]
  train_labels <- rep(0,length(index))
  test_data <-  rbind(inliers[-index,-c(58,59)],outliers[,-c(58,59)])
  test_labels <- c(rep(0,nrow(inliers[-index,])), rep(1,nrow(outliers)))
  
  if(hog_exp == TRUE){
    DB_gen(train_data, true_inliers = FALSE)
    
  }
  return(list(train_data,train_labels,test_data,test_labels))
}

load_CARDIOTOCOGRAPHY <- function(tsplit = .2, 
                        hog_exp = TRUE, seed = 47385123){
  #' @title Loads the named dataset from the dataset directories. To download the
  #' datasets follow the instructions of the repository. It outputs a list
  #' containing the train set, the test set, the train labels and the test labels
  #' Arguments:
  #' @param tsplit: Split of train-test
  #' @param hog_exp: Selector to check if the experiment is using hidden outliers.
  #' It allows extra computations (generates the DB object)
  #' @param seed: Seed for the split
  
  true_size_data = readARFF("databases/true_size/Cardiotocography_withoutdupl_norm_22.arff")
  inliers = (true_size_data %>% filter(outlier == "no"))
  outliers = (true_size_data %>% filter(outlier == "yes"))
  set.seed(seed)
  
  index = as.vector(createDataPartition(inliers[,"outlier"], p = tsplit, list = FALSE))
  train_data <- inliers[index,-c(22,23)]
  train_labels <- rep(0,length(index))
  test_data <-  rbind(inliers[-index,-c(22,23)],outliers[,-c(22,23)])
  test_labels <- c(rep(0,nrow(inliers[-index,])), rep(1,nrow(outliers)))
  
  if(hog_exp == TRUE){
    DB_gen(train_data, true_inliers = FALSE)
    
  }
  return(list(train_data,train_labels,test_data,test_labels))
}

load_ANNTHYROID <- function(tsplit = .2, 
                        hog_exp = TRUE, seed = 47385123){
  #' @title Loads the named dataset from the dataset directories. To download the
  #' datasets follow the instructions of the repository. It outputs a list
  #' containing the train set, the test set, the train labels and the test labels
  #' Arguments:
  #' @param tsplit: Split of train-test
  #' @param hog_exp: Selector to check if the experiment is using hidden outliers.
  #' It allows extra computations (generates the DB object)
  #' @param seed: Seed for the split
  
  true_size_data = readARFF("databases/true_size/Annthyroid_withoutdupl_norm_07.arff")
  inliers = (true_size_data %>% filter(outlier == "no"))
  outliers = (true_size_data %>% filter(outlier == "yes"))
  set.seed(seed)
  
  index = as.vector(createDataPartition(inliers[,"outlier"], p = tsplit, list = FALSE))
  train_data <- inliers[index,-c(22,23)]
  train_labels <- rep(0,length(index))
  test_data <-  rbind(inliers[-index,-c(22,23)],outliers[,-c(22,23)])
  test_labels <- c(rep(0,nrow(inliers[-index,])), rep(1,nrow(outliers)))
  
  if(hog_exp == TRUE){
    DB_gen(train_data, true_inliers = FALSE)
    
  }
  return(list(train_data,train_labels,test_data,test_labels))
}

load_PAGEBLOCKS <- function(tsplit = .2, 
                        hog_exp = TRUE, seed = 47385123){
  #' @title Loads the named dataset from the dataset directories. To download the
  #' datasets follow the instructions of the repository. It outputs a list
  #' containing the train set, the test set, the train labels and the test labels
  #' Arguments:
  #' @param tsplit: Split of train-test
  #' @param hog_exp: Selector to check if the experiment is using hidden outliers.
  #' It allows extra computations (generates the DB object)
  #' @param seed: Seed for the split
  
  true_size_data = readARFF("databases/true_size/PageBlocks_withoutdupl_09.arff")
  inliers = (true_size_data %>% filter(outlier == "no"))
  outliers = (true_size_data %>% filter(outlier == "yes"))
  set.seed(seed)
  
  index = as.vector(createDataPartition(inliers[,"outlier"], p = tsplit, list = FALSE))
  train_data <- inliers[index,-c(11,12)]
  train_labels <- rep(0,length(index))
  test_data <-  rbind(inliers[-index,-c(11,12)],outliers[,-c(11,12)])
  test_labels <- c(rep(0,nrow(inliers[-index,])), rep(1,nrow(outliers)))
  
  if(hog_exp == TRUE){
    DB_gen(train_data, true_inliers = FALSE)
    
  }
  return(list(train_data,train_labels,test_data,test_labels))
}

load_WILT <- function(tsplit = .2, 
                        hog_exp = TRUE, seed = 47385123){
  #' @title Loads the named dataset from the dataset directories. To download the
  #' datasets follow the instructions of the repository. It outputs a list
  #' containing the train set, the test set, the train labels and the test labels
  #' Arguments:
  #' @param tsplit: Split of train-test
  #' @param hog_exp: Selector to check if the experiment is using hidden outliers.
  #' It allows extra computations (generates the DB object)
  #' @param seed: Seed for the split
  
  true_size_data = readARFF("databases/true_size/Wilt_norm_05.arff")
  inliers = (true_size_data %>% filter(outlier == "no"))
  outliers = (true_size_data %>% filter(outlier == "yes"))
  set.seed(seed)
  
  index = as.vector(createDataPartition(inliers[,"outlier"], p = tsplit, list = FALSE))
  train_data <- inliers[index,-c(6,7)]
  train_labels <- rep(0,length(index))
  test_data <-  rbind(inliers[-index,-c(6,7)],outliers[,-c(6,7)])
  test_labels <- c(rep(0,nrow(inliers[-index,])), rep(1,nrow(outliers)))
  
  if(hog_exp == TRUE){
    DB_gen(train_data, true_inliers = FALSE)
    
  }
  return(list(train_data,train_labels,test_data,test_labels))
}

