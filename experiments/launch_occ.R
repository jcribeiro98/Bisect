source("experiments/one_class_classification.R")

datasets = c("WILT", "PIMA", "STAMPS", "PAGEBLOCKS", "HEARTDISEASE", "ANNTHYROID",
             "CARDIOTOCOGRAPHY", "PARKINSON", "IONOSPHERE", "WPBC", "SPAMBASE",
             "ARRYTHMIA")
run_all_bisect()
run_all_hidden(datasets = datasets[-4], eps = .1)
run_all_bisect(method = "KNN")
run_all_unsup(odm = "LOF")
run_all_unsup(odm = "KNN")
run_all_unsup(odm = "OCSVM")
run_all_unsup(odm = "DeepSVDD")
run_all_unsup(odm = "AnoGAN")
run_all_unsup(odm = "GAAL", datasets = datasets[-8], seed_vector = c(47385123L, 12345L, 
                                                         54321L, 010101L, 
                                                         121212L, 19L,33L))
run_all_hyperbox()

### WRITING RESULTS ---------------------------------------------------------
datasets = c("WILT", "PIMA", "STAMPS", "PAGEBLOCKS", "HEARTDISEASE", "ANNTHYROID",
             "CARDIOTOCOGRAPHY", "PARKINSON", "IONOSPHERE", "WPBC", "SPAMBASE",
             "ARRYTHMIA")

get_res_bisect = function(datasets, metric = "ROC"){
  total_res = c()
  for (m in c(0.25,.5,1)){
    res = c()
    for (dataset in datasets){
      res = rbind(res, median(get(glue("{m}_{tolower(dataset)}_bisect_results"))[metric,]))
    }
    rownames(res) = datasets
    total_res = cbind(total_res, res)
  }
  colnames(total_res) = c(0.25,.5,1)
  return(total_res)
}
  

get_res_unsup = function(datasets,odm, metric = "ROC"){
    res = c()
    for (dataset in datasets){
      res = rbind(res, median(get(glue("{tolower(dataset)}_{tolower(odm)}_results"))[metric,]))
    }
    rownames(res) = datasets
  return(res)
}  

get_res_hidden = function(datasets, metric = "ROC"){
  total_res = c()
  for (m in c(1)){
    res = c()
    for (dataset in datasets){
      res = rbind(res, median(get(glue("{m}_{tolower(dataset)}_hidden_results"))[metric,]))
    }
    rownames(res) = datasets
    total_res = cbind(total_res, res)
  }
  colnames(total_res) = c(1)
  return(total_res)
}

get_res_hyperbox = function(datasets, metric = "ROC"){
  total_res = c()
  for (m in c(0.25,0.5,1)){
    res = c()
    for (dataset in datasets){
      res = rbind(res, median(get(glue("{m}_{tolower(dataset)}_hyperbox_results"))[metric,]))
    }
    rownames(res) = datasets
    total_res = cbind(total_res, res)
  }
  colnames(total_res) = c(0.25,.5,1)
  return(total_res)
}

### WILCOX TEST -------------------------------------------------------------
#Running the wilcox test between the top 2 (2nd vs 1st).
wilcox.test(`1_wilt_bisect_results`["ROC",], `1_wilt_hidden_results`["ROC",])$p.value 
wilcox.test(`1_pima_hyperbox_results`["ROC",],pima_knn_results["ROC",])$p.value 
wilcox.test(`1_stamps_bisect_results`["ROC",], stamps_knn_results["ROC",])$p.value  
wilcox.test(`1_pageblocks_bisect_results`["ROC",],pageblocks_lof_results["ROC",])$p.value  
wilcox.test(`1_heartdisease_hidden_results`["ROC",], `1_heartdisease_hyperbox_results`["ROC",])$p.value  
wilcox.test(`1_annthyroid_hidden_results`["ROC",], annthyroid_deepsvdd_results["ROC",])$p.value  
wilcox.test(`1_cardiotocography_bisect_results`["ROC",], `1_cardiotocography_hidden_results`["ROC",])$p.value  
wilcox.test(parkinson_knn_results["ROC",], parkinson_deepsvdd_results["ROC",])$p.value  
wilcox.test(ionosphere_knn_results["ROC",], ionosphere_deepsvdd_results["ROC",])$p.value  
wilcox.test(wpbc_gaal_results["ROC",], wpbc_knn_results["ROC",])$p.value  
wilcox.test(`1_spambase_bisect_results`["ROC",], `1_spambase_hidden_results`["ROC",])$p.value  
wilcox.test(`1_spambase_hyperbox_results`["ROC",], `1_spambase_hidden_results`["ROC",])$p.value  
wilcox.test(`1_arrythmia_bisect_results`["ROC",], `1_arrythmia_hidden_results`["ROC",])$p.value  


wilcox.test(`1_wilt_bisect_results`["ROC",], wilt_knn_results["ROC",])$p.value
wilcox.test(`1_arrythmia_bisect_results`["ROC",], arrythmia_knn_results["ROC",])$p.value


baseline_comp <- function(datasets = c("WILT", "PIMA", "STAMPS", "PAGEBLOCKS", 
                                      "HEARTDISEASE", "ANNTHYROID",
                                      "CARDIOTOCOGRAPHY", "PARKINSON", 
                                      "IONOSPHERE", "WPBC", "SPAMBASE","ARRYTHMIA"),
                          hog_method="Bisect",baseline="LOF",side = "great"){
  i = 1
  p_val_matrix = matrix(nrow = length(datasets), ncol = 1)
  row.names(p_val_matrix) = datasets
  for (dataset in datasets){
    p_val_matrix[i,] = wilcox.test(get(glue("1_{tolower(dataset)}_{tolower(hog_method)}_results"))["ROC",],get(glue("{tolower(dataset)}_{tolower(baseline)}_results"))["ROC",], alternative = side)$p.value
    
    i = i + 1
  }
  return(p_val_matrix)
}



