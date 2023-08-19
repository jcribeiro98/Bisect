source("experiments/new_balancing_experiments.R")
datasets = c("WILT", "PIMA", "STAMPS", "PAGEBLOCKS" , "HEARTDISEASE", "ANNTHYROID",
             "CARDIOTOCOGRAPHY","PARKINSON","IONOSPHERE", "WPBC", "SPAMBASE",
             "ARRYTHMIA")
### Bisect
run_all_bisect(datasets)

### Hidden
run_all_hidden(datasets = "PAGEBLOCKS", eps = 0.5)

### SMOTE
run_all_synthetic_oversamp(K = 5 , 
                           datasets = c("SPAMBASE","CARDIOTOCOGRAPHY", 
                                        "ANNTHYROID", "PAGEBLOCKS", "WILT"))

### HYPERBOX
run_all_hyperbox(datasets[-8])
### Nothing --------------------------------------------------------------

seed_vector = c(47385123L, 12345L, 54321L, 010101L, 121212L, 19L, 33L)
parkinson_nothing_results <- nothing_exp("PARKINSON", seed_vector = seed_vector)
heartdisease_nothing_results <- nothing_exp("HEARTDISEASE", seed_vector = seed_vector)
pima_nothing_results <- nothing_exp("PIMA", seed_vector = seed_vector)
stamps_nothing_results <- nothing_exp("STAMPS", seed_vector = seed_vector)
WPBC_nothing_results <- nothing_exp("WPBC", seed_vector = seed_vector)
arrythmia_nothing_results <- nothing_exp("ARRYTHMIA", seed_vector = seed_vector)
cardiotocography_bisect_results <-nothing_exp("CARDIOTOCOGRAPHY", seed_vector = seed_vector)
spambase_nothing_results <- nothing_exp("SPAMBASE", seed_vector = seed_vector)
ionosphere_nothing_results <- nothing_exp("IONOSPHERE", seed_vector = seed_vector)
annthyroid_nothing_results <- nothing_exp("ANNTHYROID", seed_vector = seed_vector)
pageblocks_nothing_results <- nothing_exp("PAGEBLOCKS", seed_vector = seed_vector)
wilt_nothing_results <- nothing_exp("WILT", seed_vector = seed_vector)



### Writing data ---------------------------------------------------------
datasets = c("WILT", "PIMA", "STAMPS", "PAGEBLOCKS" , "HEARTDISEASE", "ANNTHYROID",
             "CARDIOTOCOGRAPHY","IONOSPHERE", "WPBC", "SPAMBASE",
             "ARRYTHMIA")
ho_method = "weighted"
class_method = "rf"
for(dataset in datasets){
  assign("x", get(glue("{tolower(dataset)}_bisect_results"))[glue("{ho_method}_{class_method}"),])
  assign("y", get(glue("{tolower(dataset)}_nothing_results"))[glue("{class_method}"),])
  print(glue("{dataset} {wilcox.test(x,y)$p.value} {(wilcox.test(x,y)$p.value < .05)}"))
}
for(dataset in datasets){
  print(glue('{dataset}: {round(median(get(glue("1_{tolower(dataset)}_hidden_results"))["ROC",]), 3)}'))
}
for(dataset in datasets){
  print(glue('{dataset}: {round(median(get(glue("{tolower(dataset)}_nothing_results"))["rf",]), 3)}'))
}
datasets =c("ARRYTHMIA",
            "HEARTDISEASE", 
            "HEPATITIS",
            "PIMA",
            "STAMPS", 
            "WPBC",
            "IONOSPHERE",
            "SPAMBASE",
            "CARDIOTOCOGRAPHY", 
            "ANNTHYROID", 
            "PAGEBLOCKS",
            "WILT")
ho_method = "eps0.1"
for(dataset in datasets){
  assign("x", get(glue("{tolower(dataset)}_hidden_results"))[glue("{class_method}_{ho_method}"),])
  assign("y", get(glue("{tolower(dataset)}_nothing_results"))[glue("{class_method}"),])
  print(glue("{dataset} {wilcox.test(x,y)$p.value} {(wilcox.test(x,y)$p.value < .05)}"))
}
for(dataset in datasets){
  assign("x", get(glue("1_{tolower(dataset)}_hyperbox_results"))["ROC",])
  assign("y", get(glue("{tolower(dataset)}_nothing_results"))[glue("rf"),])
  print(glue("{dataset} {wilcox.test(x,y)$p.value} {(wilcox.test(x,y)$p.value < .05)}"))
}
for(dataset in datasets){
  assign("x", get(glue("{tolower(dataset)}_adasyn_results"))["rf",])
  assign("y", get(glue("{tolower(dataset)}_nothing_results"))[glue("rf"),])
  print(glue("{dataset} {wilcox.test(x,y)$p.value} {(wilcox.test(x,y)$p.value < .05)}"))
}
for(dataset in datasets){
  assign("x", get(glue("{tolower(dataset)}_cGAN_results"))["rf",])
  assign("y", get(glue("{tolower(dataset)}_nothing_results"))[glue("rf"),])
  print(glue("{dataset} {wilcox.test(x,y)$p.value} {(wilcox.test(x,y)$p.value < .05)}"))
}


