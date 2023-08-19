
# New_exec_time -----------------------------------------------------------
library(miceadds)
source.all("src/HOGen/")
library(farff)


## Functions ---------------------------------------------------------------

load_dataset_with_dim <- function(path="databases/Arrhythmia_withoutdupl_norm_02_v10.arff",seed, size){
  #' @title Load a dataset with a specific number of features
  #' @description Loads a dataset from path and downsample the number of features
  #' Arguments:
  #' @param seed: Seed for the downsampling
  #' @param size: Number of features after downsampling
  
  set.seed(seed)
  data = readARFF(path)
  train_data = select(data,-c(id,outlier))
  
  index = sample(1:ncol(train_data), size)
  train_data = train_data[,index]
  DB_gen(train_data, true_inliers = FALSE)
}

time_exper <- function(path = "databases/Arrhythmia_withoutdupl_norm_02_v10.arff",
                       seed_vector = c(47385123, 12345,54321, 010101, 33), 
                       size_vector,eps = .1, ...){
  #' @title Time experiment with real data
  #' @description Launch the time experiment on the specified dataset (with path)
  #' and for specified feature sizes
  #' Arguments:
  #' @param path: Path to dataset
  #' @param seed_vector: 
  
  time_matrix_bisect = matrix(0, ncol = length(seed_vector), 
                       nrow = length(size_vector))
  time_matrix_hidden = time_matrix_bisect
  i = 1
  for (size in size_vector){
    j = 1
    for(seed in seed_vector){
      load_dataset_with_dim(path,seed,size)
      #hog_model = main_multibisect(gen_points = 500, method = "pyod_LOF", 
      #                 seed = seed, num_workers = 1, type = "weighted", dev_opt = F)
      #time_matrix_bisect[i, j] = as.numeric(gsub("[secelapsed ]","", hog_model$exec_time))
      hog_model = main_hidden(gen_points = 500, eps = eps, method = "pyod_LOF", num_workers = 1, until_gen = TRUE)
      time_matrix_hidden[i, j] = as.numeric(gsub("[secelapsed ]","", hog_model$exec_time))
      j = j + 1
    }
  i = i + 1
  }
  res = cbind(time_matrix_bisect,time_matrix_hidden)
  colnames(res) = c(rep("Bisect",length(seed_vector)), rep("Hidden",length(seed_vector)))
  rownames(res) = size_vector
  
  return(res)
}


time_ionosphere <- time_exper(path = "databases/true_size/Ionosphere_withoutdupl_norm.arff", 
                              size_vector = c(10,15,20,25,30), eps = .75)
save.image("experiments/new_exec_time_exp/exec_times_ionosphere_weighted_eps075.RData")
time_wpbc <- time_exper(path = "databases/WPBC_withoutdupl_norm.arff", size_vector = c(10,15,20,25,30), eps = .5)
save.image("experiments/new_exec_time_exp/exec_times_wpbc_weightedeps05.RData")
time_spambase <- time_exper(path = "databases/true_size/SpamBase_withoutdupl_norm_40.arff", size_vector = c(10,20,30,40,50), eps = .75)
save.image("experiments/new_exec_time_exp/exec_time_spambase_weightedeps075.RData")
time_arrythmia <- time_exper(size_vector = c(50,100,150,200,250))
save.image("experiments/new_exec_time_exp/exec_time_arrythmia_weighted.RData")

## Write the CSVs: ---------------------------------------------------------------

plotting_table = data.frame(matrix(0,5*5*2,3))
colnames(plotting_table) = c("Time", "Method", "Feature Size")
plotting_table[,'Method'] = c(rep("Bisect",5*5),rep("Hidden",5*5))
plotting_table[,'Time'] = c(time_arrythmia[,1],time_arrythmia[,2],time_arrythmia[,3],
                            time_arrythmia[,4],time_arrythmia[,5],time_arrythmia[,6],
                            time_arrythmia[,7],time_arrythmia[,8],time_arrythmia[,9],
                            time_arrythmia[,10])
plotting_table[,"Feature Size"] = rep(c(50,100,150,200,250),10)
write.table(plotting_table, "experiments/new_exec_time_exp/exec_time_arryhtmia.csv")

plotting_table[,'Time'] = c(time_wpbc[,1],time_wpbc[,2],time_wpbc[,3],
                            time_wpbc[,4],time_wpbc[,5],time_wpbc[,6],
                            time_wpbc[,7],time_wpbc[,8],time_wpbc[,9],
                            time_wpbc[,10])
plotting_table[,"Feature Size"] = rep(c(10,15,20,25,30),10)
write.table(plotting_table, "experiments/new_exec_time_exp/exec_time_wpbc.csv")

plotting_table[,'Time'] = c(time_spambase[,1],time_spambase[,2],time_spambase[,3],
                            time_spambase[,4],time_spambase[,5],time_spambase[,6],
                            time_spambase[,7],time_spambase[,8],time_spambase[,9],
                            time_spambase[,10])
plotting_table[,"Feature Size"] = rep(c(10,20,30,40,50),10)
write.table(plotting_table, "experiments/new_exec_time_exp/exec_time_spambase.csv")

plotting_table[,'Time'] = c(time_ionosphere[,1],time_ionosphere[,2],time_ionosphere[,3],
                            time_ionosphere[,4],time_ionosphere[,5],time_ionosphere[,6],
                            time_ionosphere[,7],time_ionosphere[,8],time_ionosphere[,9],
                            time_ionosphere[,10])
plotting_table[,"Feature Size"] = rep(c(10,15,20,25,30),10)
write.table(plotting_table, "experiments/new_exec_time_exp/exec_time_ionosphere.csv")

plotting_table = data.frame(matrix(0,5*5*2*4,4))
colnames(plotting_table) = c("Time", "Method", "Feature Size","Dataset")
plotting_table[,'Method'] = rep(c(rep("Bisect",5*5),rep("Hidden",5*5)),4)
plotting_table[,'Time'] = c(c(time_arrythmia[,1],time_arrythmia[,2],time_arrythmia[,3],
                            time_arrythmia[,4],time_arrythmia[,5],time_arrythmia[,6],
                            time_arrythmia[,7],time_arrythmia[,8],time_arrythmia[,9],
                            time_arrythmia[,10]),
                            c(time_wpbc[,1],time_wpbc[,2],time_wpbc[,3],
                              time_wpbc[,4],time_wpbc[,5],time_wpbc[,6],
                              time_wpbc[,7],time_wpbc[,8],time_wpbc[,9],
                              time_wpbc[,10]),
                            c(time_spambase[,1],time_spambase[,2],time_spambase[,3],
                              time_spambase[,4],time_spambase[,5],time_spambase[,6],
                              time_spambase[,7],time_spambase[,8],time_spambase[,9],
                              time_spambase[,10]),
                            c(time_ionosphere[,1],time_ionosphere[,2],time_ionosphere[,3],
                              time_ionosphere[,4],time_ionosphere[,5],time_ionosphere[,6],
                              time_ionosphere[,7],time_ionosphere[,8],time_ionosphere[,9],
                              time_ionosphere[,10]))
plotting_table[,'Dataset'] = c(rep("Arrythmia",5*5*2),rep("WPBC",5*5*2),
                               rep("SpamBase",5*5*2),rep("Ionosphere",5*5*2))
plotting_table[,"Feature Size"] = c(rep(c(50,100,150,200,250),10),rep(c(10,15,20,25,30),10),rep(c(10,20,30,40,50),10), rep(c(10,15,20,25,30),10))

write.table(plotting_table, "experiments/new_exec_time_exp/exec_time_all_weighted.csv")
