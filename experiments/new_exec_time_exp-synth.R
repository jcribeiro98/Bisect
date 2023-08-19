# Execution Time (Synthetic)  ---------------------------------------------
library(miceadds)
source.all("src/HOGen/")
library(farff)
library(rearrr)
library(dplyr)
library(ggplot2)

gen_clusters <- function(nsamples,n,clusters){
  #' @title Generates a normal clustered dataset
  #' Arguments:
  #' @param nsamples: Number of samples in the dataset
  #' @param n: Number of features 
  #' @param clusters: Number of clusters
  
  norMatrix <- matrix(0, nrow = nsamples, ncol = n)
  for (cluster in 1:clusters){
    for (i in 1:n) {
      if(cluster == 1){
        norMatrix[1:(nsamples / clusters), i] <- rnorm( nsamples/clusters, 
                                                        runif(1, min = 1, max = 10),
                                                        runif(1, min = 1, max = 10))
      }else{
      norMatrix[((cluster - 1)*(nsamples / clusters) + 1):(cluster * 
                                                   nsamples / clusters), i] <- rnorm(
        nsamples/clusters, 
        runif(1, min = -3, max = 10),
        runif(1, min = 1, max = 11)
      )
      }}
  }
  return(norMatrix)
}


time_exper <- function(seed_vector = c(47385123, 12345, 
                                       54321, 010101, 33),
                       size_vector = c(7,15,30,50,100,150,200),
                       clusters = 1,
                       ...){
  #' @title Time experiment with synthetic data
  #' @description Performs the time experiments with different normal 
  #' distributions obtained by 'gen_clusters'
  #' 
  #' Arguments:
  #' @param seed_vectors: Vector of seeds for the repetitions
  #' @param size_vector: Vector of features
  #' @param clusters: Number of clusters
  
  time_matrix_bisect = matrix(0, ncol = length(seed_vector), 
                              nrow = length(size_vector))
  time_matrix_hidden = time_matrix_bisect
  i = 1
  for (size in size_vector){
    j = 1
    for(seed in seed_vector){
      set.seed(seed)
      DB_gen(gen_clusters(1000,size,clusters))
      hog_model = main_multibisect(gen_points = 500, method = "pyod_LOF", 
                                   seed = seed, num_workers = 1, type = "weighted", dev_opt = TRUE)
      time_matrix_bisect[i, j] = as.numeric(gsub("[secelapsed ]","", hog_model$exec_time))
      set.seed(seed)
      hog_model = main_hidden(gen_points = 500, eps = .5, method = "pyod_LOF", num_workers = 1, until_gen = TRUE)
      time_matrix_hidden[i, j] = as.numeric(gsub("[secelapsed ]","", hog_model$exec_time))
      j = j + 1
    }
    i = i + 1
  }
  return(cbind(time_matrix_bisect,time_matrix_hidden))
}

res_cluster1 = time_exper(clusters = 1)
colnames(res_cluster1) = c("lo seed1","lo seed2","lo seed3","lo seed4",
                           "lo seed5","eps1 seed1","eps1 seed2","eps1 seed3",
                           "eps1 seed4","eps1 seed5")
rownames(res_cluster1) = c(7,15,30,50,100,150,200)
save.image("experiments/new_exec_time_exp/new_generation/new_synth_cluster1_time_exp.RData")

res_cluster2 = time_exper(clusters = 2)
colnames(res_cluster2) = c("lo seed1","lo seed2","lo seed3","lo seed4",
                           "lo seed5","eps1 seed1","eps1 seed2","eps1 seed3",
                           "eps1 seed4","eps1 seed5")
rownames(res_cluster2) = c(7,15,30,50,100,150,200)
save.image("experiments/new_exec_time_exp/new_generation/new_synth_cluster2_time_exp.RData")

res_cluster5 = time_exper(clusters = 5, size_vector = c(7,15,30,50))
colnames(res_cluster5) = c("lo seed1","lo seed2","lo seed3","lo seed4",
                           "lo seed5","eps1 seed1","eps1 seed2","eps1 seed3",
                           "eps1 seed4","eps1 seed5")
rownames(res_cluster5) = c(7,15,30,50)
save.image("experiments/new_exec_time_exp/new_generation/0.5_new_synth_cluster5_time_exp_part2.RData")


## Graphing  ---------------------------------------------------------------

res_cluster1 = read.csv("experiments/new_exec_time_exp/new_generation/res_cluster1.csv", sep = "")

plotting_table = data.frame(matrix(0,7*5*2,3))
colnames(plotting_table) = c("Time", "Method", "Feature Size")
plotting_table[,'Method'] = c(rep("Bisect",5*7),rep("Hidden",5*7))
plotting_table[,'Time'] = c(res_cluster1[,1],res_cluster1[,2],res_cluster1[,3],
                           res_cluster1[,4],res_cluster1[,5],res_cluster1[,6],
                           res_cluster1[,7],res_cluster1[,8],res_cluster1[,9],
                           res_cluster1[,10])
plotting_table[,"Feature Size"] = rep(c(7,15,30,50,100,150,200),10)
write.table(plotting_table,"experiments/new_exec_time_exp/new_generation/res_cluster1_plotting_table.csv")


res_cluster2 = read.csv("experiments/new_exec_time_exp/new_generation/res_cluster2.csv", sep = "")

plotting_table = data.frame(matrix(0,7*5*2,3))
colnames(plotting_table) = c("Time", "Method", "Feature Size")
plotting_table[,'Method'] = c(rep("Bisect",5*7),rep("Hidden",5*7))
plotting_table[,'Time'] = c(res_cluster2[,1],res_cluster2[,2],res_cluster2[,3],
                            res_cluster2[,4],res_cluster2[,5],res_cluster2[,6],
                            res_cluster2[,7],res_cluster2[,8],res_cluster2[,9],
                            res_cluster2[,10])
plotting_table[,"Feature Size"] = rep(c(7,15,30,50,100,150,200),10)
write.table(plotting_table,"experiments/new_exec_time_exp/new_generation/res_cluster2_plotting_table.csv")


res_cluster5 = read.csv("experiments/new_exec_time_exp/new_generation/res_cluster5.csv", sep = "")

plotting_table = data.frame(matrix(0,7*5*2,3))
colnames(plotting_table) = c("Time", "Method", "Feature Size")
plotting_table[,'Method'] = c(rep("Bisect",5*7),rep("Hidden",5*7))
plotting_table[,'Time'] = c(res_cluster5[,1],res_cluster5[,2],res_cluster5[,3],
                            res_cluster5[,4],res_cluster5[,5],res_cluster5[,6],
                            res_cluster5[,7],res_cluster5[,8],res_cluster5[,9],
                            res_cluster5[,10])
plotting_table[,"Feature Size"] = rep(c(7,15,30,50,100,150,200),10)
write.table(plotting_table,"experiments/new_exec_time_exp/new_generation/res_cluster5_plotting_table.csv")



plotting_table = data.frame(matrix(0,7*5*2*3,4))
colnames(plotting_table) = c("Time", "Method", "Feature Size", "Clusters")
plotting_table[,'Method'] = c(rep("Bisect",5*7),rep("Hidden",5*7))
plotting_table[,'Time'] = c(c(res_cluster1[,1],res_cluster1[,2],res_cluster1[,3],
                              res_cluster1[,4],res_cluster1[,5],res_cluster1[,6],
                              res_cluster1[,7],res_cluster1[,8],res_cluster1[,9],
                              res_cluster1[,10]),
                            c(res_cluster2[,1],res_cluster2[,2],res_cluster2[,3],
                              res_cluster2[,4],res_cluster2[,5],res_cluster2[,6],
                              res_cluster2[,7],res_cluster2[,8],res_cluster2[,9],
                              res_cluster2[,10]),
                            c(res_cluster5[,1],res_cluster5[,2],res_cluster5[,3],
                              res_cluster5[,4],res_cluster5[,5],res_cluster5[,6],
                              res_cluster5[,7],res_cluster5[,8],res_cluster5[,9],
                              res_cluster5[,10]))
plotting_table[,"Feature Size"] = rep(c(7,15,30,50,100,150,200),10*3)
plotting_table[,"Clusters"] = c(rep(1,7*5*2),rep(2,7*5*2),rep(5,7*5*2))
write.table(plotting_table,"experiments/new_exec_time_exp/new_generation/res_cluster5_plotting_table.csv")

max(plotting_table %>% filter( Method == "Hidden") %>% select(Time)/
  plotting_table %>% filter( Method == "Bisect") %>% select(Time))

max(plotting_table %>% filter( Method == "Hidden" & Clusters == 5 & `Feature Size` == 200) %>% select(Time)/
      plotting_table %>% filter( Method == "Bisect" & Clusters == 5 & `Feature Size` == 200) %>% select(Time))

