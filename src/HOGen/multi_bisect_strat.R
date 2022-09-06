source("src/ODM/inference_methods.R")
source("src/HOGen/bisect_strat.R")


interval_check <- function(L, method, x,...){
  D = 1:(ncol(DB)-2)
  segmentation_points = seq(0, L, length=5)
  check = array(0,length(segmentation_points))
  
  i = 1
  for (c in segmentation_points){
    if (inference(c*x + colMeans(DB[2:(ncol(DB) - 1)]), S = D, method,...)){
      check[i] = 1
    }
    i = i +1
   }
  
  previous = check[1]
  interval = list()
  for (i in 1:length(check)){
    if(check[i] != previous){
      interval = append(interval,list(c(segmentation_points[i-1],
                                        segmentation_points[i])))
    }
    previous = check[i]
  }
  
  return(interval)
}


multi_bisect <- function(x, L, iternum = 1000, method, verb = T,...){
  
  interval = interval_check(L,method,x,...)
  interval = sample(interval, 1)[[1]]
  
  a = interval[1]; b = interval[2]
  for (i in 1:iternum){
    c = (b+a)/2
    
    check_if_outlier = f(c*x + colMeans(DB[2:(ncol(DB) - 1)]), 
                         method = method,...)
    outlier_indicator = check_if_outlier[[1]]
    outlier_type = check_if_outlier[[2]]
    
    if (outlier_indicator > 0){b = c} #Directional vector 
    #(x*c) + an origin (colMeans) 
    else if (outlier_indicator < 0){a = c}
    else{ print(glue("x in {outlier_type}")); break}
  }
  return(list(c, outlier_type))
}


main_multibisect <-function(B=100, method="mahalanobis", seed=F,
                            verb=T, dev_opt=F,...){
  if(class(seed) == 'numeric'){
    set.seed(seed)
  }
  if(seed == T){
    seed = sample(1:10^9,1)
    print("A seed has been generated")
    dir.create("seed", showWarnings=F)
    write.table(seed, file="seed/seed.txt", sep=" ")
    set.seed(seed)
  }
  if(dev_opt == T){warning(glue("The developer option has been activated. 
                                This means that the ODM environment will not 
                                be deleted, and therefore some collusion might
                                occur with following experiments. Poceed with
                                caution
                                "))}
  
  x_list = runif_on_sphere(n=B, d=ncol(DB) - 2, r=1) #sample the 
  #directional vectors
  L = max(sqrt(rowSums(DB[2:(ncol(DB)-1)]^2)))
  hidden_x_list = matrix(0, nrow=nrow(x_list), ncol=ncol(x_list))
  hidden_x_type = matrix(0, nrow=nrow(x_list), ncol=1)
  
  tic()
  for (i in 1:nrow(x_list)){
    bisection_results = multi_bisect(L=L, x=x_list[i,], 
                                     method=method, verb=verb)
    hidden_c = bisection_results[[1]]
    outlier_type = bisection_results[[2]]
    
    if(outlier_type %in% c("H1", "H2")){ 
      hidden_x_list[i,] = hidden_c*(x_list[i,]) + colMeans(DB[2:(ncol(DB)-1)])
      hidden_x_type[i,] = outlier_type
    }
  }
  exec_time = toc()
  exec_time = exec_time$callback_msg
  
  hidden_x_list = hidden_x_list[rowSums(hidden_x_list) != 0,]
  gen_result = hog_method(DB, B, method, "multi_bisection", 
                          ODM_env, hidden_x_list, hidden_x_type, exec_time)
  
  if(dev_opt == F){
    rm(ODM_env, envir = globalenv())}
  return(gen_result)
  
  
  
  
  
}



