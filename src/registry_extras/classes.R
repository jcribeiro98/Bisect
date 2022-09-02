library(crayon)


hog_method <- function(DB, B ,od_method, hog_method, 
                       ODM_env, ho_list, ho_type, exec_tyme){
  structure(list(DB, B ,od_name, hog_name,
                 ODM_env, ho_list, ho_type, exec_tyme),
            names = c("DB", "n_data_gen" ,"odm_name", "hog_name",
                      "ODM_env", "ho_list", "ho_type", "exec_time"),
            class = "hog_method")
}

print.hog_method <- function(method){
  print(glue("{bold('Hidden Outlier Generation Method Oject')}\n\n
             Outlier Gen Method used: {method$odm_name}\n
             Synthetic HO generation method employed: {method$hog_name}.\n\n
             {cyan('Use summary() for a detailed description')}\n
             {silver('All the fitted ODMs for each subspace can be found in the 
             enviroment ODM_env (<your_method>$ODM_env).')}"))}
          
summary.hog_method <- function(method){
  print(glue("{bold('Hidden Outlier Generation Method Oject')}\n\n
             Outlier detection method used: {method$odm_name}\n
             Synthetic HO generation method employed: {method$hog_name}.
             
             \n\n Database summary:
             \n\t\t number of features: {ncol(method$DB)-2}
             \n\t\t total number of data points: {nrow(method$DB)}
             \n\t\t total amount of syntetic data generated: {method$n_data_gen}
              \n\t\t\t ...of of which hidden outliers: {nrow(method$ho_list)}
             \n\t\t number of H1 outliers: {as.data.frame(ho_type) %>%
              filter(ho_type == 'H1') %>%
              count()}
             \n\t\t number of H2 outliers: {as.data.frame(ho_type) %>%
              filter(ho_type == 'H2') %>%
              count()}.

             \n\n Total execution time: {method$exec_time}."))
}