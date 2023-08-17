library(crayon)


hog_method <- function(DB, gen_data ,od_name, hog_name, 
                       ODM_env, ho_list, ho_type, exec_time, directions = NA){
  #' @title Hidden Outlier Generation Method Class
  #' @description Definition of the HOG method class for R
  #' 
  #' Arguments:
  #' @param DB: Dataset used
  #' @param gen_data: Number of data generated
  #' @param od_name: Adersary name
  #' @param hog_name: HOG method used
  #' @param ODM_env: Environment containing all fitted detectors
  #' @param ho_list: Generated Hidden Outliers
  #' @param ho_type: Type of the generated hidden outliers (H1 or H2)
  #' @param exec_time: Execution time of the generation
  #' @param directions: Vectors sampled inside the sphere (defaults to none)
  structure(list(DB, gen_data ,od_name, hog_name,
                 ODM_env, ho_list, ho_type, exec_time, directions),
            names = c("DB", "n_data_gen" ,"odm_name", "hog_name",
                      "ODM_env", "ho_list", "ho_type", "exec_time", "directions"),
            class = "hog_method")
}

 
print.hog_method <- function(method){
  print(glue("{bold('Hidden Outlier Generation Method Object')}\n\n
             Adversary used: {method$odm_name}\n
             Synthetic HO generation method employed: {method$hog_name}.\n\n
             {cyan('Use summary() for a detailed description')}\n
             {silver('All the fitted ODMs for each subspace can be found in the 
             enviroment ODM_env (<your_method>$ODM_env).')}"))}
          
summary.hog_method <- function(method){
  print(glue("{bold('Hidden Outlier Generation Method Object')}\n\n
             Adversary used: {method$odm_name}\n
             Synthetic HO generation method employed: {method$hog_name}.
             
\n\n {underline('Database summary')}:
\n
             * {italic('Number of features:')} {ncol(method$DB)-2}
             * {italic('Total number of data points:')} {nrow(method$DB)}
             * {italic('Total amount of syntetic data generated:')} {
             nrow(method$ho_type)}
             \t          {italic('...of which hidden outliers:')} {
             if(class(method$ho_list)[1]!='matrix' && 
             length(method$ho_list)!= 0){
             1
             }  
             if(class(method$ho_list)[1]=='matrix'){
             nrow(method$ho_list)}else{
             0
             }
             }
             * {italic('Number of H1 outliers')}: {
              as.data.frame(method$ho_type) %>%
              filter(V1 == 'H1') %>%
              count()}
             * {italic('Number of H2 outliers:')} {
              as.data.frame(method$ho_type) %>%
              filter(V1 == 'H2') %>%
              count()}.

\n\n {cyan('Total execution time:')} {method$exec_time}."))
}
