library(sets)
library(glue)
library(car)
library(dplyr)
library(hrbrthemes)
library(ggplot2)
library(latex2exp)
source('src/HOGen/bisect_strat.R')
source('src/HOGen/hidden_strat.R')

#------------------------------------------------------------------------------#
#                         n-dimensional Linear Search                          #
#------------------------------------------------------------------------------#

#### 1. 3d sampling of HO with bisect ####
n = 3
nsamples = 1000


norMatrix = matrix(0,nrow = nsamples, ncol = n)
for (i in 1:n) {norMatrix[,i] = rnorm(nsamples, runif(1,min=3,max = 10), 
                                      runif(1,min = 2, max=5))}

DB = DB_gen(norMatrix)

x_bisect = main(B=10)
png("images/bisect_3d")

scp<-scatterplot3d(DB[,2:4], pch=16, box = F, angle = 60, color = 'black',)
scp$points3d(x_bisect[,1],x_bisect[,2],x_bisect[,3], pch = 16, col = "#1B9E77",
             type = 'h')
legend("topright", legend = c("Inliers", "Hidden Outliers"), fill = 
         c('black', "#1B9E77"))

dev.off()


##### 3d proportion of H1 outliers vs H2 outliers #####
x_bisect = main(B=100)

colnames(x_bisect) = c('y1', 'y2', 'y3')
x_in = matrix(0,1,nrow(x_bisect)) #Generating the colors for the plot
ih2 = 0
for (i in 1:nrow(x_bisect)){
if (f(x_bisect[i,], verb = T, h_index = T)[2] == "h1"){x_in[i] = 'steelblue'}
else{x_in[i] = "#D95F02"; ih2 = ih2 + 1}
}

png("images/bisect_H1vsH2")

scp<-scatterplot3d(DB[,2:4], pch=16, box = F, angle = 60, color = 'black',)
scp$points3d(x_bisect[,1],x_bisect[,2],x_bisect[,3], pch = 16, col = x_in ,
             type = 'h')

legend("topright", legend = c("Inliers", 
                              glue("H1 Hidden Outliers: {100-ih2}/100"), 
                              glue("H2 Hidden Outliers: {ih2}/100")), 
       fill =  c('black', "steelblue", "#D95F02"))
dev.off()


#### 2. 9d sampling of HO using bisect ####
n = 9
nsamples = 1000


norMatrix = matrix(0,nrow = nsamples, ncol = n)
for (i in 1:n) {norMatrix[,i] = rnorm(nsamples, runif(1,min=3,max = 10), 
                                      runif(1,min = 2, max=5))}

DB = DB_gen(norMatrix)

start_time = Sys.time()
x_bisect = main(B=10)
end_time = Sys.time()
total_time_taken = end_time - start_time

start_time = Sys.time()
x_bisect = main(B=100)
end_time = Sys.time()
total_time_taken2 = end_time - start_time

png("images/freq_B10")

barplot(c(1,0), names.arg = c('H1', 'H2'), col = c("steelblue","#D95F02"), 
        main = 'Frequency of H1 outliers vs H2 outliers, for B=10, d=9', sub = 
          glue('Time taken: {sprintf(total_time_taken,fmt = "%#.2f")} minutes'))
dev.off()

png("images/freq__B100")
barplot(c(1,0), names.arg = c('H1', 'H2'), col = c("steelblue","#D95F02"), 
        main = 'Frequency of H1 outliers vs H2 outliers, for B=100, d=9', sub = 
          glue('Time taken: {sprintf(total_time_taken2,fmt = "%#.2f")} minutes'))

dev.off()



#### 3. 9d sampling of HO using hidden ####


x_hidden = main_hidden(B=100, eps = 0.3)



##### The problem with eps #####
B=100
load('examp_eps.RData')
x_hidden1 = main_hidden(B=B, eps = 0.3)
x_hidden2 = main_hidden(B=B, eps = 0.7)

n_x_hidden1 = ncol(x_hidden1)
n_x_hidden2 = ncol(x_hidden2)

png("images/eps_comp")
bar = barplot(c(n_x_hidden1/B, n_x_hidden2/B), 
              names.arg = c('eps = 0.3','eps = 0.7'), 
              col = c("steelblue","#D95F02"), 
              main = 'Frequency of HO for different eps in hidden')




text(bar, c(n_x_hidden1/B - .03, n_x_hidden2/B-.01),  
     paste("Outlier Count: ", c(n_x_hidden1, n_x_hidden2), sep=""), 
     family = "URWTimes")
text(bar, c(n_x_hidden1/B - .029, n_x_hidden2/B - .0085),  
     paste("Outlier Count: ",c(n_x_hidden1, n_x_hidden2), sep=""), 
     family = "URWTimes", col="white")
dev.off()

###### Graph ######
B=100
load('examp_eps.RData')

n_vector = matrix(0, nrow=1,ncol=9)

for (i in 1:9){
  print(glue('Working on experiment: {i}\n'))
  n_vector[i] = nrow(main_hidden(B, eps = (i/10)))/B
  print(glue('Finished. \n Result: {n_vector[i]*B} outliers found for eps = 
             {i/10}\n \n'))
}

png("images/eps_comparison_9d.png")
as.data.frame(cbind(t(n_vector),(1:9)/10)) %>%
    select('n_H0'='V1', 'eps' = 'V2') %>% 
    ggplot( aes(x=eps, y=n_H0)) +
      geom_line(color = 'darkgrey') +
      geom_point(shape=16, col = "#66C2A5", size = 2.5) +
      geom_point( aes(x = .3, y = .46), col = "#FC8D62", size = 3) +
      annotate(x = .39, y = 0.46, "text", label = c(TeX("$\\epsilon= 0.3$")), 
               col = alpha("black", .3)) + 
      labs( y = "Number of Hidden Outliers", x = TeX('$\\epsilon$'),subtitle = 
            "For B = 100, and 9 features", 
            title = TeX("Number of Hidden Outliers found for each $\\epsilon$ value"))+
      theme_ipsum(axis_title_size = 17 )  
dev.off()


#### 4. Hidden vs Bisect ####

##### Execution time 9 features #####

###### eps = 0.7 ######

start_time = Sys.time()
n_bisect = nrow(main(B=100))
end_time = Sys.time()
total_time_taken_bisect = end_time - start_time

start_time = Sys.time()
n_hidden1 = nrow(main_hidden(B=100, eps = 0.7))
end_time = Sys.time()
total_time_taken_hidden1 = end_time - start_time

start_time = Sys.time()
n_hidden2 = nrow(main_hidden(B=200, eps = 0.7))
end_time = Sys.time()
total_time_taken_hidden2 = end_time - start_time

start_time = Sys.time()
n_hidden5 = nrow(main_hidden(B=500, eps = 0.7))
end_time = Sys.time()
total_time_taken_hidden5 = end_time - start_time

png("images/nout_and_runtime_graph.png")
bar = barplot(c(n_hidden5,n_bisect,n_hidden2,n_hidden1),
        names.arg = c('Hidden B=500', 'Bisect B=100', 'Hidden B=200', 
                      'Hidden B=100'),
        col = c("#66C2A5","#FC8D62","#8DA0CB","#E78AC3"),
        main = 'Number of Hidden Outliers found',
        sub = 'Hidden performed with eps = 0.7')
text(bar, c(n_hidden5 + 1.5 ,n_bisect - 2,n_hidden2 + 1.5,n_hidden1 + 1.5), 
     paste("Time taken (mins): ",
           c(sprintf(total_time_taken_hidden5,fmt = "%#.2f"),
             sprintf(total_time_taken_bisect,fmt = "%#.2f"), 
             sprintf(total_time_taken_hidden2,fmt = "%#.2f"), 
             sprintf(total_time_taken_hidden1,fmt = "%#.2f")), 
           sep = ""),
     family = "URWTimes",
     cex = .9)        
dev.off()


##### Execution time 15 features #####
###### eps = 0.7 ######
n = 15
nsamples = 1000


norMatrix = matrix(0,nrow = nsamples, ncol = n)
for (i in 1:n) {norMatrix[,i] = rnorm(nsamples, runif(1,min=3,max = 10), 
                                      runif(1,min = 2, max=5))}

DB = DB_gen(norMatrix)

start_time = Sys.time()
n_bisect = nrow(main(B=100))
end_time = Sys.time()
total_time_taken_bisect = end_time - start_time

start_time = Sys.time()
n_hidden1 = nrow(main_hidden(B=100, eps = 0.7))
end_time = Sys.time()
total_time_taken_hidden1 = end_time - start_time

start_time = Sys.time()
n_hidden2 = nrow(main_hidden(B=200, eps = 0.7))
end_time = Sys.time()
total_time_taken_hidden2 = end_time - start_time

start_time = Sys.time()
n_hidden5 = nrow(main_hidden(B=500, eps = 0.7))
end_time = Sys.time()
total_time_taken_hidden5 = end_time - start_time

png("images/15feat_nout_and_runtime_graph.png")
bar = barplot(c(n_hidden5,n_bisect,n_hidden2,n_hidden1),
              names.arg = c('Hidden B=500', 'Bisect B=100', 'Hidden B=200', 
                            'Hidden B=100'),
              col = c("#66C2A5","#FC8D62","#8DA0CB","#E78AC3"),
              main = 'Number of Hidden Outliers found',
              sub = 'Hidden performed with eps = 0.7')
text(bar, c(n_hidden5 + 1.5 ,n_bisect - 2,n_hidden2 + 1.5,n_hidden1 + 1.5), 
     paste("hours: ",
           c(sprintf(total_time_taken_hidden5,fmt = "%#.2f"),
             sprintf(total_time_taken_bisect,fmt = "%#.2f"), 
             sprintf(total_time_taken_hidden2,fmt = "%#.2f"), 
             sprintf(total_time_taken_hidden1,fmt = "%#.2f")), 
           sep = ""),
     family = "URWTimes",
     cex = .9)        
dev.off()

