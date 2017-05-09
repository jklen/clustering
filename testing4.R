
library(ggplot2)
library(fpc) # calinski-harabasz index, pamk
library(cluster) # pam, daisy (similarity matrix - gower distance)
library(plyr) # mapvalues function
library(FactoMineR) # MCA
library(homals) # MCA
library(ISLR) # college dataset
library(dplyr)
library(Rtsne) # tsne plot
library(scales)

iris_test_cats <- iris_test
iris_test_cats$Petal.Width.Sc.Cat <- ifelse(iris_test_cats$Petal.Width.Scaled >= 0 & iris_test_cats$Petal.Width.Scaled < 0.3,
                                           'pwscat1', ifelse(iris_test_cats$Petal.Width.Scaled >= 0.3 & iris_test_cats$Petal.Width.Scaled < 0.7,
                                                             'pwscat2', 'pwscat3'))
iris_test_cats$Sepal.Length.Sc.Cat <- ifelse(iris_test_cats$Sepal.Length.Scales >= 0 & iris_test_cats$Sepal.Length.Scales < 0.45,
                                            'slscat1', ifelse(iris_test_cats$Sepal.Length.Scales >= 0.45 & iris_test_cats$Sepal.Length.Scales < 0.6,
                                                              'slscat2', 'slscat3'))
iris_test_cats[5:6] <- lapply(iris_test_cats[5:6], factor)

####################################

itc_dist <- daisy(iris_test_cats) # distance !!! diamonds dataset (54000 x 10) was not able to calculate max 20000
itc_dist_matrix <- as.matrix(itc_dist) # distance to matrix

#####################################

# pam from cluster- with gower distance works well, but problem big data sets
# kmodes from klaR- worse than pam, in terms of clustering and time of processing, diamonds[1:20000,] did not finish in 5 mins
# kmeansvar from ClustOfVar - NO, clusters variables

# possibilities:
# one hot encoding
# coordinates from MCA as new variables
#   mca$ind$coord (individuals - add all dimensions to each entry and remove original categorical variables), or 
#   mca$var$coord (join with category and add as new columns, remove original categorical variables)


tablo_clust_cat <- read.table('2catspwsTabl4clust.csv', sep = ';', header = T)
colnames(tablo_clust_cat)[1] <- 'cluster'
tablo_clust_cat[4:7] <- lapply(tablo_clust_cat[4:7], 
                               function(x){as.numeric(gsub(pattern = ',', replacement = '.', as.character(x)))})


# clustering in tablo was done with vars: 2 categoricals + Petal.Width.Scaled

# conversion from json
#   I need example json data!!!
#   all strings, or JSON data types?
#   how to identify dimensions (if for example '1', '2', '3', '4', ... - factor? ordered factor? measure?)

# conversion to json


# variables validation

# testing and comparing to tableau

dataset_type <- function(df){
  
  class_vector <- sapply(df, class)
  
  numeric_cond <- (('numeric' %in% class_vector) | ('integer' %in% class_vector)) & 
    !('character' %in% class_vector) & !('factor' %in% class_vector)
  
  cat_cond <- (('character' %in% class_vector) | ('factor' %in% class_vector)) &
    !('numeric' %in% class_vector) & !('integer' %in% class_vector)
  
  if (numeric_cond){
    
    return('numeric')
    
  } else {
 
    if (cat_cond){
 
      return('categorical')
 
    } else {
  
      return('mixed')
 
    }

  }

  
}


prepare_vars <- function(x){
  
  # inputs
  #   x - dataframe
  # returns
  #   dataframe where categorical variables are converted with MCA
  #     into first 2 dimensions of individual coordinates/object scores, all variables are scaled to min-max/0-1
  
  df <- x

  df[sapply(df, is.character)] <- lapply(df[sapply(df, is.character)], as.factor) # character to factor
  df[sapply(df, is.ordered)] <- lapply(df[sapply(df, is.ordered)], as.integer) # ordered factor to integer
  df_cat <- df[sapply(df, is.factor)]

  if (length(df_cat) > 0){
    
    # conversion of factors to 2 MCA dimensions
  
    mca <- MCA(df_cat, ncp = 2, graph = F)
    df[sapply(df, is.factor)] <- NULL
    df <- cbind(df, mca$ind$coord)
    
  }
  
  df <- as.data.frame(lapply(df, rescale)) # scale all variables, min-max/0-1
  
  return(df)
  
}

kmeansHH3 <- function(x, nclusters = NULL){
  
  df <- prepare_vars(x)
  
  optimal_nclusters <- F
  
  if (is.null(nclusters)){
    
    optimal_nclusters <- T
    nclusters <- 10
    calHar <- numeric()
    dfs_list <- list()
    
  }
  
  for (i in 1:nclusters){
    
    print(paste('run: ', i))
    
    variance <- sapply(df, var) # variables variance
    print('variable variance:')
    print(variance)
    
    max_var <- variance[variance == max(variance)] # variable with maximum variance
    print('variable with highest variance:')
    print(max_var)
    
    max_var_mean <- mean(df[[names(max_var)]]) # mean of this variable
    print('mean of variable with highest variance:')
    print(max_var_mean)
    
    
    cent1 <- sapply(df[df[names(max_var)] >= 0 &
                         df[names(max_var)] < max_var_mean,], mean)
    cent2 <- sapply(df[df[names(max_var)] >= max_var_mean,], mean)
    cent <- rbind(cent1, cent2) # centroids for k-means
    
    print('centroids to get 2 clusters:')
    print(cent)
    
    km <- kmeans(x = df,
                 #algorithm = 'Lloyd',
                 centers = cent)
    
    print('clusters withinss:')
    print(km$withinss)
    
    if (i == 1){
      
      df$cluster <- km$cluster
      df$withinss <- ifelse(km$cluster == 1, km$withinss[1], km$withinss[2])
      df_toKeep <- df
      
    } else {
      
      if (i < nclusters){
        
        df$cluster <- ifelse(km$cluster == 1, max(df_toKeep$cluster) + 1, max(df_toKeep$cluster) + 2)
        print(str(df))
        
        df$withinss <- ifelse(km$cluster == 1, km$withinss[1], km$withinss[2])
        
        # replace cluster which was splitted with these clusters
        
        df_toKeep[row.names(df_toKeep) %in% row.names(df),] <- df
      }
      
    }
    
    print(paste('df_toKeep rows nr: ', nrow(df_toKeep)))

    # when finding optimum number of clusters
    # last kmeans run with all centroids like when specified number of clusters is omitted
    
    if (optimal_nclusters == T){
      
      df_toKeep$cluster <- mapvalues(df_toKeep$cluster,
                                     from = unique(df_toKeep$cluster),
                                     to = 1:length(unique(df_toKeep$cluster)))
      
      dfs_list <- append(dfs_list, list(df_toKeep))
      calHar <- append(calHar, calinhara(x = df_toKeep[, !(colnames(df_toKeep) %in% c('cluster', 'withinss'))],
                                         clustering = df_toKeep$cluster))
      print(paste('max clustering as input to calinhara: ', max(df_toKeep$cluster)))
      print('calinski-harabsz index OPT CLUSTERS:')
      print(calHar)
      
      if (i >= 2){
        
        if (calHar[i] < calHar[i - 1]){
          
          df_toKeep <- dfs_list[[i - 1]]
          
          # to have cluster numbering starting 1
          
          
          
          return(df_toKeep)
          
        }
        
      }
      
    }
    
    print('unique clusters:')
    print(unique(df_toKeep$cluster))
    print('number of unique clusters:')
    print(length(unique(df_toKeep$cluster)))
    
   
    if (i < nclusters){
      
      # get withinss of all clusters in df_toKeep
      
      withinss_byCluster <- by(df_toKeep[, colnames(df_toKeep) == 'withinss'],
                               df_toKeep[, 'cluster'],
                               mean)
      
      print('withinss by cluster:')
      print(withinss_byCluster)
      
      # get cluster with highest withinss
      
      max_cluster <- names(sapply(withinss_byCluster, max)[sapply(withinss_byCluster, max) == max(sapply(withinss_byCluster, max))])
      
      print('cluster name, which has maximum withinss')
      print(max_cluster)
      
      # this cluster will be as df into next iteration
      
      df <- df_toKeep[df_toKeep$cluster == max_cluster, !(colnames(df_toKeep) %in% c('cluster', 'withinss'))]
      
      print('structure of dataframe for next iteration')
      print(str(df))
      
    } else {
      
      # on last iteration run k-means with all centroids
      
      last_centroids <- by(df_toKeep[, !(colnames(df_toKeep) %in% c('cluster', 'withinss'))],
                           df_toKeep[, 'cluster'],
                           function(x){sapply(x, mean)})
      
      print('last centroids:')
      print(last_centroids)
      
      last_centroids <- matrix(unlist(last_centroids), ncol = length(df_toKeep) - 2, byrow = T)
      
      km <- kmeans(df_toKeep[, !(colnames(df_toKeep) %in% c('cluster', 'withinss'))],
                   centers = last_centroids)
      
      df_toKeep$cluster <- km$cluster
      
     
      
    }
    
  }
  
  # plotting
  
  # df_toPlot <- df_toKeep
  # df_toPlot$cluster <- as.factor(df_toPlot$cluster)
  # print(ggplot(aes(x = Petal.Width.Scaled, y = Petal.Length.Scaled), data = df_toPlot) +
  #         geom_point(aes(color = cluster, alpha = 0.1), size = 2) +
  #         xlim(0, 1.1) +
  #         ylim(0, 1.1))

  
  return(df_toKeep)
  
}

