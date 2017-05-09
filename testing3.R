
library(ggplot2)
library(fpc) # calinski-harabasz index
library(plyr) # mapvalues function
library(FactoMineR) # MCA

kmeansHH3 <- function(x, nclusters = NULL){
  
  df <- x
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
    
    if (optimal_nclusters == T){
      
      dfs_list <- append(dfs_list, list(df_toKeep))
      calHar <- append(calHar, calinhara(df_toKeep[, !(colnames(df_toKeep) %in% c('cluster', 'withinss'))],
                                         df_toKeep$cluster,
                                         cn = length(unique(df_toKeep$cluster))))
      print('calinski-harabsz index:')
      print(calHar)
      
      if (i >= 2){
        
        if (calHar[i] < calHar[i - 1]){
          
          df_toKeep <- dfs_list[[i - 1]]
          
          # to have cluster numbering starting 1
          
          df_toKeep$cluster <- mapvalues(df_toKeep$cluster,
                                         from = unique(df_toKeep$cluster),
                                         to = 1:length(unique(df_toKeep$cluster)))
          
          return(df_toKeep)
          
        }
        
      }
      
    }
    
    print('unique clusters:')
    print(unique(df_toKeep$cluster))
    print('number of unique clusters:')
    print(length(unique(df_toKeep$cluster)))
    
    # plotting
    
    df_toPlot <- df_toKeep
    df_toPlot$cluster <- as.factor(df_toPlot$cluster)
    print(ggplot(aes(x = Petal.Width.Scaled, y = Petal.Length.Scaled), data = df_toPlot) +
      geom_point(aes(color = cluster, alpha = 0.1), size = 2) +
      xlim(0, 1.1) +
      ylim(0, 1.1) +
      ggtitle(i))
    
    if (i < nclusters){
    
      # get withinss of all clusters in df_toKeep

      withinss_byCluster <- by(df_toKeep[, colnames(df_toKeep) == 'withinss'],
                               df_toKeep[, 'cluster'],
                               mean)
      
      print('withinss by cluster:')
      print(withinss_byCluster)
      
      # get cluster with highest withinss

      max_cluster <- names(sapply(withinss_byCluster, max)[sapply(withinss_byCluster, max) == max(sapply(withinss_byCluster, max))])
      
      print('cluster name, which contains variable with maximum variance')
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
      
      last_centroids <- matrix(unlist(last_centroids), ncol = 4, byrow = T)
      
      km <- kmeans(df_toKeep[, !(colnames(df_toKeep) %in% c('cluster', 'withinss'))],
                   centers = last_centroids)
      
      df_toKeep$cluster <- km$cluster
      
      # plotting
      
      df_toPlot <- df_toKeep
      df_toPlot$cluster <- as.factor(df_toPlot$cluster)
      print(ggplot(aes(x = Petal.Width.Scaled, y = Petal.Length.Scaled), data = df_toPlot) +
              geom_point(aes(color = cluster, alpha = 0.1), size = 2) +
              xlim(0, 1.1) +
              ylim(0, 1.1) +
              ggtitle("last iteration"))
      
    }

  }
  
  return(df_toKeep)
  
}

