
library(fpc) # calinski-harabasz index
library(plyr) # mapvalues function
library(cluster) # daisy, PAM
library(dummy)
library(scales)
library(StatMatch) # gower.dist


dataset_type <- function(df){
  
  # !!!!!!!!!!!!!!!!OPRAVIT!!!!!!!!!!!!
  
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

kmeansHH <- function(x, nclusters = NULL){
  
  # inputs
  #   x - dataframe with continous scaled variables (numeric)
  #   nclusters - specify number of clusters. In case it is not specified, function finds optimum number of clusters
  
  # returns
  #   dataframe with input variables, plus additional variable 'cluster' - entry's assigned cluster
  #     and 'withinss' - clusters within sum of squared
  
  df <- x
  
  print(dataset_type(df))
  
  # rescaling variables
  
  #df[sapply(df, is.numeric)] <- lapply(df[sapply(df, is.numeric)], rescale)
  #df[sapply(df, is.integer)] <- lapply(df[sapply(df, is.integer)], rescale)
  
  if (dataset_type(df) == 'kok'){
    
    optimal_nclusters <- F
    
    if (is.null(nclusters)){
      
      optimal_nclusters <- T
      nclusters <- 10
      calHar <- numeric()
      dfs_list <- list()
      
    }
    
    for (i in 1:nclusters){
      
      variance <- sapply(df, var) # variables variance
      
      max_var <- variance[variance == max(variance)] # variable with maximum variance
      
      max_var_mean <- mean(df[[names(max_var)]]) # mean of this variable
      
      cent1 <- sapply(df[df[names(max_var)] >= 0 &
                           df[names(max_var)] < max_var_mean,], mean)
      cent2 <- sapply(df[df[names(max_var)] >= max_var_mean,], mean)
      cent <- rbind(cent1, cent2) # centroids for k-means
      
      km <- kmeans(x = df,
                   #algorithm = 'Lloyd',
                   centers = cent)
      
      if (i == 1){
        
        df$cluster <- km$cluster
        df$withinss <- ifelse(km$cluster == 1, km$withinss[1], km$withinss[2])
        df_toKeep <- df
        
      } else {
        
        if (i < nclusters){
          
          df$cluster <- ifelse(km$cluster == 1, max(df_toKeep$cluster) + 1, max(df_toKeep$cluster) + 2)
          
          df$withinss <- ifelse(km$cluster == 1, km$withinss[1], km$withinss[2])
          
          # replace cluster which was splitted with these clusters
          
          df_toKeep[row.names(df_toKeep) %in% row.names(df),] <- df
        }
        
      }
      
      # when finding optimum number of clusters
      
      if (optimal_nclusters == T){
        
        dfs_list <- append(dfs_list, list(df_toKeep))
        calHar <- append(calHar, calinhara(df_toKeep[, !(colnames(df_toKeep) %in% c('cluster', 'withinss'))],
                                           df_toKeep$cluster,
                                           cn = length(unique(df_toKeep$cluster))))
        
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
      
      if (i < nclusters){
        
        # get withinss of all clusters in df_toKeep
        
        withinss_byCluster <- by(df_toKeep[, colnames(df_toKeep) == 'withinss'],
                                 df_toKeep[, 'cluster'],
                                 mean)
        
        # get cluster with highest withinss
        
        max_cluster <- names(sapply(withinss_byCluster, max)[sapply(withinss_byCluster, max) == max(sapply(withinss_byCluster, max))])
        
        # this cluster will be as df into next iteration
        
        df <- df_toKeep[df_toKeep$cluster == max_cluster, !(colnames(df_toKeep) %in% c('cluster', 'withinss'))]
        
      } else {
        
        # on last iteration run k-means with all centroids, without further splits
        
        last_centroids <- by(df_toKeep[, !(colnames(df_toKeep) %in% c('cluster', 'withinss'))],
                             df_toKeep[, 'cluster'],
                             function(x){sapply(x, mean)})
        
        last_centroids <- matrix(unlist(last_centroids), ncol = length(df_toKeep) - 2, byrow = T)
        
        km <- kmeans(df_toKeep[, !(colnames(df_toKeep) %in% c('cluster', 'withinss'))],
                     centers = last_centroids)
        
        df_toKeep$cluster <- km$cluster
        
        
      }
      
    }
    
    return(df_toKeep)
    
  } else {
    
    if (nrow(df) < 1000){
      
      # gower distance matrix + PAM
      
      gower_dist_matrix <- as.matrix(daisy(x = df, metric = 'gower'))
      
      # determining optimum number of clusters via silhouette width
      
      if (is.null(nclusters)){
        
        sil_width <- numeric()
        
        nclusters <- 10
        
        for (i in 2:nclusters){
          
          pam_fit <- pam(x = gower_dist_matrix, diss = T, k = i)
          sil_width <- append(sil_width, pam_fit$silinfo$avg.width)
          
          print('silhouette widht:')
          print(sil_width)
          
          if (i >= 3){
            
            if ((sil_width[i - 1]/sil_width[i -2] < 1.1) | (sil_width[i - 1] < sil_width[i - 2])){
              
              nclusters <- which(sil_width == sil_width[i - 1])
              
              pam_fit <- pam(x = gower_dist_matrix, diss = T, k = nclusters)
              df$cluster <- pam_fit$clustering
              
              return(df)
              
            }
            
          }
          
        }
        
        # just in case situation above does not happen, run PAM with nclusters where max(silhouette width)
        
        nclusters <- which(sil_width == max(sil_width)) + 1
        
      }
      
      pam_fit <- pam(x = gower_dist_matrix, diss = T, k = nclusters)
      df$cluster <- pam_fit$clustering
      
      return(df)
      
      
    } else {
      
      # one hot encoding of categorical variables
      # or sample -> gower dist -> PAM -> gower distance of not sampled entrys to medoids of sample
      
      sampled_rows <- sample.int(nrow(df), 1000)
      
      df_sampled <- data.frame(df[sampled_rows,], row.names = sampled_rows)
      df_not_sampled <- data.frame(df[!(row.names(df) %in% row.names(df_sampled)),])
      
      gower_dist_matrix <- as.matrix(daisy(x = df_sampled, metric = 'gower'))
      
      if (is.null(nclusters)){
        
        sil_width <- numeric()
        
        nclusters <- 10
        
        for (i in 2:nclusters){
          
          pam_fit <- pam(x = gower_dist_matrix, diss = T, k = i)
          sil_width <- append(sil_width, pam_fit$silinfo$avg.width)
          
          print('silhouette widht:')
          print(sil_width)
          
          if (i >= 3){
            
            if ((sil_width[i - 1]/sil_width[i -2] < 1.1) | (sil_width[i - 1] < sil_width[i - 2])){
              
              nclusters <- which(sil_width == sil_width[i - 1])
              
              pam_fit <- pam(x = gower_dist_matrix, diss = T, k = nclusters)
              df_sampled$cluster <- pam_fit$clustering
              
              # calculate gower distance between pam_fit medoids and not sampled data
              # and assign cluster number of the nearest medoid
              
              df_sampled_medoids <- df_sampled[pam_fit$medoids, ]
              print(str(df_sampled_medoids))
              print(str(df_not_sampled))
              not_sampled_dist_medoids <- data.frame(t(gower.dist(data.x = df_sampled_medoids[, colnames(df_sampled_medoids) != 'cluster'],
                                                                  data.y = df_not_sampled)),
                                                     row.names = row.names(df_not_sampled))
              colnames(not_sampled_dist_medoids) <- row.names(df_sampled_medoids)
              closest_medoids_rownr <- apply(not_sampled_dist_medoids, 1, function(x){names(x)[x == min(x)]})
              df_not_sampled$cluster <- mapvalues(x = closest_medoids_rownr,
                                                  from = row.names(df_sampled_medoids),
                                                  to = df_sampled_medoids$cluster)
              
              df <- rbind(df_sampled, df_not_sampled)
              
              return(df)
              
            }
            
          }
          
        }
        
        # just in case situation above does not happen, run PAM with nclusters where max(silhouette width)
        
        #nclusters <- which(sil_width == max(sil_width)) + 1
        
      }
      
    }
    
  }
  
}

