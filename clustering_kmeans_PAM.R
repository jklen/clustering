
library(fpc) # calinski-harabasz index
library(plyr) # mapvalues function
library(cluster) # daisy, PAM
library(dummy)
library(scales) # rescale
library(StatMatch) # gower.dist
library(ggplot2)
library(tidyr)
library(dplyr)

dataset_type <- function(x){
  
  df <- x
  
  numeric_var_count <- sum(sapply(df, function(x){is.numeric(x) | is.integer(x)}))
  dim_var_count <- sum(sapply(df, function(x){is.ordered(x) | is.factor(x)}))
  
  if (numeric_var_count > 0 & dim_var_count == 0){
    
    return('numeric')
    
  } else {
    
    if (dim_var_count > 0 & numeric_var_count == 0){
      
      return('categorical')
      
    } else {
      
      return('mixed')
      
    }
    
  }
  
}

kmPAM <- function(x, nclusters = NULL){
  
  df <- x
  df <- na.omit(df)
  
  # get min and max values for each integer or numeric variable, for rescaling back before returning final df
  
  min_values <- sapply(df[sapply(df, function(x) is.numeric(x) | is.integer(x))], min)
  max_values <- sapply(df[sapply(df, function(x) is.numeric(x) | is.integer(x))], max)
  
  # rescaling variables
  
  df[sapply(df, function(x) is.numeric(x) | is.integer(x))] <- 
    lapply(df[sapply(df, function(x) is.numeric(x) | is.integer(x))], rescale)
  
  if (dataset_type(df) == 'numeric'){
    
    optimal_nclusters <- F
    
    if (is.null(nclusters)){
      
      optimal_nclusters <- T
      nclusters <- 10
      calHar <- numeric()
      dfs_list <- list()
      
    }
    
    for (i in 1:nclusters){
      
       variance <- sapply(df, var) # variables variance
      
      # too low rows number
      
      if (is.na(variance) & (i >= 2)){
        
        # rescaling back variables (min_values and max_values have same variable names)
        
        df_toKeep[names(min_values)] <- 
          sapply(names(min_values), function(x){rescale(x = df_toKeep[x], from = c(0, 1), 
                                                        to = c(min_values[names(min_values) == x], 
                                                               max_values[names(max_values) == x]))})
        
        df_toKeep$cluster <- as.factor(df_toKeep$cluster)
        
        return(df_toKeep[, colnames(df_toKeep) != 'withinss'])
        
      }
      
      max_var <- variance[variance == max(variance)] # variable with maximum variance
      
      # if more variables have same max. variance
      
      if (length(max_var) >= 2){
        
        max_var <- max_var[1]
        
      }
      max_var_mean <- mean(df[[names(max_var)]]) # mean of this variable
      
      cent1 <- sapply(df[df[names(max_var)] >= 0 &
                           df[names(max_var)] < max_var_mean, colnames(df), drop = F], mean)
      cent2 <- sapply(df[df[names(max_var)] >= max_var_mean, colnames(df), drop = F], mean)
      cent <- rbind(cent1, cent2) # centroids for k-means
      
      # if there are no splits possible to continue (2 entries), put them in one cluster, run kmeans with all centroids
      # and finish
      
      if (nrow(df) > 2){
        
        km <- kmeans(x = df,
                     #algorithm = 'Lloyd',
                     centers = cent)
        
      } else {
        
        # 2 rows in dataset
        
        if (i == 1){
          
          df$cluster <- c(1,2)
          
          df$cluster <- as.factor(df$cluster)
          
          return(df)
          
        }
        
        df$cluster <- max(df_toKeep$cluster) + 1
        df_toKeep[row.names(df_toKeep) %in% row.names(df),] <- df
        
        
        last_centroids <- by(df_toKeep[, !(colnames(df_toKeep) %in% c('cluster', 'withinss')), drop = F],
                             df_toKeep[, 'cluster'],
                             function(x){sapply(x, mean)})
        
        last_centroids <- matrix(unlist(last_centroids), ncol = length(df_toKeep) - 2, byrow = T)
        
        km <- kmeans(df_toKeep[, !(colnames(df_toKeep) %in% c('cluster', 'withinss'))],
                     centers = last_centroids)
        
        df_toKeep$cluster <- km$cluster
        
        df_toKeep$withinss <- NULL
        
        # rescaling back variables (min_values and max_values have same variable names)
        
        df_toKeep[names(min_values)] <- 
          sapply(names(min_values), function(x){rescale(x = df_toKeep[x], from = c(0, 1), 
                                                        to = c(min_values[names(min_values) == x], 
                                                               max_values[names(max_values) == x]))})
        
        df_toKeep$cluster <- as.factor(df_toKeep$cluster)
        
        return(df_toKeep)
        
      }
      
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
            
            df_toKeep$withinss <- NULL
            
            # rescaling back variables (min_values and max_values have same variable names)
            
            df_toKeep[names(min_values)] <- 
              sapply(names(min_values), function(x){rescale(x = df_toKeep[x], from = c(0, 1), 
                                                            to = c(min_values[names(min_values) == x], 
                                                                   max_values[names(max_values) == x]))})
            
            df_toKeep$cluster <- as.factor(df_toKeep$cluster)
            
            return(df_toKeep)
            
          }
          
        }
        
      }
      
      if (i < nclusters){
        
        # when too low nubmer of rows
        
        if ((nrow(df) == nrow(df_toKeep)) & (i >= 2)){
          
          # rescaling back variables (min_values and max_values have same variable names)
          
          df_toKeep[names(min_values)] <- 
            sapply(names(min_values), function(x){rescale(x = df_toKeep[x], from = c(0, 1), 
                                                          to = c(min_values[names(min_values) == x], 
                                                                 max_values[names(max_values) == x]))})
          df_toKeep$cluster <- as.factor(df_toKeep$cluster)
          
          return(df_toKeep[, colnames(df_toKeep) != 'withinss'])
          
        }
        
        # get withinss of all clusters in df_toKeep
        
        withinss_byCluster <- by(df_toKeep[, colnames(df_toKeep) == 'withinss'],
                                 df_toKeep[, 'cluster'],
                                 mean)
        
        # get cluster with highest withinss
        
        max_cluster <- names(sapply(withinss_byCluster, max)[sapply(withinss_byCluster, max) == max(sapply(withinss_byCluster, max))])
   
        # this cluster will be as df into next iteration
        
        df <- df_toKeep[df_toKeep$cluster == max_cluster, !(colnames(df_toKeep) %in% c('cluster', 'withinss')), drop = F]
        
      } else {
        
        # on last iteration run k-means with all centroids, without further splits
        
        last_centroids <- by(df_toKeep[, !(colnames(df_toKeep) %in% c('cluster', 'withinss')), drop = F],
                             df_toKeep[, 'cluster'],
                             function(x){sapply(x, mean)})
        
        last_centroids <- matrix(unlist(last_centroids), ncol = length(df_toKeep) - 2, byrow = T)
        
        km <- kmeans(df_toKeep[, !(colnames(df_toKeep) %in% c('cluster', 'withinss'))],
                     centers = last_centroids)
        
        df_toKeep$cluster <- km$cluster
        
        df_toKeep$withinss <- NULL
        
        
      }
      
    }
    
    # rescaling back variables (min_values and max_values have same variable names)
    
    df_toKeep[names(min_values)] <- 
      sapply(names(min_values), function(x){rescale(x = df_toKeep[x], from = c(0, 1), 
                                                    to = c(min_values[names(min_values) == x], 
                                                           max_values[names(max_values) == x]))})
    
    df_toKeep$cluster <- as.factor(df_toKeep$cluster)
    
    return(df_toKeep)
    
  } else {
    
    if (nrow(df) < 1000){
      
      # gower distance matrix + PAM
      
      gower_dist_matrix <- as.matrix(daisy(x = df, metric = 'gower'))
      
      # determining optimum number of clusters via silhouette width
      
      if (is.null(nclusters)){
        
        sil_width <- numeric()
        
        nclusters <- 10
        
        # too low row number
        
        nclusters <- ifelse(nclusters >= nrow(df), nrow(df) - 1, nclusters)
        
        for (i in 2:nclusters){
          
          pam_fit <- pam(x = gower_dist_matrix, diss = T, k = i)
          sil_width <- append(sil_width, pam_fit$silinfo$avg.width)
          
          print('siluette width:')
          print(sil_width)
          
          if (i >= 3){
            
            if ((sil_width[i - 1]/sil_width[i -2] < 1.1) | (sil_width[i - 1] < sil_width[i - 2])){
              
              nclusters <- which(sil_width == sil_width[i - 1])
              
              pam_fit <- pam(x = gower_dist_matrix, diss = T, k = nclusters)
              df$cluster <- pam_fit$clustering
              
              # rescaling back variables (min_values and max_values have same variable names)
              
              df[names(min_values)] <- 
                sapply(names(min_values), function(x){rescale(x = df[x], from = c(0, 1), 
                                                              to = c(min_values[names(min_values) == x], 
                                                                     max_values[names(max_values) == x]))})
              
              df$cluster <- as.factor(df$cluster)
              
              return(df)
              
            }
            
          }
          
        }
        
        # just in case situation above does not happen, run PAM with nclusters where max(silhouette width)
        
        nclusters <- which(sil_width == max(sil_width)) + 1
        
      }
      
      # too low row number
      
      nclusters <- ifelse(nclusters >= nrow(df), nrow(df) - 1, nclusters)

      pam_fit <- pam(x = gower_dist_matrix, diss = T, k = nclusters)
      df$cluster <- pam_fit$clustering
      
      # rescaling back variables (min_values and max_values have same variable names)
      
      df[names(min_values)] <- 
        sapply(names(min_values), function(x){rescale(x = df[x], from = c(0, 1), 
                                                      to = c(min_values[names(min_values) == x], 
                                                             max_values[names(max_values) == x]))})
      
      df$cluster <- as.factor(df$cluster)
      
      return(df)
      
      
    } else {
      
      # sample -> gower dist -> PAM -> gower distance of not sampled entrys to medoids of sample
      
      set.seed(888)
      
      sampled_rows <- sample.int(nrow(df), 1000)
      
      df_sampled <- data.frame(df[sampled_rows,], 
                               row.names = sampled_rows)
      colnames(df_sampled) <- colnames(df)
      
      df_not_sampled <- data.frame(df[!(row.names(df) %in% row.names(df_sampled)),], 
                                   row.names = rownames(df)[!(rownames(df) %in% rownames(df_sampled))])
      colnames(df_not_sampled) <- colnames(df)
      
      gower_dist_matrix <- as.matrix(daisy(x = df_sampled, metric = 'gower'))
      
      if (is.null(nclusters)){
        
        sil_width <- numeric()
        
        nclusters <- 10
        
        # too low rows number
        
        nclusters <- ifelse(nclusters >= nrow(df), nrow(df) - 1, nclusters)
        
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
              not_sampled_dist_medoids <- data.frame(t(gower.dist(data.x = df_sampled_medoids[, colnames(df_sampled_medoids) != 'cluster', drop = F],
                                                                  data.y = df_not_sampled)),
                                                     row.names = row.names(df_not_sampled))
              colnames(not_sampled_dist_medoids) <- row.names(df_sampled_medoids)
              closest_medoids_rownr <- apply(not_sampled_dist_medoids, 1, function(x){names(x)[x == min(x)]})
              
              df_not_sampled$cluster <- mapvalues(x = closest_medoids_rownr,
                                                  from = row.names(df_sampled_medoids),
                                                  to = df_sampled_medoids$cluster)
              
              df <- rbind(df_sampled, df_not_sampled)
              df <- df[order(as.integer(row.names(df))),]
              
              # rescaling back variables (min_values and max_values have same variable names)
              
              df[names(min_values)] <- 
                sapply(names(min_values), function(x){rescale(x = df[x], from = c(0, 1), 
                                                              to = c(min_values[names(min_values) == x], 
                                                                     max_values[names(max_values) == x]))})
              
              df$cluster <- as.factor(df$cluster)
              
              return(df)
              
            }
            
          }
          
        }
        
        # just in case situation above does not happen, run PAM with nclusters where max(silhouette width)
        
        nclusters <- which(sil_width == max(sil_width)) + 1
        
      }
      
      # too low rows number
      
      nclusters <- ifelse(nclusters >= nrow(df), nrow(df) - 1, nclusters)
      
      pam_fit <- pam(x = gower_dist_matrix, diss = T, k = nclusters)
      df_sampled$cluster <- pam_fit$clustering
      
      df_sampled_medoids <- df_sampled[pam_fit$medoids, ]
      not_sampled_dist_medoids <- data.frame(t(gower.dist(data.x = df_sampled_medoids[, colnames(df_sampled_medoids) != 'cluster', drop = F],
                                                          data.y = df_not_sampled)),
                                             row.names = row.names(df_not_sampled))
      colnames(not_sampled_dist_medoids) <- row.names(df_sampled_medoids)
      
      closest_medoids_rownr <- apply(not_sampled_dist_medoids, 1, function(x){names(x)[x == min(x)]})
      
      # in case an entry has same distance to more medoids, take first one and change to vector
      
      if (class(closest_medoids_rownr) == 'list'){
        
        closest_medoids_rownr <- sapply(closest_medoids_rownr, '[[', 1)
        
      }
      
      df_not_sampled$cluster <- mapvalues(x = closest_medoids_rownr,
                                          from = row.names(df_sampled_medoids),
                                          to = df_sampled_medoids$cluster)
      
      df <- rbind(df_sampled, df_not_sampled)
      
      df <- df[order(as.integer(row.names(df))),]
      
      # rescaling back variables (min_values and max_values have same variable names)
      
      df[names(min_values)] <- 
        sapply(names(min_values), function(x){rescale(x = df[x], from = c(0, 1), 
                                                      to = c(min_values[names(min_values) == x], 
                                                             max_values[names(max_values) == x]))})
      
      df$cluster <- as.factor(df$cluster)
      
      return(df)
      
    }
    
  }
  
}

