
library(ggplot2)
library(colorspace)
library(scales)
library(fpc)

data("iris")
str(iris)

iris <- cbind(iris, sapply(iris[1:4], rescale, to = c(0,1)))
colnames(iris)[6:9] <- c('Sepal.Length.Scales', 'Sepal.Width.Scaled',
                         'Petal.Length.Scaled', 'Petal.Width.Scaled')
iris_test <- iris[6:9]

kmeansHH <- function(x, nclusters = NULL){
  
  df <- x
  nvars <- length(df)
  optimal_nclusters <- F
  nclusters <- ifelse(nclusters == 1, 2, nclusters)
  kms <- list()

  if (is.null(nclusters)){
    
    optimal_nclusters <- T
    nclusters <- 25
    calHar <- numeric()
    
  }
  
  for (i in 1:(nclusters - 1)){
    
    print(paste('run:', i))
    print(df)
    
    variance <- sapply(df, var) # variables variance
    print('variables variance')
    print(variance)
    
    max_var <- variance[variance == max(variance)] # variable with maximum variance
    print('maximum variance and variable')
    print(max_var)
    
    max_var_mean <- mean(df[[names(max_var)]]) # mean of this variable
    print('mean of this variable')
    print(max_var_mean)
    
    cent1 <- sapply(df[df[names(max_var)] >= 0 &
                                df[names(max_var)] < max_var_mean,], mean)
    print('cent1')
    print(cent1)
    
    cent2 <- sapply(df[df[names(max_var)] >= max_var_mean,], mean) # based on its value, calculate 1st and 2nd centroid
    print('cent2')
    print(cent2)
    
    if (i == 1){
      
      cent <- rbind(cent1, cent2)
      
    } else {
      
      cent <- rbind(cent, cent1)
      cent <- rbind(cent, cent2)
      
    }
    
    print(cent)
    kms <- append(kms, list(kmeans(x = x, centers = cent, algorithm = 'Lloyd')))
    print('unique clusters')
    print(unique(kms[[i]]$cluster))
    
    if (optimal_nclusters == T){

      calHar <- append(calHar, calinhara(x, kms[[i]]$cluster))
      
      print('calinski-harabasz index:')
      print(calHar)
      
      if (i >= 2){
        
        if (calHar[i] < calHar[i - 1]){
          
          km <- kms[[i -1]]
          
          return(km)
          
        }

        
      }
      
    }
    
    if (i != (nclusters - 1)){
      
      print('clusters withinss')
      print(kms[[i]]$withinss)
    
      high_var_cluster <- which(kms[[i]]$withinss == max(kms[[i]]$withinss)) # cluster with highest withinss
      print(paste('cluster nr. with highest withinss:', high_var_cluster))
      
      low_var_cluster <- which(kms[[i]]$withinss != max(kms[[i]]$withinss)) # remaining clusters
      print('clusters with not maximum withinss')
      print(low_var_cluster)
      
      
      cent <- kms[[i]]$centers[low_var_cluster,] # keep only centroids of clusters with not highest withinss
      print('centroids to remain')
      print(cent)
      
      print(kms[[i]]$cluster == high_var_cluster)
      df <- x[kms[[i]]$cluster == high_var_cluster,] # keep only data of cluster with highest withinss for next iteration
      
      print('% of data of cluster with maximum withinss')
      print(nrow(df)/nrow(x))
      print(str(df))
      
    }
    
    print('-----------')
    
  }
  
  km <- kms[[nclusters - 1]]
  
  return(km)
  
}

# rozdelit dataset na 2 casti podla max_var
# vyratat centroidy tychto 2 casti
# spustit k-means s tymito centroidmi
# vybrat cluster, ktory ma vacsiu within-cluster variance
# opakovat, pokial nedosiahnem specifikovany pocet clustrov, napr. 4
