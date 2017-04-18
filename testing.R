
library(ggplot2)
library(colorspace)
library(scales)

data("iris")
str(iris)

iris <- cbind(iris, sapply(iris[1:4], rescale, to = c(0,1)))
colnames(iris)[6:9] <- c('Sepal.Length.Scales', 'Sepal.Width.Scaled',
                         'Petal.Length.Scaled', 'Petal.Width.Scaled')
iris_test <- iris[6:9]

kmeansHH <- function(x, nclusters){
  
  df <- x
  nvars <- length(df)
  
  for (i in 1:(nclusters - 1)){
    
    print(paste('run:', i))
    print(df)
    
    variance <- sapply(df, var) # vyratm rozptyl premennych
    print('variables variance')
    print(variance)
    
    max_var <- variance[variance == max(variance)] # ziskam maximalny rozptyl a premennu
    print('maximum variance and variable')
    print(max_var)
    
    max_var_mean <- mean(df[[names(max_var)]]) # ziskam priemer tejto premennej
    print('mean of this variable')
    print(max_var_mean)
    
    cent1 <- sapply(df[df[names(max_var)] >= 0 &
                                df[names(max_var)] < max_var_mean,], mean)
    print('cent1')
    print(cent1)
    
    cent2 <- sapply(df[df[names(max_var)] >= max_var_mean,], mean) # podla jej priemeru vyratam 1. a 2. centroid
    print('cent2')
    print(cent2)
    
    if (i == 1){
      
      cent <- matrix(rbind(cent1, cent2), ncol = nvars)
      
    } else {
      
      cent <- matrix(rbind(cent, cent1), ncol = nvars)
      cent <- matrix(rbind(cent, cent2), ncol = nvars)
      
    }
    
    print(cent)
    km <- kmeans(x = x, centers = cent, algorithm = 'Lloyd')
    print('unique clusters')
    print(unique(km$cluster))
    
    if (i != (nclusters - 1)){
      
      print('clusters withinss')
      print(km$withinss)
    
      high_var_cluster <- which(km$withinss == max(km$withinss)) # index (cislo) clustra s najvacsim withinss
      print(paste('cluster nr. with highest withinss:', high_var_cluster))
      
      low_var_cluster <- which(km$withinss != max(km$withinss)) # index (cisla) ostatnych clustrov
      print('clusters with not maximum withinss')
      print(low_var_cluster)
      
      
      cent <- km$centers[low_var_cluster,] # ponecham len centroidy s mensou ako maximalnou withinss
      print('centroids to remain')
      print(cent)
      
      print(km$cluster == high_var_cluster)
      df <- x[km$cluster == high_var_cluster,] # ponecham data len clustru s maximalnou withinss
      
      print('% of data of cluster with maximum withinss')
      print(nrow(df)/nrow(x))
      print(str(df))
      
    }
    
    print('-----------')
    
  }
  
  return(km)
  
}

# rozdelit dataset na 2 casti podla max_var
# vyratat centroidy tychto 2 casti
# spustit k-means s tymito centroidmi
# vybrat cluster, ktory ma vacsiu within-cluster variance
# opakovat, pokial nedosiahnem specifikovany pocet clustrov, napr. 4
