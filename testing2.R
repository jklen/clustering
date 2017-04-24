
library(ggplot2)


kmeansHH2 <- function(x, nclusters){
  
  df <- x
  df_toKeep <- data.frame()
  
  for (i in (1:(nclusters))){
    
    print(paste('run: ', i))
    print(paste('row number df_toKeep: ', nrow(df_toKeep)))
    
    if (i == 1){
    
      variance <- sapply(df, var) # variables variance
      max_var <- variance[variance == max(variance)] # variable with maximum variance
    
    }
    
    df$cluster <- NULL
    
    max_var_mean <- mean(df[[names(max_var)]]) # mean of this variable
    
    cent1 <- sapply(df[df[names(max_var)] >= 0 &
                         df[names(max_var)] < max_var_mean,], mean)
    cent2 <- sapply(df[df[names(max_var)] >= max_var_mean,], mean)
    
    cent <- rbind(cent1, cent2)
    
    print(colnames(df))
    
    km <- kmeans(x = df, centers = cent, algorithm = 'Lloyd')
    
    df$cluster <- km$cluster
    
    #print(paste('df_toKeep rows number:', nrow(df_toKeep)))
    
    df1 <- df[df$cluster == 1,]
    print(paste('df1 rows number: ', nrow(df1)))
    
    df2 <- df[df$cluster == 2,]
    print(paste('df2 rows number: ', nrow(df2)))
    
    if (nrow(df1) == 1 | nrow(df2) == 1){
      
      df_last <- rbind(df1, df2)
      df_last$cluster <- i
      df_toKeep <- rbind(df_toKeep, df_last)
      
      return (df_toKeep)

    }
    
    variance1 <- sapply(df1, var)
    print('variable variance 1')
    print(variance1)
    
    variance2 <- sapply(df2, var)
    print('variable variance 2')
    print(variance2)
    
    max_var1 <- variance1[variance1 == max(variance1)]
    print('maxvar1')
    print(max_var1)
    
    max_var2 <- variance2[variance2 == max(variance2)]
    print('maxvar2')
    print(max_var2)
    
    if ((max_var1 > max_var2) & i != nclusters){
      
      df <- df1
      max_var <- max_var1
      df2$cluster <- i
      df_toKeep <- rbind(df_toKeep, df2)

    } else {
      
      if ((max_var1 < max_var2) & i != nclusters){
        
        df <- df2
        max_var <- max_var2
        df1$cluster <- i
        df_toKeep <- rbind(df_toKeep, df1)
        
      } else {
        
        df_last <- rbind(df1, df2)
        df_last$cluster <- i
        df_toKeep <- rbind(df_toKeep, df_last)
        
      }
      
    }
    
    print(paste('rows number df_toKeep: ', nrow(df_toKeep)))
    
    df_toKeep$cluster <- as.character(df_toKeep$cluster)
    
    p <- ggplot(aes(x = Petal.Length.Scaled, y = Petal.Width.Scaled), data = df_toKeep) +
      geom_point(aes(color = cluster)) +
      ggtitle(i) +
      ylim(0, 1.1) +
      xlim(0, 1.1)
    print(p)
    
    print('--------------------------------------------')
    
  }
  
  #df_toKeep$cluster <- as.factor(df_toKeep$cluster)
  
  return(df_toKeep)
  
}