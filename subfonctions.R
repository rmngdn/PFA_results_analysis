# Subfonctions for main_clustering.R:


#libraries:

require(SNFtool)
library(SNFtool)
library(ggplot2)
library(cluster)

# SNFCalc <- function(cancerTypeName = "GLIO_", path.to.data ="~/Documents/Local_Work/GBM/", n, Knn, alpha, t) {
#   # K is the number of nearest neighboors: N/10 is a good value 
#   # alpha is the hyperparameter ~  pick it in 0.3 and 0.8 but the method isn't really sensitive to it
#   # t is the number of iterations for the dispersion process
#   
#   
#   listName <- c("Gene","Mirna","Methy")
#   
#   #Import data
#   listName <- unlist(lapply(listName, function(p) {return(paste0(cancerTypeName, p, "_Expression.txt"))}))
#   datasets <- lapply(listName, function(name) {return(read.table(paste0(path.to.data, name), header = TRUE))})
#   
#   #Rename patients because patient IDs don't match between different data types in the GBM data set (and other test datasets) ...
#   # We can suppose that it's ordered however and just rename by patient_1, patient_2,...patient_215 
#   # in both SNF and PCA to enable the comparison and make sure we have the same order (to compare partition results)
#   
#   pID = unlist(lapply((1:n), function(p) {return(paste0("patient_",p))})) 
#   A <- function(data, pID) {colnames(data) <- pID;return(data)}
#   datasets <- lapply(datasets, function(data) A(data, pID))
#   
#   # SNF
#   
#   dist_list <- lapply(datasets, function(p) dist2(as.matrix(t(p)),as.matrix(t(p))))
#   affinityList <- lapply(dist_list, function(dist) affinityMatrix(dist, Knn, alpha ) )
#   fusedAffinity <- SNF(affinityList, Knn, t)
#   return(fusedAffinity)
# }

SNFCalc2 <- function(datasets, Knn, alpha, t) {
  # K is the number of nearest neighboors: N/10 is a good value 
  # alpha is the hyperparameter ~  pick it in 0.3 and 0.8 but the method isn't really sensitive to it
  # t is the number of iterations for the dispersion process
  
  # SNF
  dist_list <- lapply(datasets, function(p) dist2(as.matrix(t(p)),as.matrix(t(p))))
  affinityList <- lapply(dist_list, function(dist) affinityMatrix(dist, Knn, alpha ) )
  fusedAffinity <- SNF(affinityList, Knn, t)
  return(fusedAffinity)
}



plotPC <- function(Y_PFA,l){
  plot(setNames(data.frame(Y_PFA[,l]),(1:(dim(Y_PFA)[2]))[l]))
}



plotDistMatrix <- function(dist, cluster, return.index = FALSE) {
  library(plyr)
  library(SNFtool)
  n <- dim(dist)[1]
  occ <- count(df = data.frame("c" = cluster),vars =  "c")$freq
  C <- length(occ)
  counter <- matrix(1,1,C)
  starter <- matrix(0,1,C)
  for (x in 2:C) {
    starter[1,x] <- sum(occ[1:x-1]) 
  }
  index <- vector('integer', n)
  for (i in 1:n) {
    clustNumb <- cluster[i]
    index[i] <- starter[1,clustNumb] + counter[1,clustNumb] 
    counter[clustNumb] <- counter[clustNumb] + 1
  }
  dist <- dist[index,index]
  image(dist, axes = FALSE)
  title(main = "Distance matrix ordered by clusters", xlab = "Patients", ylab = "Patients")
  
  #heat <- heatmap(x = dist,symm = TRUE, distfun = daisy, verbose = FALSE )
  #return(heat)
}



spectral <- function(x, K) {
  return(data.frame("cluster" = spectralClustering(W, K)))
}



sd_mean <- function(dataList) {
  for (dataType in dataList) {
    # dataType is feature x sample
    n <- dim(dataType)[2]
    h <- dim(dataType)[1]
    means <- data.frame('mean' = 1:h)
    sds <- data.frame('sd' = 1:h)
    for (i in 1:h) {
      feat <- as.numeric(dataType[i,])
      means$mean[i] <- mean(feat) 
      sds$sd[i] <- sd(feat)
    }
    plot.means <- ggplot(means) + aes(x="", y=mean) + geom_boxplot()
    plot.sds <- ggplot(sds) + aes(x="", y = sd) + geom_boxplot()
    par(mfrow=c(2,2))
    plot(plot.means, main = "Boxplot of means of features")
    plot(plot.sds, main = "Boxplot of standard deviations of features")
  }
}

clusterQuality <- function(Y, Klist, distanceMatrix = FALSE) {
  # Y have to be in samples x features format
  i <- 1
  width <- vector('numeric', length(Klist))
  for (k in Klist) {
    clarax <- clara(x = Y, k = k)
    sil <- silhouette(clarax, full = TRUE)
    plot(sil, main = paste0("Silhouette plot of clusters with clara function with ", k," clusters."))
    width[i] <- summary(sil)$avg.width
    i <- i + 1
  }
  plot(width, main = "Total average width of clusters for each number of cluster", xlab = "Number of clusters", ylab = "Total average width")
}

pam.clusterQuality <- function(Y, Klist) {
  # Y have to be in samples x features format
  dist <- daisy(Y)
  i <- 1
  width <- data.frame("K" = 1:length(Klist), "width" = 1:length(Klist))
  for (k in Klist) {
    pamx <- pam(x = dist, k = k, diss = TRUE)
    sil <- silhouette(x = pamx$clustering, dist = dist)
    plot(sil, main = paste0("Silhouette plot of clusters with clara function with ", k," clusters."))
    width$width[i] <- summary(sil)$avg.width
    width$K[i] <- k
    i <- i + 1
  }
  plot(width)#, main = "Total average width of clusters for each number of cluster", xlab = "Number of clusters", ylab = "Total average width")
}










