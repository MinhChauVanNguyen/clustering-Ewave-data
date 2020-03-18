library(reshape2)
library(ggplot2)
library(lattice)
library(viridisLite)
library(RColorBrewer)
library(grid)
library(vcd)
library(igraph)
library(dplyr)



NCols <- 5
NRows <- 20

## fixed values of starting parameters
  empty.vect <- c()
  fixed.par <- c(mu1=0,mu2=0.6,mu3=0.8,mu4=1.2,
                 phi1=0,phi2=0.5,phi3=0.7,phi4=0.3,
                 alpha1=0.5,alpha2=-0.5)
  fixed.par <- as.vector(fixed.par)
  vect <- c(empty.vect,fixed.par)
  

## calculate p 
## p is the response prob of y having category k=1,2,3,4 belongs to cluster r=1,2
  
  # OSM model
  sum.p.cluster1 <- exp(vect[1]-vect[5]*vect[9])+exp(vect[2]-vect[6]*vect[9])+exp(vect[3]-vect[7]*vect[9])+exp(vect[4]-vect[8]*vect[9])
  p1.clust1 <- exp(vect[1]-vect[5]*vect[9])/sum.p.cluster1
  p2.clust1 <- exp(vect[2]-vect[6]*vect[9])/sum.p.cluster1
  p3.clust1 <- exp(vect[3]-vect[7]*vect[9])/sum.p.cluster1
  p4.clust1 <- exp(vect[4]-vect[8]*vect[9])/sum.p.cluster1
  
  sum.p.cluster2 <- exp(vect[1]-vect[5]*vect[10])+exp(vect[2]-vect[6]*vect[10])+exp(vect[3]-vect[7]*vect[10])+exp(vect[4]-vect[8]*vect[10])
  p1.clust2 <- exp(vect[1]-vect[5]*vect[10])/sum.p.cluster2
  p2.clust2 <- exp(vect[2]-vect[6]*vect[10])/sum.p.cluster2
  p3.clust2 <- exp(vect[3]-vect[7]*vect[10])/sum.p.cluster2
  p4.clust2 <- exp(vect[4]-vect[8]*vect[10])/sum.p.cluster2
  
  theta1 <- c(p1.clust1, p2.clust1, p3.clust1, p4.clust1)
  theta2 <- c(p1.clust2, p2.clust2, p3.clust2, p4.clust2)
  
  
  cluster1 <- matrix(sample(1:4,size=NCols*NRows,replace=TRUE,prob=theta1),nrow=NRows)
  cluster2 <- matrix(sample(1:4,size=NCols*NRows,replace=TRUE,prob=theta2),nrow=NRows)
  
  ### OPTION 1 ###
  data5 <- rbind(cluster1,cluster2)
  my.data <- melt(data5)
  my.data <- data.frame(my.data)
  colnames(my.data) <- c("ROW","COL","Y")
  myplot <- ggplot(my.data, aes(x=COL, y=ROW, fill=factor(Y))) + geom_tile()
  myplot
  
  ### OPTION 2 ###
  my.data2 <- as.matrix(data5)
  myplot2 <- levelplot(t(my.data2[c(nrow(my.data2):1),]), col.regions = terrain.colors(100), 
                       xlab = "COL", ylab = "ROW")
  coul <- colorRampPalette(brewer.pal(8, "PiYG"))(25)
  myplot3 <- levelplot(my.data2, col.regions = coul, 
                       xlab = "ROW", ylab = "COL")
  myplot2
  myplot3
  
  ### OPTION 3 ###
  n <- 40
  NRows <- 20
  NCols <- 5
  phi.vect <- c(0, 0.5, 0.7, 1)
  R <- 2
  Labels <- list(categ = c("1","2","3","4"),
                 cluster = paste("R", 1:R, sep=""),
                 row = paste("r", 1:n, sep=""),
                 col = paste("c", 1:NCols, sep=""))
  data.matrix <- rbind(cluster1,cluster2)
  data.matrix
  is.matrix(data.matrix)
  ClusterRowY <- array(NA,n)
  for(i in 1:n){
    ClusterRowY[i] <- sample(1:R, 1, replace = TRUE)
  }
  rownames(ClusterRowY) <- Labels$row
  
  mosa <- spaced.mosaic.plot(y.mat = data.matrix,
                                v.phi.est = phi.vect,
                                R.best = R,
                                ClusterRowY = ClusterRowY,
                                labels = Labels)
  mosa
  is.matrix(mosa)
  names(mosa)
  plot(mosa, main = "Row Clustering Results")
  dim(mosa)
  dimnames(mosa) <- list
  xyplot(Labels~value, data=mosa)
  
  
  
  ### OPTION 4 ###
  my.plot5 <- ggplot(my.data, aes(x = COL, y = ROW)) + 
    geom_raster(aes(fill=Y)) + 
    scale_fill_gradient(low="grey90", high="red") +
    labs(x="letters", y="LETTERS", title="Matrix") +
    theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                       axis.text.y=element_text(size=9),
                       plot.title=element_text(size=11))
  my.plot5
  
  ### OPTION 5 ###
  cluster1 <- matrix(sample(1:4,size=NCols*NRows,replace=TRUE,prob=theta1),nrow=NRows)
  cluster2 <- matrix(sample(1:4,size=NCols*NRows,replace=TRUE,prob=theta2),nrow=NRows)
  
  data5 <- rbind(cluster1,cluster2)
  is.matrix(data5)
  clusters <- lapply(split(data5, rep(1:2, each = nrow(data5)/2)),
         function(a) matrix(a, ncol = ncol(data5)))
  cluster1 <- clusters[[1]]
  cluster2 <- clusters[[2]]
  data5$colour <- ifelse(cluster1==cluster2, cluster1, 0)
  data5$colour <- factor(data5$colour)
  data5$ROW <- factor(data5$ROW, levels=unique(cluster1[,1]))
  data5$COL <- factor(data5$COL, levels=unique(cluster2[,2]))
  
  
  ggplot(data5, aes(x = COL, y = ROW, fill=colour)) + 
    geom_raster() + 
    scale_fill_manual(values=c("grey80", "#B40404", "#0B6121", "#FFBF00")) +
    labs(x="letters", y="LETTERS", title="Matrix") +
    theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                       axis.text.y=element_text(size=9),
                       plot.title=element_text(size=11),
                       legend.text=element_text(size=7))
  
  

  
  ############################## BOX PLOT ######################################
  
  ## true value
  true.par <- c(mu1=0,mu2=0.6,mu3=0.8,mu4=1.2,u1=0,u2=0.1,u3=0.2,u4=1,alpha1=0.5,alpha2=-0.5)
  
  
  ## generated value
  gen.par <- rnorm(n=length(true.par))
  gen.par
  
  first <- boxplot(true.par,
                   main = "Range of true parameter values",
                   xlab = "Range",
                   ylab = "Fixed Value",
                   col = "orange",
                   border = "brown",
                   horizontal = TRUE,
                   notch = FALSE)
  
  second <- boxplot(gen.par,
                    main = "Range of generated parameter values",
                    xlab = "Range",
                    ylab = "Generated Value",
                    col = "green",
                    border = "blue",
                    horizontal = TRUE,
                    notch = FALSE)
  
  third <- boxplot(true.par,gen.par,
                   main = "Generated vs. Fixed values",
                   names = c("Fixed", "Rnorm"),
                   las = 2,
                   col = c("green","yellow"),
                   border = "blue",
                   horizontal = TRUE,
                   notch = FALSE)
  
  gen.mu.vect <- rnorm(n=100)
  fixed.mu.vect <- c(0,0.6,0.8,1)
  hist(gen.mu.vect)
  abline(v=fixed.mu.vect, col="red")
  
  
