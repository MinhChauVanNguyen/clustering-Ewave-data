library(clustord)
library(reshape2)

# True parameters
fixed.par <- c(mu1=0, mu2=0.6, mu3=0.8, mu4=1.2,
               phi1=0, phi2=0.4, phi3=0.6, phi4=1,
               alpha1=0.5, alpha2=-0.5)


################################ ELBOW PLOT #################################
data.row <- function(NRows, NCols, fixed.par, model="OSM"){
  ## fixed values of starting parameters
  empty.vect <- c()
  fixed.par <- as.vector(fixed.par)
  vect <- c(empty.vect, fixed.par)
  
  ## calculate p 
  ## p is the response prob of y having category k=1,2,3,4 belongs to cluster r=1,2
  
  # OSM model
  sum.p.cluster1 <- exp(vect[1]-vect[5]*vect[9])+
    exp(vect[2]-vect[6]*vect[9])+
    exp(vect[3]-vect[7]*vect[9])+
    exp(vect[4]-vect[8]*vect[9])
  p1.clust1 <- exp(vect[1]-vect[5]*vect[9])/sum.p.cluster1
  p2.clust1 <- exp(vect[2]-vect[6]*vect[9])/sum.p.cluster1
  p3.clust1 <- exp(vect[3]-vect[7]*vect[9])/sum.p.cluster1
  p4.clust1 <- exp(vect[4]-vect[8]*vect[9])/sum.p.cluster1
  
  sum.p.cluster2 <- exp(vect[1]-vect[5]*vect[10])+
    exp(vect[2]-vect[6]*vect[10])+
    exp(vect[3]-vect[7]*vect[10])+
    exp(vect[4]-vect[8]*vect[10])
  p1.clust2 <- exp(vect[1]-vect[5]*vect[10])/sum.p.cluster2
  p2.clust2 <- exp(vect[2]-vect[6]*vect[10])/sum.p.cluster2
  p3.clust2 <- exp(vect[3]-vect[7]*vect[10])/sum.p.cluster2
  p4.clust2 <- exp(vect[4]-vect[8]*vect[10])/sum.p.cluster2
  
  theta1 <- c(p1.clust1, p2.clust1, p3.clust1, p4.clust1)
  theta2 <- c(p1.clust2, p2.clust2, p3.clust2, p4.clust2)

  cluster1 <- matrix(sample(1:4,size=NCols*NRows,replace=TRUE,prob=theta1),nrow=NRows)
  cluster2 <- matrix(sample(1:4,size=NCols*NRows,replace=TRUE,prob=theta2),nrow=NRows)
 
  
  ## put cluster 1 & 2 together to produce a data frame with three factor columns
  data <- rbind(cluster1,cluster2)
  pseudo.data <- melt(data)
  colnames(pseudo.data) <- c("ROW","COL","Y")
  pseudo.data$Y <- factor(pseudo.data$Y)
  return(data.frame(pseudo.data))
}


AIC.score.func <- function(N, npar, NRows=100, NCols=10, formula="Y~row", model="OSM", nclus.row){
  
  initvect <- matrix(NA, nrow=N, ncol=npar)
  cluster.result <- list()
  
  # generate pseudo-data 
  data <- data.row(NRows=NRows, NCols=NCols, fixed.par=fixed.par, model=model)
  # initialize maximum
  Max.AIC <- -Inf                                   
  
  for(i in 1:N){
    initvect[i,] <- gen.row(npar=npar, seed=i)
    cluster <- rowclustering(formula = formula, 
                             model = model, 
                             nclus.row = nclus.row, 
                             long.df = data, 
                             initvect = initvect[i,])
    
    cluster.result[[i]] <- list()
    cluster.result[[i]][[1]] <- cluster$RowClusters
    cluster.result[[i]][[2]] <- cluster$EM.status$best.lli
    cluster.result[[i]][[3]] <- cluster$criteria
    
    # name the required elements for each list
    names(cluster.result[[i]]) <- c("Row Clusters", "best.lli", "criteria")

    # return the maximum value of EM.status$best.lli of N iterations
    AIC.score <- cluster.result[[i]]$criteria$AIC
    if (AIC.score > Max.AIC){
      Max.AIC <- AIC.score
    }
  }
  return(Max.AIC)
}

aic2 <- AIC.score.func(N=100,npar=6,nclus.row=2)
aic3 <- AIC.score.func(N=100,npar=7,nclus.row=3)
aic4 <- AIC.score.func(N=100,npar=8,nclus.row=4)
aic5 <- AIC.score.func(N=100,npar=9,nclus.row=5)
aic6 <- AIC.score.func(N=100,npar=10,nclus.row=6)
aic7 <- AIC.score.func(N=100,npar=11,nclus.row=7)
aic8 <- AIC.score.func(N=100,npar=12,nclus.row=8)


x <- c(2,3,4,5,6,7,8)
y <- c(aic2,aic3,aic4,aic5,aic6,aic7,aic8)
plot(x,y, xlim = range(x), ylim = range(y), xlab = "Number of Row Clusters", 
     ylab = "AIC score", xaxt = "n", 
     main = "Elbow Plot for pseudo-data of 2 clusters",
     pch = 8, cex=2, lwd=2)
axis(1, at = 2:8)

lines(x[order(x)], y[order(x)], xlim=range(x), ylim=range(y), pch=16)
symbols(x=2,y=aic2, circles=0.2, add=T, inches=F, fg="red")


#####################################################################################

dat.row <- function(NRows, NCols, fixed.par, model="OSM"){
  ## fixed values of starting parameters
  empty.vect <- c()
  fixed.par <- as.vector(fixed.par)
  vect <- c(empty.vect, fixed.par)
  
  # OSM model
  sum.p.cluster1 <- exp(vect[1]-vect[5]*vect[9])+
    exp(vect[2]-vect[6]*vect[9])+
    exp(vect[3]-vect[7]*vect[9])+
    exp(vect[4]-vect[8]*vect[9])
  p1.clust1 <- exp(vect[1]-vect[5]*vect[9])/sum.p.cluster1
  p2.clust1 <- exp(vect[2]-vect[6]*vect[9])/sum.p.cluster1
  p3.clust1 <- exp(vect[3]-vect[7]*vect[9])/sum.p.cluster1
  p4.clust1 <- exp(vect[4]-vect[8]*vect[9])/sum.p.cluster1
  
  sum.p.cluster2 <- exp(vect[1]-vect[5]*vect[10])+
    exp(vect[2]-vect[6]*vect[10])+
    exp(vect[3]-vect[7]*vect[10])+
    exp(vect[4]-vect[8]*vect[10])
  p1.clust2 <- exp(vect[1]-vect[5]*vect[10])/sum.p.cluster2
  p2.clust2 <- exp(vect[2]-vect[6]*vect[10])/sum.p.cluster2
  p3.clust2 <- exp(vect[3]-vect[7]*vect[10])/sum.p.cluster2
  p4.clust2 <- exp(vect[4]-vect[8]*vect[10])/sum.p.cluster2
  
  sum.p.cluster3 <- exp(vect[1]-vect[5]*vect[11])+
    exp(vect[2]-vect[6]*vect[11])+
    exp(vect[3]-vect[7]*vect[11])+
    exp(vect[4]-vect[8]*vect[11])
  p1.clust3 <- exp(vect[1]-vect[5]*vect[11])/sum.p.cluster2
  p2.clust3 <- exp(vect[2]-vect[6]*vect[11])/sum.p.cluster2
  p3.clust3 <- exp(vect[3]-vect[7]*vect[11])/sum.p.cluster2
  p4.clust3 <- exp(vect[4]-vect[8]*vect[11])/sum.p.cluster2
  
  sum.p.cluster4 <- exp(vect[1]-vect[5]*vect[12])+
    exp(vect[2]-vect[6]*vect[12])+
    exp(vect[3]-vect[7]*vect[12])+exp(vect[4]-vect[8]*vect[12])
  p1.clust4 <- exp(vect[1]-vect[5]*vect[12])/sum.p.cluster2
  p2.clust4 <- exp(vect[2]-vect[6]*vect[12])/sum.p.cluster2
  p3.clust4 <- exp(vect[3]-vect[7]*vect[12])/sum.p.cluster2
  p4.clust4 <- exp(vect[4]-vect[8]*vect[12])/sum.p.cluster2
  
  theta1 <- c(p1.clust1, p2.clust1, p3.clust1, p4.clust1)
  theta2 <- c(p1.clust2, p2.clust2, p3.clust2, p4.clust2)
  theta3 <- c(p1.clust3, p2.clust3, p3.clust3, p4.clust3)
  theta4 <- c(p1.clust4, p2.clust4, p3.clust4, p4.clust4)
  
  cluster1 <- matrix(sample(1:4,size=NCols*NRows,replace=TRUE,prob=theta1),nrow=NRows)
  cluster2 <- matrix(sample(1:4,size=NCols*NRows,replace=TRUE,prob=theta2),nrow=NRows)
  cluster3 <- matrix(sample(1:4,size=NCols*NRows,replace=TRUE,prob=theta3),nrow=NRows)
  cluster4 <- matrix(sample(1:4,size=NCols*NRows,replace=TRUE,prob=theta4),nrow=NRows)
  
  ## put cluster 1 & 2 together to produce a data frame with three factor columns
  data <- rbind(cluster1,cluster2,cluster3,cluster4)
  pseudo.data <- melt(data)
  colnames(pseudo.data) <- c("ROW","COL","Y")
  pseudo.data$Y <- factor(pseudo.data$Y)
  return(data.frame(pseudo.data))
}


fixed.par2 <- c(mu1=0,mu2=0.6,mu3=0.8,mu4=1.2,
               phi1=0,phi2=0.4,phi3=0.6,phi4=1,
               alpha1=0.5,alpha2=-0.5,alpha3=0.25,alpha4=-0.25)


AIC.score.func2 <- function(N, npar, NRows=100, NCols=20, formula="Y~row", model="OSM", nclus.row){
  
  initvect <- matrix(NA, nrow=N, ncol=npar)
  cluster.result <- list()
  
  # generate pseudo-data 
  data <- dat.row(NRows=NRows, NCols=NCols, fixed.par=fixed.par2, model=model)
  # initialize maximum
  Max.AIC <- -Inf                                   
  
  for(i in 1:N){
    initvect[i,] <- gen.row(npar=npar, seed=i)
    cluster <- rowclustering(formula = formula, 
                             model = model, 
                             nclus.row = nclus.row, 
                             long.df = data, 
                             initvect = initvect[i,])
    
    cluster.result[[i]] <- list()
    cluster.result[[i]][[1]] <- cluster$RowClusters
    cluster.result[[i]][[2]] <- cluster$EM.status$best.lli
    cluster.result[[i]][[3]] <- cluster$criteria
    
    # name the required elements for each list
    names(cluster.result[[i]]) <- c("Row Clusters", "best.lli", "criteria")
    
    # return the maximum value of EM.status$best.lli of N iterations
    AIC.score <- cluster.result[[i]]$criteria$AIC
    if (AIC.score > Max.AIC){
      Max.AIC <- AIC.score
    }
  }
  return(Max.AIC)
}

aic22 <- AIC.score.func2(N=100,npar=6,nclus.row=2)
aic33 <- AIC.score.func2(N=100,npar=7,nclus.row=3)
aic44 <- AIC.score.func2(N=100,npar=8,nclus.row=4)
aic55 <- AIC.score.func2(N=100,npar=9,nclus.row=5)
aic66 <- AIC.score.func2(N=100,npar=10,nclus.row=6)
aic77 <- AIC.score.func2(N=100,npar=11,nclus.row=7)
aic88 <- AIC.score.func2(N=100,npar=12,nclus.row=8)

x3 <- c(2,3,4,5,6,7,8)
y3 <- c(aic22,aic33,aic44,aic55,aic66,aic77,aic88)
plot(x3,y3,xlim = range(x3), ylim = range(y3), xlab = "Number of Row Clusters", 
     ylab = "AIC score", xaxt = "n", 
     main = "Elbow Plot for pseudo-data of 4 clusters",
     pch = 8, cex=2, lwd=2)
axis(1, at = 2:8)

lines(x3[order(x3)], y3[order(x3)], xlim=range(x3), ylim=range(y3), pch=16)
symbols(x=2,y=aic22, circles=0.2, add=T, inches=F, fg="red")



