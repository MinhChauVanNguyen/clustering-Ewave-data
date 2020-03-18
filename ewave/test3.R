#Set up a simulation study to test how reliably, 
#in a diverse range of scenarios (different starting values combinations), 
#we were able to estimate the parameters of stereotype models 
#using the EM algorithm. 


N <- 5
#npar <- 4
npar <- 6
npar <- 2
formula <- "Y~row"
model <- "Binary"
#model <- "POM"
nclus.row <- 2
n.r <- 20
NCols <- 10

library(clustord)
library(reshape2)

gen.row <- function(npar,seed){
  set.seed(seed)
  initvect <- rnorm(n=npar)
  initvect
}

data.row <- function(n.r, NCols, fixed.par, model){
  ## fixed values of starting parameters
  empty.vect <- c()
  fixed.par <- as.vector(fixed.par)
  vect <- c(empty.vect, fixed.par)
  
  if (!is.character(model) || !is.vector(model) || length(model)!=1) 
    stop("model must be a string, 'OSM' or 'POM'.")
  if (!(model %in% c("OSM","POM","Binary"))) 
    stop("model must be either 'OSM' or POM' or 'Binary'.")
  ## calculate p 
  ## p is the response prob of y having category k=1,2,3,4 belongs to cluster r=1,2
  
  if(model=="OSM"){
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
    
    cluster1 <- matrix(sample(1:4,size=NCols*n.r,replace=TRUE,prob=theta1),nrow=n.r)
    cluster2 <- matrix(sample(1:4,size=NCols*n.r,replace=TRUE,prob=theta2),nrow=n.r)
    
  }else if(model=="POM"){
    p1 <- expit(vect[1]-vect[9])
    p2 <- expit(vect[2]-vect[9]) - p1
    p3 <- -expit(p2) + expit(p1)
    p4 <- 1 - (p1 + p2 + p3)
    theta.PO <- c(p1, p2, p3, p4)
    
    cluster1 <- matrix(sample(1:4,size=NCols*n.r,replace=TRUE,prob=theta.PO),nrow=n.r)
    cluster2 <- matrix(sample(1:4,size=NCols*n.r,replace=TRUE,prob=theta.PO),nrow=n.r)
    
  }else{
    # Binary model
    theta1.Bin <- expit(vect[2]+vect[9])
    theta0.Bin <- 1 - theta1.Bin
    
    cluster1 <- matrix(sample(c(0,1),size=NCols*n.r,replace=TRUE,prob=c(theta0.Bin,theta1.Bin)),nrow=n.r)
    cluster2 <- matrix(sample(c(0,1),size=NCols*n.r,replace=TRUE,prob=c(theta0.Bin,theta1.Bin)),nrow=n.r)
  }
  
  ## put cluster 1 & 2 together to produce a data frame with three factor columns
  data <- rbind(cluster1,cluster2)
  pseudo.data <- melt(data)
  colnames(pseudo.data) <- c("ROW","COL","Y")
  pseudo.data$Y <- factor(pseudo.data$Y)
  return(data.frame(pseudo.data))
}


fixed.par <- c(mu1=0, mu2=0.6, mu3=0.8, mu4=1.2,
               phi1=0, phi2=0.4, phi3=0.6, phi4=1,
               alpha1=0.5, alpha2=-0.5)


N <- 5
est.func <- function(N, npar, n.r, NCols, formula="Y~row", model, nclus.row){
  initvect <- matrix(NA, nrow=N, ncol=npar)
  cluster.results <- list()
  alpha.est <- c()
  
  # generate pseudo-data 
  data <- data.row(n.r=n.r, NCols=NCols, fixed.par=fixed.par, model=model)
  
  for(i in 1:N){
      initvect[i,] <- gen.row(npar=npar, seed=i)
      cluster <- rowclustering(formula = formula, 
                               model = model, 
                               nclus.row = nclus.row, 
                               long.df = data, 
                               initvect = initvect[i,])
  }  
      cluster.results[[i]] <- list()
      for(r in 1:nclus.row){
        if(model=="OSM"){
          cluster.results[[i]][[1]] <- cluster$parlist.out$mu
          cluster.results[[i]][[2]] <- cluster$parlist.out$phi
          cluster.results[[i]][[3]] <- cluster$parlist.out$alpha
          
          mu2.est <- mean(sapply(cluster.results, function(res) res[[1]][2]))
          mu3.est <- mean(sapply(cluster.results, function(res) res[[1]][3]))
          mu4.est <- mean(sapply(cluster.results, function(res) res[[1]][4]))
          
          phi2.est <- mean(sapply(cluster.results, function(res) res[[2]][2]))
          phi3.est <- mean(sapply(cluster.results, function(res) res[[2]][3]))
          
          alpha.est[r] <- mean(sapply(cluster.results, function(res) res[[3]][r]))
          result <- c(mu2.est,mu3.est,mu4.est,phi2.est,phi3.est,alpha.est)
          
        }else if(model=="POM"){
          cluster.results[[i]][[1]] <- cluster$parlist.out$mu
          cluster.results[[i]][[3]] <- cluster$parlist.out$alpha
          
          mu2.est <- mean(sapply(cluster.results, function(res) res[[1]][1]))
          mu3.est <- mean(sapply(cluster.results, function(res) res[[1]][2]))
          mu4.est <- mean(sapply(cluster.results, function(res) res[[1]][3]))
          alpha.est[r] <- mean(sapply(cluster.results, function(res) res[[3]][r]))
          
          result <- c(mu2.est,mu3.est,mu4.est, alpha.est)
        }else{
          cluster.results[[i]][[1]] <- cluster$parlist.out$mu
          cluster.results[[i]][[3]] <- cluster$parlist.out$alpha
        
          mu.vect <- unlist(cluster.results)
          mu.est <- mean(mu.vect)
          alpha.est[r] <- mean(sapply(cluster.results, function(res) res[[3]][r]))
          
          result <- c(mu.est, alpha.est)
        }
    }
  }
  return(result)
}

####### OSM ### 
a <- est.func(N=100, npar=6, n.r=15, NCols=10, formula="Y~row", model="OSM", nclus.row=2)
b <- est.func(N=100, npar=6, n.r=50, NCols=10, formula="Y~row", model="OSM", nclus.row=2)
c <- est.func(N=100, npar=6, n.r=100, NCols=10, formula="Y~row", model="OSM", nclus.row=2)

d <- est.func(N=100, npar=6, n.r=15, NCols=20, formula="Y~row", model="OSM", nclus.row=2)
e <- est.func(N=100, npar=6, n.r=50, NCols=20, formula="Y~row", model="OSM", nclus.row=2)
f <- est.func(N=100, npar=6, n.r=100, NCols=20, formula="Y~row", model="OSM", nclus.row=2)

g <- est.func(N=100, npar=6, n.r=15, NCols=30, formula="Y~row", model="OSM", nclus.row=2)
h <- est.func(N=100, npar=6, n.r=50, NCols=30, formula="Y~row", model="OSM", nclus.row=2)
i <- est.func(N=100, npar=6, n.r=100, NCols=30, formula="Y~row", model="OSM", nclus.row=2)


######## POM ######
j <- est.func(N=100, npar=4, n.r=15, NCols=10, formula="Y~row", model="POM", nclus.row=2)
k <- est.func(N=100, npar=4, n.r=50, NCols=10, formula="Y~row", model="POM", nclus.row=2)
l <- est.func(N=100, npar=4, n.r=100, NCols=10, formula="Y~row", model="POM", nclus.row=2)

m <- est.func(N=100, npar=4, n.r=15, NCols=20, formula="Y~row", model="POM", nclus.row=2)
n <- est.func(N=100, npar=4, n.r=50, NCols=20, formula="Y~row", model="POM", nclus.row=2)
o <- est.func(N=100, npar=4, n.r=100, NCols=20, formula="Y~row", model="POM", nclus.row=2)

p <- est.func(N=100, npar=4, n.r=15, NCols=30, formula="Y~row", model="POM", nclus.row=2)
q <- est.func(N=100, npar=4, n.r=50, NCols=30, formula="Y~row", model="POM", nclus.row=2)
r <- est.func(N=100, npar=4, n.r=100, NCols=30, formula="Y~row", model="POM", nclus.row=2)



###### Binary###
s <- est.func(N=100, npar=2, n.r=15, NCols=10, formula="Y~row", model="Binary", nclus.row=2)
t <- est.func(N=100, npar=2, n.r=50, NCols=10, formula="Y~row", model="Binary", nclus.row=2)
u <- est.func(N=100, npar=2, n.r=100, NCols=10, formula="Y~row", model="Binary", nclus.row=2)

v <- est.func(N=100, npar=2, n.r=15, NCols=20, formula="Y~row", model="Binary", nclus.row=2)
w <- est.func(N=100, npar=2, n.r=50, NCols=20, formula="Y~row", model="Binary", nclus.row=2)
x <- est.func(N=100, npar=2, n.r=100, NCols=20, formula="Y~row", model="Binary", nclus.row=2)

y <- est.func(N=100, npar=2, n.r=15, NCols=30, formula="Y~row", model="Binary", nclus.row=2)
z <- est.func(N=100, npar=2, n.r=50, NCols=30, formula="Y~row", model="Binary", nclus.row=2)
zz <- est.func(N=100, npar=2, n.r=100, NCols=30, formula="Y~row", model="Binary", nclus.row=2)




a2 <- SE.func(N=100, npar=6, n.r=15, NCols=10, formula="Y~row", model="OSM", nclus.row=2)
b2 <- SE.func(N=100, npar=6, n.r=50, NCols=10, formula="Y~row", model="OSM", nclus.row=2)
c2 <- SE.func(N=100, npar=6, n.r=100, NCols=10, formula="Y~row", model="OSM", nclus.row=2)


d2 <- SE.func(N=100, npar=6, n.r=15, NCols=20, formula="Y~row", model="OSM", nclus.row=2)
e2 <- SE.func(N=100, npar=6, n.r=50, NCols=20, formula="Y~row", model="OSM", nclus.row=2)
f2 <- SE.func(N=100, npar=6, n.r=100, NCols=20, formula="Y~row", model="OSM", nclus.row=2)

g2 <- SE.func(N=100, npar=6, n.r=15, NCols=30, formula="Y~row", model="OSM", nclus.row=2)
h2 <- SE.func(N=100, npar=6, n.r=50, NCols=30, formula="Y~row", model="OSM", nclus.row=2)
i2 <- SE.func(N=100, npar=6, n.r=100, NCols=30, formula="Y~row", model="OSM", nclus.row=2)


######## POM ######
j2 <- SE.func(N=100, npar=4, n.r=15, NCols=10, formula="Y~row", model="POM", nclus.row=2)
k2 <- SE.func(N=100, npar=4, n.r=50, NCols=10, formula="Y~row", model="POM", nclus.row=2)
l2 <- SE.func(N=1000, npar=4, n.r=100, NCols=10, formula="Y~row", model="POM", nclus.row=2)

j2 <- SE.func(N=100, npar=4, n.r=15, NCols=20, formula="Y~row", model="POM", nclus.row=2)
k2 <- SE.func(N=100, npar=4, n.r=50, NCols=20, formula="Y~row", model="POM", nclus.row=2)
l2 <- SE.func(N=100, npar=4, n.r=100, NCols=20, formula="Y~row", model="POM", nclus.row=2)

m2 <- SE.func(N=100, npar=4, n.r=15, NCols=30, formula="Y~row", model="POM", nclus.row=2)
n2 <- SE.func(N=100, npar=4, n.r=50, NCols=30, formula="Y~row", model="POM", nclus.row=2)
o2 <- SE.func(N=100, npar=4, n.r=100, NCols=30, formula="Y~row", model="POM", nclus.row=2)


###### Binary###
p2 <- SE.func(N=100, npar=2, n.r=15, NCols=10, formula="Y~row", model="Binary", nclus.row=2)
q2 <- SE.func(N=100, npar=2, n.r=50, NCols=10, formula="Y~row", model="Binary", nclus.row=2)
r2 <- SE.func(N=100, npar=2, n.r=100, NCols=10, formula="Y~row", model="Binary", nclus.row=2)

s2 <- SE.func(N=100, npar=2, n.r=15, NCols=20, formula="Y~row", model="Binary", nclus.row=2)
t2 <- SE.func(N=100, npar=2, n.r=50, NCols=20, formula="Y~row", model="Binary", nclus.row=2)
u2 <- SE.func(N=100, npar=2, n.r=100, NCols=20, formula="Y~row", model="Binary", nclus.row=2)

v2 <- SE.func(N=100, npar=2, n.r=15, NCols=30, formula="Y~row", model="Binary", nclus.row=2)
w2 <- SE.func(N=100, npar=2, n.r=50, NCols=30, formula="Y~row", model="Binary", nclus.row=2)
y2 <- SE.func(N=100, npar=2, n.r=100, NCols=30, formula="Y~row", model="Binary", nclus.row=2)


