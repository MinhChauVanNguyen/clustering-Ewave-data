library(clustord)
library(reshape2)

# generate initvect
gen.row <- function(npar,seed){
  if(!is.numeric(npar)){
    stop("npar must be a positive integer")
  }else if(length(npar)>1){ 
    warning("npar has more than one element, use only the first element")
    npar <- npar[1]
    if(npar<=0 | round(npar)!=npar){
      stop("npar must be a positive integer")
    }
  }else if(npar<=0 | round(npar)!=npar){
    stop("npar must be a positive integer")
  }
  if(!is.numeric(seed)){
    stop("seed must be a real number")
  }else if(length(seed)>1){
    warning("Seed has more than one element, use only the first element.")
    seed <- seed[1]    
  }
  set.seed(seed)
  initvect <- rnorm(n=npar)
  initvect
}

fixed.par <- c(mu1=0, mu2=0.6, mu3=0.8, mu4=1.2,
               phi1=0, phi2=0.4, phi3=0.6, phi4=1,
               alpha1=0.5, alpha2=-0.5)

# generate data 

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

    cluster1 <- matrix(sample(1:4,size=NCols*NRows,replace=TRUE,prob=theta1),nrow=NRows)
    cluster2 <- matrix(sample(1:4,size=NCols*NRows,replace=TRUE,prob=theta2),nrow=NRows)
    
  }else if(model=="POM"){
    p1 <- expit(vect[1]-vect[9])
    p2 <- expit(vect[2]-vect[9]) - p1
    p3 <- -expit(p2) + expit(p1)
    p4 <- 1 - (p1 + p2 + p3)
    theta.PO <- c(p1, p2, p3, p4)
    
    cluster1 <- matrix(sample(1:4,size=NCols*NRows,replace=TRUE,prob=theta.PO),nrow=NRows)
    cluster2 <- matrix(sample(1:4,size=NCols*NRows,replace=TRUE,prob=theta.PO),nrow=NRows)
    
  }else{
    # Binary model
    theta1.Bin <- expit(vect[1]+vect[9])
    theta0.Bin <- 1 - theta1.Bin
    
    cluster1 <- matrix(sample(c(0,1),size=NCols*NRows,replace=TRUE,prob=c(theta0.Bin,theta1.Bin)),nrow=NRows)
    cluster2 <- matrix(sample(c(0,1),size=NCols*NRows,replace=TRUE,prob=c(theta0.Bin,theta1.Bin)),nrow=NRows)
  }
  
  ## put cluster 1 & 2 together to produce a data frame with three factor columns
  data <- rbind(cluster1,cluster2)
  pseudo.data <- melt(data)
  colnames(pseudo.data) <- c("ROW","COL","Y")
  pseudo.data$Y <- factor(pseudo.data$Y)
  return(data.frame(pseudo.data))
}


###################################################################################
calc.SE.rowcluster <- function(long.df, clust.out,
                               optim.control=default.optim.control()) {
  optim.control$fnscale=-1
  
  y.mat <- df2mat(long.df)
  outvect <- clust.out$outvect
  
  optim.hess <- optimHess(par=outvect,
                          fn=calc.ll,
                          long.df=long.df,
                          y.mat=y.mat,
                          model=clust.out$model,
                          submodel=clust.out$submodel,
                          ppr.m=clust.out$ppr,
                          pi.v=clust.out$pi.out,
                          RG=clust.out$info["R"],
                          constraint.sum.zero=clust.out$constraint.sum.zero,
                          SE.calc=TRUE,
                          control=optim.control)
  
  SE <- sqrt(diag(solve(-optim.hess)))
  SE
}

###################################################################################
N <- 5
npar <- 6
formula <- "Y~row"
model <- "OSM"
nclus.row <- 2
NRows <- 20
NCols <- 10

fixed.par <- c(mu1=0,mu2=0.6,mu3=0.8,mu4=1.2,
               phi1=0,phi2=0.4,phi3=0.6,phi4=1,
               alpha1=0.5,alpha2=-0.5)

# generate pseudo-data 

SE.func <- function(N, npar, n.r, NCols, formula="Y~row", model, nclus.row){
  data <- data.row(n.r=n.r, NCols=NCols, fixed.par=fixed.par, model=model)
  initvect <- matrix(NA, nrow=N, ncol=npar)
  cluster <- list()
  SE.est <- matrix(NA, nrow=N, ncol=npar)
# inner loop is executed N- times for every execution of outer loop.  
  for(i in 1:N){
    for(j in 1:N){
      initvect[i,] <- gen.row(npar=npar, seed=i)
      cluster[[i]] <- rowclustering(formula = formula, 
                                    model = model, 
                                    nclus.row = nclus.row, 
                                    long.df = data, 
                                    initvect = initvect[i,])
      SE.est[j,] <- calc.SE.rowcluster(long.df = data, clust.out = cluster[[i]])
    }
  }
    SE.mean <- colMeans(SE.est)
    return(SE.mean)
}

####### OSM ### 
a2 <- SE.func(N=5, npar=6, n.r=15, NCols=10, formula="Y~row", model="OSM", nclus.row=2)
b2 <- SE.func(N=1000, npar=6, n.r=50, NCols=10, formula="Y~row", model="OSM", nclus.row=2)
c2 <- SE.func(N=1000, npar=6, n.r=100, NCols=10, formula="Y~row", model="OSM", nclus.row=2)


d2 <- SE.func(N=1000, npar=6, n.r=15, NCols=20, formula="Y~row", model="OSM", nclus.row=2)
e2 <- SE.func(N=1000, npar=6, n.r=50, NCols=20, formula="Y~row", model="OSM", nclus.row=2)
f2 <- SE.func(N=1000, npar=6, n.r=100, NCols=20, formula="Y~row", model="OSM", nclus.row=2)

g2 <- SE.func(N=1000, npar=6, n.r=15, NCols=30, formula="Y~row", model="OSM", nclus.row=2)
h2 <- SE.func(N=1000, npar=6, n.r=50, NCols=30, formula="Y~row", model="OSM", nclus.row=2)
i2 <- SE.func(N=1000, npar=6, n.r=100, NCols=30, formula="Y~row", model="OSM", nclus.row=2)


######## POM ######
j2 <- SE.func(N=1000, npar=4, n.r=15, NCols=10, formula="Y~row", model="POM", nclus.row=2)
k2 <- SE.func(N=1000, npar=4, n.r=50, NCols=10, formula="Y~row", model="POM", nclus.row=2)
l2 <- SE.func(N=1000, npar=4, n.r=100, NCols=10, formula="Y~row", model="POM", nclus.row=2)

j2 <- SE.func(N=1000, npar=4, n.r=15, NCols=20, formula="Y~row", model="POM", nclus.row=2)
k2 <- SE.func(N=1000, npar=4, n.r=50, NCols=20, formula="Y~row", model="POM", nclus.row=2)
l2 <- SE.func(N=1000, npar=4, n.r=100, NCols=20, formula="Y~row", model="POM", nclus.row=2)

m2 <- SE.func(N=1000, npar=4, n.r=15, NCols=30, formula="Y~row", model="POM", nclus.row=2)
n2 <- SE.func(N=1000, npar=4, n.r=50, NCols=30, formula="Y~row", model="POM", nclus.row=2)
o2 <- SE.func(N=1000, npar=4, n.r=100, NCols=30, formula="Y~row", model="POM", nclus.row=2)


###### Binary###
p2 <- SE.func(N=5, npar=2, n.r=15, NCols=10, formula="Y~row", model="Binary", nclus.row=2)
q2 <- SE.func(N=1000, npar=2, n.r=50, NCols=10, formula="Y~row", model="Binary", nclus.row=2)
r2 <- SE.func(N=1000, npar=2, n.r=100, NCols=10, formula="Y~row", model="Binary", nclus.row=2)

s2 <- SE.func(N=1000, npar=2, n.r=15, NCols=20, formula="Y~row", model="Binary", nclus.row=2)
t2 <- SE.func(N=1000, npar=2, n.r=50, NCols=20, formula="Y~row", model="Binary", nclus.row=2)
u2 <- SE.func(N=1000, npar=2, n.r=100, NCols=20, formula="Y~row", model="Binary", nclus.row=2)

v2 <- SE.func(N=1000, npar=2, n.r=15, NCols=30, formula="Y~row", model="Binary", nclus.row=2)
w2 <- SE.func(N=1000, npar=2, n.r=50, NCols=30, formula="Y~row", model="Binary", nclus.row=2)
y2 <- SE.func(N=1000, npar=2, n.r=100, NCols=30, formula="Y~row", model="Binary", nclus.row=2)



