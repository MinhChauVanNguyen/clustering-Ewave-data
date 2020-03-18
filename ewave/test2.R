library(clustord)
library(reshape2)

# generate starting point combinations
gen.row <- function(npar, seed){
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

gen.row(npar=6,seed=123)

fixed.par <- c(mu1=0, mu2=0.6, mu3=0.8, mu4=1.2,
               phi1=0, phi2=0.4, phi3=0.6, phi4=1,
               alpha1=0.5, alpha2=-0.5)

fixedPO.par <- c(mu1=0.4, mu2=0.6, mu3=0.8, mu4=1.2,
               phi1=0, phi2=0.4, phi3=0.6, phi4=1,
               alpha1=0.5, alpha2=-0.5)

# generate pseudo-data
data.row <- function(NRows, NCols, fixed.par, model){
  ## fixed values of starting parameters
  empty.vect <- c()
  fixed.par <- as.vector(fixed.par)
  vect <- c(empty.vect, fixed.par)
  
  if (!is.character(model) || !is.vector(model) || length(model)!=1) 
    stop("model must be a string, 'OSM' or 'POM' or 'Binary'.")
  if (!(model %in% c("OSM","POM","Binary"))) 
    stop("model must be either 'OSM' or POM' or 'Binary'.")
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
  
  # POM model
  p1.clust12 <- exp(vect[1]-vect[9])/(1+exp(vect[1]-vect[9]))
  p2.clust12 <- exp(vect[2]-vect[9])/p1.clus12
  p3.clust12 <- exp(vect[3]-vect[9])/p2.clus12
  p4.clust12 <- 1 - (p1.clus12+p2.clus12+p3.clust12)
  
  p1.clust22 <- exp(vect[1]-vect[10])/(1+exp(vect[1]-vect[10]))
  p2.clust22 <- exp(vect[2]-vect[10])/p1.clus22
  p3.clust22 <- exp(vect[3]-vect[10])/p2.clus22
  p4.clust22 <- 1 - (p1.clus22+p2.clus22+p3.clust22)
  
  theta12 <- c(p1.clust12,p2.clust12,p3.clust12,p4.clust12)
  theta22 <- c(p1.clust22,p2.clust22,p3.clust22,p4.clust22)
  
  # Binary model
  binary.p1 <- exp(vect[2]+vect[9])/(1+exp(vect[2]+vect[9]))
  binary.p0 <- 1 - binary.theta1
  theta <- c(binary.theta0,binary.theta1)
  
  ## form clusters based on their response probabilities
  if(model=="OSM"){
    cluster1 <- matrix(sample(1:4,size=NCols*NRows,replace=TRUE,prob=theta1),nrow=NRows)
    cluster2 <- matrix(sample(1:4,size=NCols*NRows,replace=TRUE,prob=theta2),nrow=NRows)
  }else if(model=="POM"){
    cluster1 <- matrix(sample(1:4,size=NCols*NRows,replace=TRUE,prob=theta12),nrow=NRows)
    cluster2 <- matrix(sample(1:4,size=NCols*NRows,replace=TRUE,prob=theta22),nrow=NRows)
  }else{
    cluster1 <- matrix(sample(c(0,1),
                              size=NCols*NRows,
                              replace=TRUE,
                              prob=theta),
                              nrow=NRows)
    cluster2 <- matrix(sample(c(0,1),
                              size=NCols*NRows,
                              replace=TRUE,
                              prob=theta),
                              nrow=NRows)
  }
  
  ## put cluster 1 & 2 together to produce a data frame with three factor columns
  data <- rbind(cluster1,cluster2)
  pseudo.data <- melt(data)
  colnames(pseudo.data) <- c("ROW","COL","Y")
  pseudo.data$Y <- factor(pseudo.data$Y)
  return(data.frame(pseudo.data))
}

fake.data <- data.row(NRows=15, NCols=10, fixed.par=fixed.par, model="OSM")

fixed.par <- c(mu1=0,mu2=0.6,mu3=0.8,mu4=1.2,
               phi1=0,phi2=0.4,phi3=0.6,phi4=1,
               alpha1=0.5,alpha2=-0.5)

N <- 5
npar <- 2
formula <- "Y~row"
model <- "Binary"
nclus.row <- 2
NRows <- 20
NCols <- 10

sim.row2 <- function(N, npar, NRows, NCols=10, formula="Y~row", model, nclus.row)
  {
  initvect <- matrix(NA, nrow=N, ncol=npar)
  # cluster.results <- list()
  cluster <- list()
  # generate pseudo-data 
  data <- data.row(NRows=NRows, NCols=NCols, fixed.par=fixed.par, model=model)
  # initialize maximum
  # MAX <- -Inf                                   
  
  for(i in 1:N){
    initvect[i,] <- gen.row(npar=npar, seed=i)
    cluster[[i]] <- rowclustering(formula = formula, 
                             model = model, 
                             nclus.row = nclus.row, 
                             long.df = data, 
                             initvect = initvect[i,])
  
    cluster.results[[i]] <- list()
    cluster.results[[i]][[1]] <- cluster$RowClusters
    cluster.results[[i]][[2]] <- cluster$EM.status$best.lli
    cluster.results[[i]][[3]] <- cluster$criteria
    
    # name the required elements for each list
    names(cluster.results[[i]]) <- c("Row Clusters", "best.lli", "criteria")
    
    # return the maximum value of EM.status$best.lli of N iterations
    # lli.vals <- cluster.results[[i]]$best.lli
    # if (lli.vals > MAX){
    #  MAX <- lli.vals
    # }
    # return the list with the maximum best.lli 
    best.list <- cluster.results[[which.max(cluster.results[[i]]$best.lli)]]
  }
  # return the average value of EM.status$best.lli of N iterations
  # where EM.statust$best.lli is the second element on each list
  # best.lli.total <-  sum(sapply(cluster.results, function(x) x[["best.lli"]]))
  # best.lli.average <- best.lli.total/N
  
  final.result <- cbind(MAX,best.lli.average)
  colnames(final.result) <- c("Max.lli","Average.lli")
  return(best.list)
  #return(final.result)
}





