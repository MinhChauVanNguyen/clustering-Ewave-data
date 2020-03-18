
#Set up a simulation study to test how reliably, 
#in a diverse range of scenarios (different starting values combinations), 
#we were able to estimate the parameters of stereotype models 
#using the EM algorithm. 

library(clustord)
library(reshape2)

## Constraints:
# mu, u and alpha/beta can take values from -Inf to Inf
# mu_1 = 0
# 0=phi_1<=phi_2<=...<=phi_q=1
# phi_1 = 0, phi_4 = 1
# sum(alpha_r) = 0, sum(beta_r) = 0
# sum(alpha_i) = 0, sum(beta_i) = 0

## fixed
# mu2,mu3,mu4,u2,u3


# perform row clustering on 1 generated data set 10 times
# with chosen values of 10 sets of starting parameter values

# Y~ROW

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
###################################### gen.row ####################################
gen.row <- function(npar,seed){
  set.seed(seed)
  ## produce a set of random starting parameter values
  initvect <- rnorm(n=npar)
  #names(initvect) <- paste0(c("mu2","mu3","mu4","u2","u3","alpha1","alpha2"))
  initvect
}

# npar is the number of free initial parameters, excluding pi
gen.row(npar=10, seed=1234)


################################### data.row ######################################

# fixed parameter values of mu1,mu2,mu3,mu4,u1,u2,u3,u4,alpha1 & alpha2
fixed.par <- c(0.2,0.6,0.8,1.2,-5,-2,0.2,0.3,1.2,2.5)  
# apply constraints on fixed.par
fixed.par2 <- c(mu1=0,mu2=0.6,mu3=0.8,mu4=1.2,u1=0,u2=0.1,u3=0.2,u4=1,alpha1=0.5,alpha2=-0.5)

# data.row generate data of TWO clusters
data.row <- function(NRows, NCols, fixed.par){
  ## fixed values of starting parameters
  empty.vect <- c()
  fixed.par <- as.vector(fixed.par)
  vect <- c(empty.vect,fixed.par)
  #names(vect) <- c("mu1","mu2","mu3","mu4","u1","u2","u3","u4","alpha1","alpha2")
  
  ## calculate p 
  ## p is the response prob of y having category k=1,2,3,4 belongs to cluster r=1,2
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
  
  theta1 <- c(p1.clust1,p2.clust1,p3.clust1,p4.clust1)
  names(theta1) <- c("p1.clust1","p2.clust1","p3.clust1","p4.clust1")
  theta2 <- c(p1.clust2,p2.clust2,p3.clust2,p4.clust2)
  names(theta2) <- c("p1.clust2","p2.clust2","p3.clust2","p4.clust2")

  ## check if sum(theta1)=1 & sum(theta2)=1
  sum(theta1); sum(theta2)
  
  ## form clusters based on their response probabilities
  ## exp(mu_k-phi_k*a_r)/sum(ep(mu_k-phi_k*a_r))
  cluster1 <- matrix(sample(1:4,size=NCols*NRows,replace=TRUE,prob=theta1),nrow=NRows)
  cluster2 <- matrix(sample(1:4,size=NCols*NRows,replace=TRUE,prob=theta2),nrow=NRows)
  
  ## put cluster 1 & 2 together to produce a data frame with three factor columns
  data <- rbind(cluster1,cluster2)
  data4 <- melt(data)
  #data4
  colnames(data4) <- c("ROW","COL","Y")
  data4$ROW <- 1:nrow(data4)
  data4[,c(1,2,3)] <- lapply(data4[,c(1,2,3)], factor)
  return(data.frame(data4))
  
  #data2 <- as.data.frame(data)
  #ROW <- rownames(data2)
  #data2 <- cbind(ROW=ROW, data2)
  #data2
  #data3 <- gather(data = data2, key = COL, value = Y, V1:V5, factor_key = TRUE)
  #data3
  #data3$ROW <- 1:nrow(data3)
  #data3[,c(1,2,3)] <- lapply(data3[,c(1,2,3)], factor)
  #return(data3)
}

fixed.par <- c(0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,0.5,-0.5)
fixed.par2 <- c(mu1=0,mu2=0.6,mu3=0.8,mu4=1.2,u1=0,u2=-2,u3=0.2,u4=1,alpha1=0.5,alpha2=-0.5)


data.row <- function(NRows, NCols, fixed.par, model){
  ## fixed values of starting parameters
  empty.vect <- c()
  fixed.par <- as.vector(fixed.par)
  vect <- c(empty.vect,fixed.par)

  if (!is.character(model) || !is.vector(model) || length(model)!=1) 
    stop("model must be a string, 'OSM' or 'POM'.")
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
    sum.p.cluster12 <- exp(vect[1]-vect[9])+exp(vect[2]-vect[9])+exp(vect[3]-vect[9])+exp(vect[4]-vect[9])
    p1.clust12 <- exp(vect[1]-vect[9])/sum.p.cluster12
    p2.clust12 <- exp(vect[2]-vect[9])/sum.p.cluster12
    p3.clust12 <- exp(vect[3]-vect[9])/sum.p.cluster12
    p4.clust12 <- exp(vect[4]-vect[9])/sum.p.cluster12
    
    sum.p.cluster22 <- exp(vect[1]-vect[10])+exp(vect[2]-vect[10])+exp(vect[3]-vect[10])+exp(vect[4]-vect[10])
    p1.clust22 <- exp(vect[1]-vect[10])/sum.p.cluster22
    p2.clust22 <- exp(vect[2]-vect[10])/sum.p.cluster22
    p3.clust22 <- exp(vect[3]-vect[10])/sum.p.cluster22
    p4.clust22 <- exp(vect[4]-vect[10])/sum.p.cluster22
    
    theta12 <- c(p1.clust12,p2.clust12,p3.clust12,p4.clust12)
    theta22 <- c(p1.clust22,p2.clust22,p3.clust22,p4.clust22)
  
  # Binary model
    theta13 <- exp(vect[1]+vect[9])/(1+exp(vect[1]+vect[9]))
    #theta23 <- exp(vect[1]+vect[10])/(1+exp(vect[1]+vect[10]))
    theta23 <- 1 - exp(vect[1]+vect[9])/(1+exp(vect[1]+vect[9]))
    
  ## form clusters based on their response probabilities
  ## exp(mu_k-phi_k*a_r)/sum(ep(mu_k-phi_k*a_r))
  
  if(model=="OSM"){
    cluster1 <- matrix(sample(1:4,size=NCols*NRows,replace=TRUE,prob=theta1),nrow=NRows)
    cluster2 <- matrix(sample(1:4,size=NCols*NRows,replace=TRUE,prob=theta2),nrow=NRows)
  }else if(model=="POM"){
    cluster1 <- matrix(sample(1:4,size=NCols*NRows,replace=TRUE,prob=theta12),nrow=NRows)
    cluster2 <- matrix(sample(1:4,size=NCols*NRows,replace=TRUE,prob=theta22),nrow=NRows)
  }else{
    cluster1 <- matrix(sample(c(0,1),size=NCols*NRows,replace=TRUE,prob=c(theta13,theta23)),nrow=NRows)
    cluster2 <- matrix(sample(c(0,1),size=NCols*NRows,replace=TRUE,prob=c(theta13,theta23)),nrow=NRows)
  }
    
  ## put cluster 1 & 2 together to produce a data frame with three factor columns
  data <- rbind(cluster1,cluster2)
  data4 <- melt(data)
  data4
  colnames(data4) <- c("ROW","COL","Y")
  data4$Y <- factor(data4$Y)
  return(data.frame(data4))
}

fixed.par <- c(0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,0.5,-0.5)

data.row(NRows=20, NCols=5, fixed.par=fixed.par, model="POM")
data.row(NRows=20, NCols=5, fixed.par=fixed.par, model="OSM")
data.row(NRows=20, NCols=5, fixed.par=fixed.par, model="Binary")


  
######################################## sim.row #####################################

# perform simulation on a fixed data set 10 times 
# using 10 different initial vector parameters (10 different seeds)
# return a list of 10 lists 
# each list should contain 3 elements:
# RowCluster, best.lli and criteria

### FIRST 
N <- 10                                         # number of iterations
npar <- 6
initvect <- matrix(NA, nrow=N, ncol=npar)       # store initvect in a matrix
fixed.par <- c(0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,0.5,-0.5)
data <- data.row(NRows=20, NCols=5, fixed.par=fixed.par, model="OSM")
cluster <- list()
cluster.results <- rep(list(list()),N)

for(i in 1:nrow(initvect)){
  initvect[i,] <- gen.row(npar=npar, seed=i)
  for(j in i){
  cluster[[j]] <- rowclustering(formula = "Y~row", 
                           model = "OSM", 
                           nclus.row = 2, 
                           long.df = data, 
                           initvect = initvect[i,])
  cluster.results[[j]][[1]] <- cluster[[j]]$RowClusters
  cluster.results[[j]][[2]] <- cluster[[j]]$EM.status$best.lli
  cluster.results[[j]][[3]] <- cluster[[j]]$criteria
  }
  names(cluster.results[[j]]) <- c("Row Clusters", "best.lli", "criteria")
}

length(cluster)
length(cluster.results)
cluster[[2]]$ppr
cluster[[2]]$info
cluster[[2]][[1]]
cluster[[2]][[2]]
cluster[[2]][[3]]
cluster[[2]][[10]]
cluster[[2]][[14]]
length(cluster[[2]])
is.list(cluster[[2]])



### SECOND
# CAN'T NAME THE LIST ELEMENTS OF THE LIST
# n DOES NOT EQUAL TO NRows
# best.lli IS THE SAME AMONG THE LISTS??????

N <- 10
npar <- 6                                     #(does not include pi)
cluster <- list()
initvect <- matrix(NA, nrow=N, ncol=npar)     # store initvect in a matrix

# extract three elements (RowCluster, best.lli, criteria) of every list 
# from cluster (a list) and store them in a list of 10 lists called 
# cluster.results

cluster.results <-  rep(list(list()), N)

for(i in 1:N){
  initvect[i,] <- gen.row(npar=npar, seed=i)
  for(j in 1:N){
    cluster[[j]] <- rowclustering(formula = "Y~row", 
                                  model = "OSM", 
                                  nclus.row = 2, 
                                  long.df = data, 
                                  initvect = initvect[i,])
    
  cluster.results[[j]][[1]] <- cluster[[j]]$RowClusters
  cluster.results[[j]][[2]] <- cluster[[j]]$EM.status$best.lli
  cluster.results[[j]][[3]] <- cluster[[j]]$criteria
  }
  names(cluster.results[[j]]) <- c("Row Clusters", "best.lli", "criteria")
  # Obtain the average of EM.status$best.lli of N iterations
  total.lli <- sapply(cluster.results, function(x) x[["best.lli"]])
  best.lli.average <- best.lli.total/N
  print(total.lli)
  #print(best.lli.average)
}

length(cluster.results)
cluster.results[[1]][[2]]
cluster.results[[2]][[2]]
cluster.results[[3]][[2]]


### THIRD

# calculate the sum of EM.status$best.lli for 10 iterations
# lapply returns a list
# sapply returns a vector
total.lli <- numeric()
for(i in 1:length(cluster.results)){
  # total.lli <- sum(lapply(cluster.results, function(x) x[[i]][["best.lli"]]))
  # total.lli <- sum(sapply(cluster.results[[i]], function(x) x[["best.lli"]]))
  total.lli <- sum(sapply(cluster.results, function(x) x[[i]][["best.lli"]]))
}

total.lli
for(i in 1:length(cluster.results)){
  # total.lli <- sum(lapply(cluster.results, function(x) x[[i]][["best.lli"]]))
  total.lli[[i]] <- sapply(cluster.results, function(x) x[["best.lli"]])
  total.lli <- sum(total.lli[[i]])
}




### FOURTH
## Return the maximum value of best.lli element from the list of N lists
cluster.results <-  list()
MAX <- -Inf        ## initialize maximum

for(i in 1:N){
  initvect[i,] <- gen.row(npar=npar, seed=i)
  cluster <- rowclustering(formula = "Y~row", 
                                  model = "OSM", 
                                  nclus.row = 2, 
                                  long.df = data, 
                                  initvect = initvect[i,])
    cluster.results[[i]] <- list()
    cluster.results[[i]][[1]] <- cluster$RowClusters
    cluster.results[[i]][[2]] <- cluster$EM.status$best.lli
    cluster.results[[i]][[3]] <- cluster$criteria
    
    names(cluster.results[[i]]) <- c("Row Clusters", "best.lli", "criteria")
    
    lli.vals <- cluster.results[[i]]$best.lli
    if (lli.vals > MAX) 
      MAX <- lli.vals
    # return the list with the maximum best.lli 
      best.list <- cluster.results[[which.max(cluster.results[[i]]$best.lli)]]
    
    print(best.list)
    #print(MAX)
}

#names(cluster.results) <- LETTERS[1:N]
#names(cluster.results) < - assign(paste0("List",j),as.integer(j))
#names(cluster.results) <- paste0("List", j)
#max.lli <- lapply(cluster.results[[j]], function(list) max(list[[2]]))
#print(max.lli)
sim.row <- function(N, npar, NRows, NCols, formula="Y~row", model, nclus.row){
  # store initial parameters in a matrix with ncol equivalent to
  # the number of free parameters (excluding pi) & nrow equivalent to 
  # the number of iterations N
  # hence produce a matrix of N initial parameter vectors
  initvect <- matrix(NA, nrow=N, ncol=npar)
  # create an empty list to store the N lists of row clustering results
  cluster <- list()
  # create another empty list of N lists of row clustering results 
  # showing only the three required elements for each list
  cluster.results <- list()
  # so length(cluster) = length(cluster.results) = N
  length(cluster) ; length(cluster.results)
  
  # generate pseudo-data 
  data <- data.row(NRows=NRows, NCols=NCols, fixed.par=fixed.par, model=model)
  # initialize maximum
  MAX <- -Inf                                   
  
  for(i in 1:N){
    initvect[i,] <- gen.row(npar=npar, seed=i)

    cluster <- rowclustering(formula = formula, 
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
    lli.vals <- cluster.results[[i]]$best.lli
    if (lli.vals > MAX) 
      MAX <- lli.vals
    }
    # return the average value of EM.status$best.lli of N iterations
    # where EM.statust$best.lli is the second element on each list
    best.lli.total <-  sum(sapply(cluster.results, function(x) x[["best.lli"]]))
    best.lli.average <- best.lli.total/N
    final.result <- cbind(MAX,best.lli.average)
    colnames(final.result) <- c("Max.lli","Average.lli")
    return(final.result)
}

sim.row(N=5, npar=6,NRows=20, NCols=10, formula="Y~row", model="OSM", nclus.row=2)
  
# sim.row can only work for formula='Y~row', 
# needs more flexibility
set.seed(123)
NCols <- 5
q <- 4

mu <- c(0.2,0.3,0.4,0.5)
u <- c(0.1,0.2,0.3,0.9)
beta <- rnorm(n=NCols)
alpha1 <- 0.5
alpha2 <- 1 - alpha1
exp(mu[1] - alpha1 - beta[1])
exp(mu[2] - alpha1 - beta[2])

m <- sum(exp(mu[1] - alpha1 - beta[1]):exp(mu[q] - alpha1 - beta[NCols]))
length(exp(mu[1] - alpha1 - beta[1]):exp(mu[q] - alpha1 - beta[NCols]))
for(i in 1:NCols){
  for(j in 1:q){
    # cluster 1
    # should return NCols x q combinations
    # p.cluster1 should be a vector with NCols x q elements
    p.cluster1 <- (exp(mu[j] - alpha2 - beta[i]))/m
    # should return NCols x q combinations
    # cluster 2
    total2 <- sum(exp(mu[j] - alpha2 - beta[i]))
    p.cluster2 <- (exp(mu[j] - alpha2 - beta[i]))/sum(total2)
  }
  print(p.cluster1)
}
sum(m)

exp(mu[1]-alpha1-beta[1])


  #name <- paste("lis", j , sep = "_")
  #assign(name, cluster.results)
  #names(cluster.results[[j]][[1]]) <- "RowClusters"
  #names(cluster.results[[j]][[1]][[1]]) <- "Cluster.1"
  #names(cluster.results[[j]][[1]][[2]]) <- "Cluster.2"
  #names(cluster.results[[j]][[2]]) <- "Best.Likelihood.Value"
  #names(cluster.results[[j]][[3]]) <- "Criteria"
  #names(cluster.results[[j]][[3]][[1]]) <- "Res.Dev"
  #names(cluster.results[[j]][[3]][[2]]) <- "AIC"
  #names(cluster.results[[j]][[3]][[3]]) <- "AIC.c"
  #names(cluster.results[[j]][[3]][[4]]) <- "BIC"
  #names(cluster.results[[j]][[3]][[5]]) <- "ICL"
#  }
  #return(cluster.results)
#}


N <- 5
npar <- 6
data <- data.row(NRows=20,NCols=10,fixed.par=fixed.par,model="OSM")
cluster <- list()
for(i in 1:N){
  initvect[i,] <- gen.row(npar=npar, seed=i)
  cluster[[i]] <- rowclustering(formula = "Y~row", 
                                model = "OSM", 
                                nclus.row = 2, 
                                long.df = data, 
                                initvect = initvect[i,])
}


MAX <- -Inf        ## initialize maximum

for(i in 1:N){
  initvect[i,] <- gen.row(npar=npar, seed=i)
  cluster <- rowclustering(formula = "Y~row", 
                           model = "OSM", 
                           nclus.row = 2, 
                           long.df = data, 
                           initvect = initvect[i,])
  cluster.results[[i]] <- list()
  cluster.results[[i]][[1]] <- cluster$RowClusters
  cluster.results[[i]][[2]] <- cluster$EM.status$best.lli
  cluster.results[[i]][[3]] <- cluster$criteria
  
  names(cluster.results[[i]]) <- c("Row Clusters", "best.lli", "criteria")
  
  lli.vals <- cluster.results[[i]]$best.lli
  if (lli.vals > MAX)
    MAX <- lli.vals
  # return the list with the maximum best.lli 
  best.list <- cluster.results[[which.max(cluster.results[[i]]$best.lli)]]
  
  print(best.list)
  #print(MAX)
}

