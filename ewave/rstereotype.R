#####################################################################################
# Ordered stereotype functions
#####################################################################################
#' Check parameters in the ordered stereotype model
#' @param muvec A vector of intensities, with muvec[1]=0 but no other constraints
#' @param phivec A vector of scores, with phivec[1]=0, phivec[q]=1 
#' @return "OK" if all is well but the text of an error message otherwise
#' @export
check_stereotype_parameters <- function(muvec, phivec) {
  q <- length(muvec)
  if(length(phivec)!=q) return("muvec and phivec must be the same length")
  if(muvec[1]!=0) return("muvec[1] must be zero")
  if(phivec[1]!=0 || phivec[q]!=1 || min(diff(phivec))<0) {
    return("phivec not correctly specified")
  }
  return("OK")
}

#' Summary statistics for the ordered stereotype model
#' @param muvec A vector of intensities, with muvec[1]=0 but no other constraints
#' @param phivec A vector of scores, with phivec[1]=0, phivec[q]=1 
#' @param eta A scalar slump value
#' @return the mean, variance and entropy of the distribution
#' @export
summary_stereotype <- function(muvec, phivec, eta=1) {
  msg <- check_stereotype_parameters(muvec,phivec)
  if(msg!="OK") stop(msg)
  if(length(eta)!=1) stop("eta must have length 1")
  q <- length(muvec)
  pvec <- dstereotype(1:q, muvec, phivec, eta)
  meanval <- sum((1:q)*pvec)
  varval <- sum(((1:q)^2)*pvec) - meanval^2
  entropy <- -sum(pvec*log(pvec))
  retval <- c(mean=meanval, var=varval, entropy=entropy)
  return(retval)
}

#####################################################################################
#' Calculate probabilities associated with the ordered stereotype model
#'
#' @param kvec A scalar or vector of integers between 1 and q
#' @param muvec A vector of intensities, with muvec[1]=0 but no other constraints
#' @param phivec A vector of scores, with phivec[1]=0, phivec[q]=1 
#' and phivec[1]<=phivec[2]<= ... <=phivec[q]
#' @param eta A scalar incorporating any covariates and associated effects 
#' (e.g. eta = x^T beta)
#' @param log.p If TRUE then log probabilities will be returned
#' @param .useCpp if TRUE (the default) then a cpp version of the code is used 
#' @return The probability/probabilities that ordered stereotype variable takes 
#' the value(s) \code{k}
#' @details The ordered stereotype model has a q-level ordinal random variable Y with
#' log(P(Y=k|eta)/P(Y=1|eta)) = mu_k + phi_k*eta.
#' @examples
#' dstereotype(1, muvec=c(0,4,2), phivec=c(0,0.8,1), eta=8)
#' 
#' @export
dstereotype <- function(kvec, muvec, phivec, eta, log.p=FALSE, .useCpp=TRUE) {
  msg <- check_stereotype_parameters(muvec,phivec)
  if(msg!="OK") stop(msg)
  if(length(eta)!=1 && length(eta)!=length(kvec)) {
    stop("eta must be of length 1 or the same length as kvec")
  }
  if(.useCpp) return(dstereotype_cpp(kvec, muvec, phivec, eta, log.p))
  q <- length(muvec)
  if(length(eta)==1) {
     omega <- exp(muvec + phivec*eta)
     pvec <- omega/sum(omega)
     retval <- pvec[kvec]
     if(log.p) retval <- log(retval)
     return(retval)
  } else if(length(eta)==length(kvec)) {
     n <- length(eta)
     omega <- exp(muvec + outer(phivec,eta))  # object is q x n
     pmat <- apply(omega,2,function(x) x/sum(x))  # object is q x n
     retval <- pmat[cbind(kvec,1:n)]
     if(log.p) retval <- log(retval)
     return(retval)
  } else {
    stop("eta must be of length 1 or the same length as kvec")
  }
}

#####################################################################################
#' Calculate cumulative probabilities associated with the ordered stereotype model
#'
#' @param kvec A scalar or vector of integers between 1 and q
#' @param muvec A vector of intensities, with muvec[1]=0 but no other constraints
#' @param phivec A vector of scores, with phivec[1]=0, phivec[q]=1 
#' and phivec[1]<=phivec[2]<= ... <=phivec[q]
#' @param eta A scalar incorporating any covariates and associated effects 
#' (e.g. eta = x^T beta)
#' @param lower.tail if TRUE then cumulative probabilities are P(Y<=kvec) otherwise
#' P(Y>kvec) 
#' @param log.p If TRUE then log cumulative probabilities will be returned
#' @param .useCpp if TRUE (the default) then a cpp version of the code is used 
#' @return The probability/probabilities that ordered stereotype variable takes 
#' the value(s) \code{k}
#' @details The ordered stereotype model has a q-level ordinal random variable Y with
#' log(P(Y=k|eta)/P(Y=1|eta)) = mu_k + phi_k*eta.  This function returns Prob(Y<=k|eta).
#' @examples
#' pstereotype(1, muvec=c(0,4,2), phivec=c(0,0.8,1), eta=8)
#' 
#' @export
pstereotype <- function(kvec, muvec, phivec, eta, lower.tail=TRUE, log.p=FALSE,
                        .useCpp=TRUE) {
  msg <- check_stereotype_parameters(muvec,phivec)
  if(msg!="OK") stop(msg)
  if(length(eta)!=1 && length(eta)!=length(kvec)) {
    stop("eta must be of length 1 or the same length as kvec")
  }
  if(.useCpp) return(pstereotype_cpp(kvec, muvec, phivec, eta, lower.tail, log.p))
  q <- length(muvec)
  if(length(eta)==1) {
    omega <- exp(muvec + phivec*eta)
    pvec <- cumsum(omega/sum(omega))
    retval <- pvec[kvec]
    if(!lower.tail) retval <- 1-retval
    if(log.p) retval <- log(retval)
    return(retval)
  } else if(length(eta)==length(kvec)) {
    n <- length(eta)
    omega <- exp(muvec + outer(phivec,eta))  # object is q x n
    pmat <- apply(omega,2,function(x) cumsum(x)/sum(x))  # object is q x n
    retval <- pmat[cbind(kvec,1:n)]
    if(!lower.tail) retval <- 1-retval
    if(log.p) retval <- log(retval)
    return(retval)
  } else {
    stop("eta must be of length 1 or the same length as kvec")
  }
}

#####################################################################################
#' Calculate quantiles associated with the ordered stereotype model
#'
#' @param p A scalar or vector of cumulative probabilities between 0 and 1
#' @param muvec A vector of intensities, with muvec[1]=0 but no other constraints
#' @param phivec A vector of scores, with phivec[1]=0, phivec[q]=1 
#' and phivec[1]<=phivec[2]<= ... <=phivec[q]
#' @param eta A scalar incorporating any covariates and associated effects 
#' (e.g. eta = x^T beta) 
#' @param lower.tail if TRUE then quantiles k satisfy P(Y<=k) otherwise
#' P(Y>k) 
#' @param log.p If TRUE then p is treated as (a vector of) log probabilities
#' @param .useCpp if TRUE (the default) then a cpp version of the code is used 
#' @return The quantiles k associated with the cumulative probabilities p 
#' @details The ordered stereotype model has a q-level ordinal random variable Y with
#' log(P(Y=k|eta)/P(Y=1|eta)) = mu_k + phi_k*eta.  This function returns k where
#' k = 1+argmax(j) Prob(Y<=j|eta) < p
#' @examples
#' qstereotype(0.4, muvec=c(0,4,2), phivec=c(0,0.8,1), eta=8)
#' 
#' @export
qstereotype <- function(p, muvec, phivec, eta, lower.tail=TRUE, log.p=FALSE,
                        .useCpp=TRUE) {
  msg <- check_stereotype_parameters(muvec,phivec)
  if(msg!="OK") stop(msg)
  if(length(eta)!=1 && length(eta)!=length(p)) {
    stop("eta must be of length 1 or the same length as p")
  }
  if(.useCpp) return(qstereotype_cpp(p, muvec, phivec, eta, lower.tail, log.p))
  q <- length(muvec)
  if(log.p) p <- exp(p)
  if(!lower.tail) p <- 1-p
  p <- pmin(1,p)
  if(length(eta)==1) {
     omega <- exp(muvec + phivec*eta)
     pvec <- cumsum(omega/sum(omega))
     pvec[q] <- 1
     kvec <- pmin(q,apply(outer(p,pvec,"<="),1,which.max))
     return(kvec)
  } else if(length(eta)==length(p)) {
     n <- length(eta)
     omega <- exp(muvec + outer(phivec,eta))  # object is q x n
     pmat <- apply(omega,2,function(x) cumsum(x)/sum(x))  # object is q x n
     pmat[q,] <- 1
     dmat <- t(pmat)-p>0 # object is n x q
     kvec <- pmin(q,apply(dmat,1,which.max))
     return(kvec)
  } else {
     stop("eta must be of length 1 or the same length as p")
  }
}

#####################################################################################
#' Generate random draws from the ordered stereotype model
#'
#' @param n Number of draws required
#' @param muvec A vector of intensities, with muvec[1]=0 but no other constraints
#' @param phivec A vector of scores, with phivec[1]=0, phivec[q]=1 
#' and phivec[1]<=phivec[2]<= ... <=phivec[q]
#' @param eta A scalar or vector of length n incorporating any covariates and 
#' associated effects (e.g. eta = x^T beta) 
#' @param .useCpp if TRUE (the default) then a cpp version of the code is used 
#' @return The n draws 
#' @details The ordered stereotype model has a q-level ordinal random variable Y with
#' log(P(Y=k|eta)/P(Y=1|eta)) = mu_k + phi_k*eta.  This function generates
#' from this distribution.
#' @examples
#' rstereotype(4, muvec=c(0,4,2), phivec=c(0,0.8,1), eta=8)
#' rstereotype(4, muvec=c(0,4,2), phivec=c(0,0.8,1), eta=c(1,2,3,8)
#' 
#' @export
rstereotype <- function(n, muvec, phivec, eta, .useCpp=TRUE) {
  msg <- check_stereotype_parameters(muvec,phivec)
  if(msg!="OK") stop(msg)
  if(length(eta)!=1 && length(eta)!=n) {
    stop("eta must be of length 1 or n")
  }
  if(.useCpp) return(rstereotype_cpp(n, muvec, phivec, eta))
  q <- length(muvec)
  pvec <- runif(n)
  kvec <- qstereotype(pvec, muvec, phivec, eta, .useCpp=.useCpp)
  return(kvec)
}

#####################################################################################
logit <- function(p) log(p/(1-p))
expit <- function(x) 1/(1+exp(-x))
#####################################################################################
# convert omega to probability
#' @export
omega.as.probability <- function(x) {
  x <- x-max(x)
  p <- exp(x)
  return(p/sum(p))
}
# convert probability to omega
#' @export
probability.as.omega <- function(p) {
  omega <- log(p/p[1])
  return(omega)
}
#####################################################################################
# Score parameter transformation
# Version 1
phivec.to.uvec <- function(phivec) {
  q <- length(phivec)
  uvec <- logit((phivec[2:(q-1)]-phivec[1:(q-2)])/(1-phivec[1:(q-2)]))
  return(uvec)
}
uvec.to.phivec <- function(uvec) {
  q <- length(uvec)+2
  phivec <- c(rep(0,q-1),1)
  for(k in 2:(q-1)) {
    phivec[k] <- (1 + phivec[k-1]*exp(-uvec[k-1]))/(1 + exp(-uvec[k-1]))
  }
  return(phivec)
}
# Version 2
phivec.to.wvec <- function(phivec) {
  q <- length(phivec)
  nuvec <- logit(phivec[-c(1,q)])
  wvec <- c(nuvec[1], log(diff(nuvec)))
  return(wvec)
}
wvec.to.phivec <- function(wvec) {
  q <- length(wvec)+2
  phivec <- c(rep(0,q-1),1)
  if(q>2) {
    phivec[2] <- expit(wvec[1])
    if(q>3) {
      evec <- cumsum(exp(wvec[-1]))
      for(k in 3:(q-1)) {
        phivec[k] <- expit(wvec[1] + evec[k-2])
      }
    }
  }
  return(phivec)
}
#####################################################################################
# Simulate from the model muvec and phivec are the standard parameters
# and xmat is an n x p design matrix
sim.ostereotype <- function(muvec, phivec, betavec, xmat, .useCpp=TRUE) {
  msg <- check_stereotype_parameters(muvec,phivec)
  if(msg!="OK") stop(msg)
  if(.useCpp) return(sim_ostereotype_cpp(muvec, phivec, betavec, xmat))
  n <- nrow(xmat) # number of observations
  p <- ncol(xmat) # number of predictors
  eta <- as.vector(xmat%*%betavec)
  yvec <- rstereotype(n, muvec, phivec, eta, .useCpp=.useCpp)
  return(yvec)
}

# Log likelihood for a set of observations
llike.ostereotype <- function(kvec, muvec, phivec, betavec, xmat, .useCpp=TRUE) {
   if(.useCpp) return(llike_ostereotype_cpp(kvec, muvec, phivec, betavec, xmat))
   eta <- xmat%*%betavec
   llvals <- dstereotype(kvec, muvec, phivec, eta, log.p=TRUE, .useCpp=.useCpp)
   retval <- sum(llvals)
   return(retval)
}
# Function to be optimised using optim()
load.optfunc.parvec <- function(q, muvec, phivec, betavec) {
  # Put the parameters into a single vector for optim() to use
  return(c(muvec[-1], phivec.to.wvec(phivec), betavec))
}
unload.optfunc.parvec <- function(q, p, parvec) {
  # Extract the parameters from the parameter vector
  muvec <- c(0,parvec[1:(q-1)]) # q-1
  wvec <- parvec[(q-1)+(1:(q-2))] # q-1 + q-2 = 2q-3
  phivec <- wvec.to.phivec(wvec)
  betavec <- parvec[2*q-3 + (1:p)] # 2q-3+p
  return(list(muvec=muvec, phivec=phivec, betavec=betavec))
}
optfunc.llike.ostereotype <- function(parvec, q, kvec, xmat, 
                                      verbose=FALSE, .useCpp=TRUE) {
   # Call the log likelihood function
   if(.useCpp) {
     retval <- optfunc_llike_ostereotype_cpp(parvec, q, kvec, xmat)
   } else {
     n <- length(kvec)
     p <- ncol(xmat)
     parlist <- unload.optfunc.parvec(q, p, parvec)
     retval <- llike.ostereotype(kvec, parlist$muvec, parlist$phivec, parlist$betavec, 
                               xmat, .useCpp)
   }
   if(verbose) print(c(parvec,retval))
   return(retval)
}
fit.ostereotype <- function(q, kvec, xmat, verbose=FALSE, twice=TRUE, parvec0=NULL, 
                            method="BFGS", control=list(maxit=100),
                            .useCpp=TRUE) {
  # estimate the parameters of an ordered stereotype model
  control$fnscale <- -1
  n <- nrow(xmat)
  p <- ncol(xmat)
  if(is.null(parvec0)) {
     tt <- table(kvec)
     fk <- rep(0,q)
     names(fk) <- 1:q
     fk[names(tt)] <- tt
     if(fk[1]>0) {
        muvec0 <- log(fk/fk[1])
     } else {
        muvec0 <- rep(0,q)
     }
     phivec0 <- (0:(q-1))/(q-1)
     betavec0 <- rep(0,p)
     parvec0 <- load.optfunc.parvec(q, muvec0, phivec0, betavec0)
  }
  if(twice) {
    o1 <- optim(parvec0, optfunc.llike.ostereotype,
                q=q, kvec=kvec, xmat=xmat, .useCpp=.useCpp, verbose=verbose,
                method=method, control=control)
    parvec0 <- o1$par
  }
  o2 <- optim(parvec0, optfunc.llike.ostereotype,
              q=q, kvec=kvec, xmat=xmat, .useCpp=.useCpp, verbose=verbose,
              method=method, control=control,
              hessian=TRUE)
  o2$se.par <- try(sqrt(diag(solve(-o2$hessian))))
  o2$estimates <- unload.optfunc.parvec(q, p, o2$par)
  o2$par0 <- parvec0
  o2$value0 <- optfunc.llike.ostereotype(parvec0, q, kvec, xmat, verbose, .useCpp)
  return(o2)
}

#####################################################################################
