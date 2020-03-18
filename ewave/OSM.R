setwd("~/Documents/STAT487/reformatted")


# Load the libraries

library(tidyr)
library(tidyverse)
library(clustord)


# Read in the files

## ewaveOrdinal vs. ewaveBinary
# (n x m) matrix
ewave.ordinal <- read.csv("ewaveOrdinal.csv", header = TRUE)  # n=76, m=235
binary.ewave <- read.csv("ewaveBinary.csv", header = TRUE)    # n=76, m=235



# Process the data

## Ordinal data
## 1st case
ewave.data <- function(ewave.data.input){
  ewave.data.input <- gather(data = ewave.data.input, 
                             key = COL, 
                             value = Y, 
                             E001:E235,
                             factor_key = TRUE)
  names(ewave.data.input)[1] <- "ROW"
  ewave.data.input$COL <- str_remove(ewave.data.input$COL, "E")
  ewave.data.input[c("ROW","COL","Y")] <- lapply(ewave.data.input[c("ROW","COL","Y")], factor)
  ewave.data.input <- ewave.data.input[complete.cases(ewave.data.input),]
  return(ewave.data.input)
}

ewave.ordinal <- ewave.data(ewave.ordinal)
binary.ewave <- ewave.data(binary.ewave)





# Ordered Stereotype Model for ordinal ewave data : Row clustering : Y~row

## 2 row clusters
init.OSM.2r3 <- c(mu2=0.15,mu3=-0.15,mu4=0.15,u2=0.25,u3=0.5,
                  alpha1=3.5) 

OSM.2rows.cluster3 <- rowclustering(formula = "Y~row", 
                                    model = "OSM", 
                                    nclus.row = 2, 
                                    long.df = ewave.ordinal,
                                    initvect = init.OSM.2r3)
# OSM.2rows.cluster2$criteria 


## 3 row clusters
init.OSM.3r3 <- c(mu2=0.15,mu3=-0.15,mu4=0.15,u2=0.25,u3=0.5,
                  alpha1=3.5,alpha2=3.0) 

OSM.3rows.cluster3 <- rowclustering(formula = "Y~row", 
                                    model = "OSM", 
                                    nclus.row = 3, 
                                    long.df = ewave.ordinal,
                                    initvect = init.OSM.3r3)
OSM.3rows.cluster3$criteria


# 4 row clusters
init.OSM.4r3 <- c(mu2=0.15,mu3=-0.15,mu4=0.15,u2=0.25,u3=0.5,
                  alpha1=3.5,alpha2=3.0,alpha3=2.5) 
OSM.4rows.cluster3 <- rowclustering(formula = "Y~row", 
                                    model = "OSM", 
                                    nclus.row = 4, 
                                    long.df = ewave.ordinal,
                                    initvect = init.OSM.4r3)
# OSM.4rows.cluster2$criteria


# 5 row clusters
init.OSM.5r3 <- c(mu2=0.15,mu3=-0.15,mu4=0.15,u2=0.25,u3=0.5,
                  alpha1=3.5,alpha2=3.0,alpha3=2.5,alpha4=4.5) 
OSM.5rows.cluster3 <- rowclustering(formula = "Y~row", 
                                    model = "OSM", 
                                    nclus.row = 5, 
                                    long.df = ewave.ordinal,
                                    initvect = init.OSM.5r3)
# OSM.5rows.cluster2$criteria




#Ordered Stereotype Model: Column Clustering : Y~column

## 2 column clusters
initvect3.2col <- c(mu2=0.15,mu3=-0.15,mu4=0.15,u2=0.25,u3=0.5,
                    beta1=3.5,alphai = rnorm(75)) 
OSM3.2col <- columnclustering(formula = "Y~row+column",
                               model = "OSM",
                               nclus.column = 2,
                               long.df = ewave.ordinal,
                               initvect = initvect3.2col)
OSM2.2col$criteria


## 3 column clusters

initvect3.3col <- c(mu2=0.15,mu3=-0.15,mu4=0.15,u2=0.25,u3=0.5,
                    beta1=3.5,beta2=3.0) 

OSM3.3col <- columnclustering(formula = "Y~column",
                               model = "OSM",
                               nclus.column = 3,
                               long.df = ewave.ordinal,
                               initvect = initvect3.3col)

#OSM2.3col$criteria


## 4 column clusters

initvect3.4col <- c(mu2=0.15,mu3=-0.15,mu4=0.15,u2=0.25,u3=0.5,
                    beta1=3.5,beta2=3.0,beta3=2.5) 

OSM3.4col <- columnclustering(formula = "Y~column",
                               model = "OSM",
                               nclus.column = 4,
                               long.df = ewave.ordinal,
                               initvect = initvect3.4col)

# OSM2.4col$criteria


## 5 column clusters

initvect3.5col <- c(mu2=0.15,mu3=-0.15,mu4=0.15,u2=0.25,u3=0.5,
                    beta1=3.5,beta2=3.0,beta3=2.5,beta4=4.5) 

OSM3.5col <- columnclustering(formula = "Y~column",
                              model = "OSM",
                              nclus.column = 5,
                              long.df = ewave.ordinal,
                              initvect = initvect3.5col)

# OSM2.5col$criteria


## 6 columns

initvect3.6col <- c(mu2=0.15,mu3=-0.15,mu4=0.15,u2=0.25,u3=0.5,
                    beta1=3.5,beta2=3.0,beta3=2.5,beta4=4.5,beta5=-2.0) 
OSM3.6col <- columnclustering(formula = "Y~column",
                              model = "OSM",
                              nclus.column = 6,
                              long.df = ewave.ordinal,
                              initvect = initvect3.6col)
# OSM2.6col$criteria


## 7 column clusters
initvect3.7col <- c(mu2=0.15,mu3=-0.15,mu4=0.15,u2=0.25,u3=0.5,
                    beta1=3.5,beta2=3.0,beta3=2.5,beta4=4.5,beta5=-2.0,beta6=-5.0) 

OSM3.7col <- columnclustering(formula = "Y~column",
                              model = "OSM",
                              nclus.column = 7,
                              long.df = ewave.ordinal,
                              initvect = initvect3.7col)
# OSM3.2col$criteria
# OSM3.3col$criteria
# OSM3.4col$criteria
# OSM3.5col$criteria
# OSM3.6col$criteria
# OSM3.7col$criteria
OSM.2rows.cluster3$EM.status$best.lli
OSM.3rows.cluster3$EM.status$best.lli
OSM.4rows.cluster3$EM.status$best.lli
OSM.5rows.cluster3$EM.status$best.lli

OSM.2rows.cluster3$criteria$AIC
OSM.3rows.cluster3$criteria$AIC
OSM.4rows.cluster3$criteria$AIC
OSM.5rows.cluster3$criteria$AIC


OSM3.2col$EM.status$best.lli
OSM3.3col$EM.status$best.lli
OSM3.4col$EM.status$best.lli
OSM3.5col$EM.status$best.lli
OSM3.6col$EM.status$best.lli
OSM3.7col$EM.status$best.lli

OSM3.2col$criteria
OSM3.3col$criteria
OSM3.4col$criteria
OSM3.5col$criteria
OSM3.6col$criteria
OSM3.7col$criteria