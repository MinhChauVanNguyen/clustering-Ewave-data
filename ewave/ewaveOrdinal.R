
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
init.OSM.2r <- c(mu2=0.4,mu3=0.6,mu4=0.8,u2=1.8,u3=2,
                 alpha1=0.4) 
init.OSM.2r2 <- c(mu2=0.25,mu3=0.1,mu4=1.1,u2=0.7,u3=-0.5,
                 alpha1=-2) 
init.OSM.2r3 <- c(mu2=0.15,mu3=-0.15,mu4=0.15,u2=0.25,u3=0.25,
                  alpha1=3.5) 

OSM.2rows.cluster2 <- rowclustering(formula = "Y~row", 
                                   model = "OSM", 
                                   nclus.row = 2, 
                                   long.df = ewave.ordinal,
                                   initvect = init.OSM.2r2)
minh <- as.matrix(OSM.2rows.cluster2$RowClusters[[1]])
connor <- as.matrix(OSM.2rows.cluster2$RowClusters[[1]])
sex <- cbind(minh,connor)
sex

## 3 row clusters
init.OSM.3r <- c(mu2=0.4,mu3=0.6,mu4=0.8,u2=1.8,u3=2,
                 alpha1=0.4,alpha2=0.6)
init.OSM.3r2 <- c(mu2=0.25,mu3=0.1,mu4=1.1,u2=0.7,u3=-0.5,
                  alpha1=-2,alpha2=0.08) 
init.OSM.3r3 <- c(mu2=0.15,mu3=-0.15,mu4=0.15,u2=0.25,u3=0.25,
                  alpha1=3.5,alpha2=3.0) 

OSM.3rows.cluster3 <- rowclustering(formula = "Y~row", 
                                   model = "OSM", 
                                   nclus.row = 3, 
                                   long.df = ewave.ordinal,
                                   initvect = init.OSM.3r3)
OSM.3rows.cluster3$criteria


# 4 row clusters
init.OSM.4r <- c(mu2=0.4,mu3=0.6,mu4=0.8,u2=1.8,u3=2,alpha1=0.4,alpha2=0.6,alpha3=0.8)
init.OSM.4r2 <- c(mu2=0.25,mu3=0.1,mu4=1.1,u2=0.7,u3=-0.5,
                  alpha1=-2,alpha2=0.08,alpha3=0.01) 
init.OSM.4r3 <- c(mu2=0.15,mu3=-0.15,mu4=0.15,u2=0.25,u3=0.25,
                  alpha1=3.5,alpha2=3.0,alpha3=2.5) 
OSM.4rows.cluster3 <- rowclustering(formula = "Y~row", 
                                   model = "OSM", 
                                   nclus.row = 4, 
                                   long.df = ewave.ordinal,
                                   initvect = init.OSM.4r3)
# OSM.4rows.cluster2$criteria


# 5 row clusters
init.OSM.5r <- c(mu2=0.4,mu3=0.6,mu4=0.8,u2=1.8,u3=2,alpha1=0.4,alpha2=0.6,alpha3=0.8, alpha4=1.0)  
init.OSM.5r2 <- c(mu2=0.25,mu3=0.1,mu4=1.1,u2=0.7,u3=-0.5,
                  alpha1=-2,alpha2=0.08,alpha3=0.01,alpha4=-0.15)
init.OSM.5r3 <- c(mu2=0.15,mu3=-0.15,mu4=0.15,u2=0.25,u3=0.25,
                  alpha1=3.5,alpha2=3.0,alpha3=2.5,alpha4=4.5) 
OSM.5rows.cluster2 <- rowclustering(formula = "Y~row", 
                                   model = "OSM", 
                                   nclus.row = 5, 
                                   long.df = ewave.ordinal,
                                   initvect = init.OSM.5r2)
OSM.5rows.cluster2$RowClusters




#Ordered Stereotype Model: Column Clustering : Y~column

## 2 column clusters
initvect1.2col <- c(mu2=0.4,mu3=0.6,mu4=0.8,u2=1.8,u3=2,beta1=0.4)
initvect2.2col <- c(mu2=0.25,mu3=0.1,mu4=1.1,u2=0.7,u3=-0.5,
                    beta1=-2)
initvect3.2col <- c(mu2=0.15,mu3=-0.15,mu4=0.15,u2=0.25,u3=0.25,
                    beta1=3.5) 
OSM2.2col3 <- columnclustering(formula = "Y~column",
                             model = "OSM",
                             nclus.column = 2,
                             long.df = ewave.ordinal,
                             initvect = initvect3.2col)
OSM2.2col$criteria


## 3 column clusters
initvect1.3col <- c(mu2=0.4,mu3=0.6,mu4=0.8,u2=1.8,u3=2,
                    beta1=0.4,beta2=0.6)
initvect2.3col <- c(mu2=0.25,mu3=0.1,mu4=1.1,u2=0.7,u3=-0.5,
                    beta1=-2,beta2=0.08)

initvect3.3col <- c(mu2=0.15,mu3=-0.15,mu4=0.15,u2=0.25,u3=0.25,
                    beta1=3.5,beta2=3.0) 

OSM2.3col3 <- columnclustering(formula = "Y~column",
                             model = "OSM",
                             nclus.column = 3,
                             long.df = ewave.ordinal,
                             initvect = initvect3.3col)

#OSM2.3col$criteria


## 4 column clusters
initvect1.4col <- c(mu2=0.4,mu3=0.6,mu4=0.8,u2=1.8,u3=2,
                    beta1=0.4,beta2=0.6,beta3=0.8)
initvect2.4col <- c(mu2=0.25,mu3=0.1,mu4=1.1,u2=0.7,u3=-0.5,
                    beta1=-2,beta2=0.08,beta3=0.01)

initvect3.4col <- c(mu2=0.15,mu3=-0.15,mu4=0.15,u2=0.25,u3=0.25,
                    beta1=3.5,beta2=3.0,beta3=2.5) 

OSM2.4col3 <- columnclustering(formula = "Y~column",
                             model = "OSM",
                             nclus.column = 4,
                             long.df = ewave.ordinal,
                             initvect = initvect3.4col)

# OSM2.4col$criteria


## 5 column clusters
initvect1.5col <- c(mu2=0.4,mu3=0.6,mu4=0.8,u2=1.8,u3=2,
                    beta1=0.4, beta2=0.6,beta3=0.8,beta4=1.0)
initvect2.5col <- c(mu2=0.25,mu3=0.1,mu4=1.1,u2=0.7,u3=-0.5,
                    beta1=-2,beta2=0.08,beta3=0.01,beta4=-0.15)

initvect3.5col <- c(mu2=0.15,mu3=-0.15,mu4=0.15,u2=0.25,u3=0.25,
                    beta1=3.5,beta2=3.0,beta3=2.5,beta4=4.5) 

OSM2.5col <- columnclustering(formula = "Y~column",
                             model = "OSM",
                             nclus.column = 5,
                             long.df = ewave.ordinal,
                             initvect = initvect2.5col)

# OSM2.5col$criteria


## 6 columns
initvect1.6col <- c(mu2=0.4,mu3=0.6,mu4=0.8,u2=1.5,u3=2,
                    beta1=0.4,beta2=0.6,beta3=0.8, beta4=1.0,beta5=1.2)
initvect2.6col <- c(mu2=0.25,mu3=0.1,mu4=1.1,u2=0.7,u3=-0.5,
                    beta1=-2,beta2=0.08,beta3=0.01,beta4=-0.15,beta5=0.64)
initvect3.6col <- c(mu2=0.15,mu3=-0.15,mu4=0.15,u2=0.25,u3=0.25,
                    beta1=3.5,beta2=3.0,beta3=2.5,beta4=4.5,beta5=-2.0) 
OSM2.6col <- columnclustering(formula = "Y~column",
                            model = "OSM",
                            nclus.column = 6,
                            long.df = ewave.ordinal,
                            initvect = initvect2.6col)
# OSM2.6col$criteria


## 7 column clusters
initvect1.7col <- c(mu2=0.4,mu3=0.6,mu4=0.8,u2=1.5,u3=2,
                    beta1=0.4,beta2=0.6,beta3=0.8, beta4=1.0,beta5=1.2,beta6=1.4)
initvect2.7col <- c(mu2=0.25,mu3=0.1,mu4=1.1,u2=0.7,u3=-0.5,
                    beta1=-2,beta2=0.08,beta3=0.01,beta4=-0.15,beta5=0.64,beta6=-1.2)
initvect3.7col <- c(mu2=0.15,mu3=-0.15,mu4=0.15,u2=0.25,u3=0.25,
                    beta1=3.5,beta2=3.0,beta3=2.5,beta4=4.5,beta5=-2.0,beta6=-5.0) 

OSM2.7col <- columnclustering(formula = "Y~column",
                              model = "OSM",
                              nclus.column = 7,
                              long.df = ewave.ordinal,
                              initvect = initvect2.7col)
OSM2.7col$ColumnCluster




# Proportional Odds Model for ordinal ewave data : Row clustering : Y~row

## 2 row clusters
PO.initvect1.2row <- c(mu2=0.4,mu3=0.6,mu4=0.8,alpha1=0.4)
PO.initvect2.2row <- c(mu2=0.25,mu3=0.1,mu4=1.1,
                       alpha1=-2)
PO.initvect3.2row <- c(mu2=0.15,mu3=-0.15,mu4=0.15,
                       alpha1=0.002) 

PO.tworow <- rowclustering(formula = "Y~row", 
                           model = "POM", 
                           nclus.row = 2, 
                           long.df = ewave.ordinal,
                           initvect = PO.initvect1.2row)



## 3 row clusters
PO.initvect1.3row <- c(mu2=0.4,mu3=0.6,mu4=0.8,alpha1=0.4,alpha2=0.6)
PO.initvect2.3row <- c(mu2=0.25,mu3=0.1,mu4=1.1,
                       alpha1=-2,alpha2=0.08)
PO.initvect3.3row <- c(mu2=0.15,mu3=-0.15,mu4=0.15,
                       alpha1=0.002,alpha2=0.002) 

PO.threerow <- rowclustering(formula = "Y~row", 
                             model = "POM", 
                             nclus.row = 3, 
                             long.df = ewave.ordinal,
                             initvect = PO.initvect1.3row)

## 4 row clusters
PO.initvect1.4row <- c(mu2=0.4,mu3=0.6,mu4=0.8,
                       alpha1=0.4,alpha2=0.6,alpha3=0.8)
PO.initvect2.4row <- c(mu2=0.25,mu3=0.1,mu4=1.1,
                        alpha1=-2,alpha2=0.08,alpha3=0.01)
PO.initvect3.4row <- c(mu2=0.15,mu3=-0.15,mu4=0.15,
                        alpha1=0.002,alpha2=0.002,alpha3=0.0015) 
PO.fourrow <- rowclustering(formula = "Y~row", 
                            model = "POM", 
                            nclus.row = 4, 
                            long.df = ewave.ordinal,
                            initvect = PO.initvect1.4row)


## 5 row clusters
PO.initvect1.5row <- c(mu2=0.4,mu3=0.6,mu4=0.8,
                       alpha1=0.4,alpha2=0.6,alpha3=0.8,alpha4=1)
PO.initvect2.5row <- c(mu2=0.25,mu3=0.1,mu4=1.1,
                  alpha1=-2,alpha2=0.08,alpha3=0.01,alpha4=-0.15)
PO.initvect3.5row <- c(mu2=0.15,mu3=-0.15,mu4=0.15,
                  alpha1=0.002,alpha2=0.002,alpha3=0.0015,alpha4=0.0015) 
PO.fiverow <- rowclustering(formula = "Y~row", 
                            model = "POM", 
                            nclus.row = 5, 
                            long.df = ewave.ordinal,
                            initvect = PO.initvect1.5row)

PO.fiverow$RowClusters
PO2.row <- PO.tworow$EM.status$best.lli
PO3.row <- PO.threerow$EM.status$best.lli
PO4.row <- PO.fourrow$EM.status$best.lli
PO5.row <- PO.fiverow$EM.status$best.lli



#Proportional Odds Model : Column clustering : Y~column

## 2 column clusters
PO.initvect1.2col <- c(mu2=0.4,mu3=0.6,mu4=0.8,
                    beta1=0.4)
PO.initvect2.2col <- c(mu2=0.25,mu3=0.1,mu4=1.1,
                    beta1=-2)
PO.initvect3.2col <- c(mu2=0.15,mu3=-0.15,mu4=0.15,
                    beta1=0.002)
PO.twocol <- columnclustering(formula = "Y~column",
                               model = "POM",
                               nclus.column = 2,
                               long.df = ewave.ordinal,
                               initvect = PO.initvect1.2col)




## 3 column clusters
PO.initvect1.3col <- c(mu2=0.4,mu3=0.6,mu4=0.8,
                       beta1=0.4,beta2=0.6)
PO.initvect2.3col <- c(mu2=0.25,mu3=0.1,mu4=1.1,
                       beta1=-2,beta2=0.08)
PO.initvect3.3col <- c(mu2=0.15,mu3=-0.15,mu4=0.15,
                       beta1=0.002,beta2=0.002)
PO.threecol <- columnclustering(formula = "Y~column",
                                 model = "POM",
                                 nclus.column = 3,
                                 long.df = ewave.ordinal,
                                 initvect = PO.initvect1.3col)




## 4 column clusters
PO.initvect1.4col <- c(mu2=0.4,mu3=0.6,mu4=0.8,
                       beta1=0.4,beta2=0.6,beta3=0.8)
PO.initvect2.4col <- c(mu2=0.25,mu3=0.1,mu4=1.1,
                       beta1=-2,beta2=0.08,beta3=0.01)
PO.initvect3.4col <- c(mu2=0.15,mu3=-0.15,mu4=0.15,
                       beta1=0.002,beta2=0.002,beta3=0.0015)
PO.fourcol <- columnclustering(formula = "Y~column",
                                model = "POM",
                                nclus.column = 4,
                                long.df = ewave.ordinal,
                                initvect = PO.initvect1.4col)


## 5 column clusters
PO.initvect1.5col <- c(mu2=0.4,mu3=0.6,mu4=0.8,
                       beta1=0.4,beta2=0.6,beta3=0.8,beta4=1.0)
PO.initvect2.5col <- c(mu2=0.25,mu3=0.1,mu4=1.1,
                       beta1=-2,beta2=0.08,beta3=0.01,beta4=-0.15)
PO.initvect3.5col <- c(mu2=0.15,mu3=-0.15,mu4=0.15,
                       beta1=0.002,beta2=0.002,beta3=0.0015,beta4=0.0015)
PO.fivecol <- columnclustering(formula = "Y~column",
                                model = "POM",
                                nclus.column = 5,
                                long.df = ewave.ordinal,
                                initvect = PO.initvect1.5col)



## 6 column clusters
PO.initvect1.6col <- c(mu2=0.4,mu3=0.6,mu4=0.8,beta1=0.4,beta2=0.6,
                       beta3=0.8,beta4=1.0,beta5=1.2)
PO.initvect2.6col <- c(mu2=0.25,mu3=0.1,mu4=1.1,beta1=-2,beta2=0.08,
                       beta3=0.01,beta4=-0.15,beta5=0.64)
PO.initvect3.6col <- c(mu2=0.15,mu3=-0.15,mu4=0.15,beta1=0.002,beta2=0.002,
                       beta3=0.0015,beta4=0.0015,beta5=0.0001)
PO.sixcol <- columnclustering(formula = "Y~column",
                               model = "POM",
                               nclus.column = 6,
                               long.df = ewave.ordinal,
                               initvect = PO.initvect1.6col)


## 7 column clusters
PO.initvect1.7col <- c(mu2=0.4,mu3=0.6,mu4=0.8,
                       beta1=0.4,beta2=0.6,beta3=0.8,beta4=1.0,beta5=1.2,beta6=1.4)
PO.initvect2.7col <- c(mu2=0.25,mu3=0.1,mu4=1.1,
                       beta1=-2,beta2=0.08,beta3=0.01,beta4=-0.15,beta5=0.64,beta6=-1.2)
PO.initvect3.7col <- c(mu2=0.15,mu3=-0.15,mu4=0.15,
                       beta1=0.002,beta2=0.002,beta3=0.0015,beta4=0.0015,beta5=0.0001,beta6=0.0001)
PO.sevencol <- columnclustering(formula = "Y~column",
                                 model = "POM",
                                 nclus.column = 7,
                                 long.df = ewave.ordinal,
                                 initvect = PO.initvect1.7col)
PO.sevencol$ColumnClusters
PO.twocol$criteria$AIC
PO.threecol$criteria$AIC
PO.fourcol$criteria$AIC
PO.fivecol$criteria$AIC
PO.sixcol$criteria$AIC
PO.sevencol$criteria$AIC



# Binary Model : Row Clustering : Y ~ row
## 2 row clusters
initvect1.2bin <- c(mu=0.6, alpha1=0.4)
initvect2.2bin <- c(mu=-1.2, alpha1=-2)
initvect3.2bin <- c(mu=0.0001, alpha1=0.002)

tworow.binary <- rowclustering(formula = "Y~row",
                               model = "Binary",
                               nclus.row = 2,
                               long.df = binary.ewave,
                               initvect = initvect1.2bin)

## 3 row clusters
initvect1.3bin <- c(mu=0.6, alpha1=0.4, alpha2=0.6)
initvect2.3bin <- c(mu=-1.2, alpha1=-2, alpha2=0.08)
initvect3.3bin <- c(mu=0.0001, alpha1=0.002,alpha2=0.002)
threerow.binary <- rowclustering(formula = "Y~row",
                                 model = "Binary",
                                 nclus.row = 3,
                                 long.df = binary.ewave,
                                 initvect = initvect1.3bin)

## 4 row clusters
initvect1.4bin <- c(mu=0.6, alpha1=0.4, alpha2=0.6, alpha3=0.8)
initvect2.4bin <- c(mu=-1.2, alpha1=-2, alpha2=0.08, alpha3=0.01)
initvect3.4bin <- c(mu=0.0001, alpha1=0.002,alpha2=0.002, alpha3=0.0015)
fourrow.binary <- rowclustering(formula = "Y~row",
                                model = "Binary",
                                nclus.row = 4,
                                long.df = binary.ewave,
                                initvect = initvect1.4bin)


## 5 row clusters
initvect1.5bin <- c(mu=0.6, alpha1=0.4, alpha2=0.6, alpha3=0.8, alpha4=1.0)
initvect2.5bin <- c(mu=-1.2, alpha1=-2, alpha2=0.08, alpha3=0.01, alpha4=-0.15)
initvect3.5bin <- c(mu=0.0001, alpha1=0.002,alpha2=0.002, alpha3=0.0015, alpha4=0.0015)
fiverow.binary <- rowclustering(formula = "Y~row",
                             model = "Binary",
                             nclus.row = 5,
                             long.df = binary.ewave,
                             initvect = initvect1.5bin)
fiverow.binary$RowClusters




# Binary Model for binary ewave data : Column Clustering : Y ~ column

## 2 column clusters
initvect1.2col.bin <- c(mu=0.4,beta1=0.4)
initvect2.2col.bin <- c(mu=0.25,beta1=-2)
initvect3.2col.bin <- c(mu=0.15,beta1=0.002)
twocol.binary <- columnclustering(formula = "Y~column",
                                    model = "Binary",
                                    nclus.column = 2,
                                    long.df = binary.ewave,
                                    initvect = initvect1.2col.bin)


## 3 column clusters
initvect1.3col.bin <- c(mu=0.4,beta1=0.4,beta2=0.6)
initvect2.3col.bin <- c(mu=0.25,beta1=-2,beta2=0.08)
initvect3.3col.bin <- c(mu=0.15,beta1=0.002,beta2=0.002)
threecol.binary <- columnclustering(formula = "Y~column",
                                    model = "Binary",
                                    nclus.column = 3,
                                    long.df = binary.ewave,
                                    initvect = initvect1.3col.bin)

## 4 column clusters
initvect1.4col.bin <- c(mu=0.4,beta1=0.4,beta2=0.6,beta3=0.8)
initvect2.4col.bin <- c(mu=0.25,beta1=-2,beta2=0.08,beta3=0.01)
initvect3.4col.bin <- c(mu=0.15,beta1=0.002,beta2=0.002,beta3=0.0015)
fourcol.binary <- columnclustering(formula = "Y~column",
                                   model = "Binary",
                                   nclus.column = 4,
                                   long.df = binary.ewave,
                                   initvect = initvect1.4col.bin)

## 5 column clusters
initvect1.5col.bin <- c(mu=0.4,beta1=0.4,beta2=0.6,beta3=0.8,beta4=1.0)
initvect2.5col.bin <- c(mu=0.25,beta1=-2,beta2=0.08,beta3=0.01,beta4=-0.15)
initvect3.5col.bin <- c(mu=0.15,beta1=0.002,beta2=0.002,beta3=0.0015,beta4=0.0015)
fivecol.binary <- columnclustering(formula = "Y~column",
                                   model = "Binary",
                                   nclus.column = 5,
                                   long.df = binary.ewave,
                                   initvect = initvect1.5col.bin)



## 6 column clusters
initvect1.6col.bin <- c(mu=0.4,beta1=0.4,beta2=0.6,beta3=0.8,beta4=1.0,beta5=1.2)
initvect2.6col.bin <- c(mu=0.25,beta1=-2,beta2=0.08,beta3=0.01,beta4=-0.15,beta5=0.64)
initvect3.6col.bin <- c(mu=0.15,beta1=0.002,beta2=0.002,beta3=0.0015,beta4=0.0015,beta5=0.0001)
sixcol.binary <- columnclustering(formula = "Y~column",
                                  model = "Binary",
                                  nclus.column = 6,
                                  long.df = binary.ewave,
                                  initvect = initvect1.6col.bin)


## 7 column clusters
initvect1.7col.bin <- c(mu=0.4,beta1=0.4,beta2=0.6,beta3=0.8,beta4=1.0,beta5=1.2,beta6=1.4)
initvect2.7col.bin <- c(mu=0.25,beta1=-2,beta2=0.08,beta3=0.01,beta4=-0.15,beta5=0.64,beta6=-1.2)
initvect3.7col.bin <- c(mu=0.15,beta1=0.002,beta2=0.002,beta3=0.0015,beta4=0.0015,beta5=0.0001,beta6=0.0001)

sevencol.binary <- columnclustering(formula = "Y~column",
                                    model = "Binary",
                                    nclus.column = 7,
                                    long.df = binary.ewave,
                                    initvect = initvect1.7col.bin)

sevencol.binary$ColumnClusters
twocol.binary$criteria
threecol.binary$criteria
fourcol.binary$criteria
fivecol.binary$criteria
sixcol.binary$criteria
sevencol.binary$criteria


