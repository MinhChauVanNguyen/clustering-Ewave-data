# Binary Model : Row Clustering : Y ~ row


## 2 row clusters
initvect3.2bin <- c(mu=0.0001,alpha1=3.5) 

tworow.binary3 <- rowclustering(formula = "Y~row",
                               model = "Binary",
                               nclus.row = 2,
                               long.df = binary.ewave,
                               initvect = initvect3.2bin)

## 3 row clusters
initvect3.3bin <- c(mu=0.0001,alpha1=3.5,alpha2=3.0)
threerow.binary3 <- rowclustering(formula = "Y~row",
                                 model = "Binary",
                                 nclus.row = 3,
                                 long.df = binary.ewave,
                                 initvect = initvect3.3bin)

## 4 row clusters
initvect3.4bin <- c(mu=0.0001,alpha1=3.5,alpha2=3.0,alpha3=2.5) 
fourrow.binary3 <- rowclustering(formula = "Y~row",
                                model = "Binary",
                                nclus.row = 4,
                                long.df = binary.ewave,
                                initvect = initvect3.4bin)


## 5 row clusters
initvect3.5bin <- c(mu=0.0001, 
                    alpha1=3.5,alpha2=3.0,alpha3=2.5,alpha4=4.5) 

fiverow.binary3 <- rowclustering(formula = "Y~row",
                                model = "Binary",
                                nclus.row = 5,
                                long.df = binary.ewave,
                                initvect = initvect3.5bin)

fiverow.binary3$EM.status$best.lli
fiverow.binary3$criteria

# Binary Model for binary ewave data : Column Clustering : Y ~ column

## 2 column clusters
initvect3.2col.bin <- c(mu=0.15,beta1=3.5) 

twocol.binary3 <- columnclustering(formula = "Y~column",
                                  model = "Binary",
                                  nclus.column = 2,
                                  long.df = binary.ewave,
                                  initvect = initvect3.2col.bin)
twocol.binary3$EM.status$best.lli
twocol.binary3$criteria


## 3 column clusters
initvect3.3col.bin <- c(mu=0.15,beta1=3.5,beta2=3.0) 

threecol.binary3 <- columnclustering(formula = "Y~column",
                                    model = "Binary",
                                    nclus.column = 3,
                                    long.df = binary.ewave,
                                    initvect = initvect3.3col.bin)
threecol.binary3$EM.status$best.lli
threecol.binary3$criteria

## 4 column clusters

initvect3.4col.bin <- c(mu=0.15,
                        beta1=3.5,beta2=3.0,beta3=2.5)
fourcol.binary3 <- columnclustering(formula = "Y~column",
                                   model = "Binary",
                                   nclus.column = 4,
                                   long.df = binary.ewave,
                                   initvect = initvect3.4col.bin)
fourcol.binary3$EM.status$best.lli
fourcol.binary3$criteria


## 5 column clusters
initvect3.5col.bin <- c(mu=0.15,
                        beta1=3.5,beta2=3.0,beta3=2.5,beta4=4.5) 

fivecol.binary3 <- columnclustering(formula = "Y~column",
                                   model = "Binary",
                                   nclus.column = 5,
                                   long.df = binary.ewave,
                                   initvect = initvect3.5col.bin)
fivecol.binary3$EM.status$best.lli
fivecol.binary3$criteria


## 6 column clusters
initvect3.6col.bin <- c(mu=0.15,
                        beta1=3.5,beta2=3.0,beta3=2.5,beta4=4.5,beta5=-2.0) 

sixcol.binary3 <- columnclustering(formula = "Y~column",
                                  model = "Binary",
                                  nclus.column = 6,
                                  long.df = binary.ewave,
                                  initvect = initvect3.6col.bin)
sixcol.binary3$EM.status$best.lli
sixcol.binary3$criteria


## 7 column clusters

initvect3.7col.bin <- c(mu=0.15,
                      beta1=3.5,beta2=3.0,beta3=2.5,beta4=4.5,beta5=-2.0,beta6=-5.0) 

sevencol.binary3 <- columnclustering(formula = "Y~column",
                                    model = "Binary",
                                    nclus.column = 7,
                                    long.df = binary.ewave,
                                    initvect = initvect3.7col.bin)

sevencol.binary3$criteria
sevencol.binary3$EM.status$best.lli
