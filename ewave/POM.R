
# Proportional Odds Model for ordinal ewave data : Row clustering : Y~row

## 2 row clusters
PO.initvect3.2row <- c(mu2=0.15,mu3=-0.15,mu4=0.15,
                       alpha1=3.5) 

PO.tworow3 <- rowclustering(formula = "Y~row", 
                           model = "POM", 
                           nclus.row = 2, 
                           long.df = ewave.ordinal,
                           initvect = PO.initvect3.2row)


## 3 row clusters
PO.initvect3.3row <- c(mu2=0.15,mu3=-0.15,mu4=0.15,
                       alpha1=3.5,alpha2=3.0)

PO.threerow3 <- rowclustering(formula = "Y~row", 
                             model = "POM", 
                             nclus.row = 3, 
                             long.df = ewave.ordinal,
                             initvect = PO.initvect3.3row)

## 4 row clusters
PO.initvect3.4row <- c(mu2=0.15,mu3=-0.15,mu4=0.15,
                       alpha1=3.5,alpha2=3.0,alpha3=2.5) 

PO.fourrow3 <- rowclustering(formula = "Y~row", 
                            model = "POM", 
                            nclus.row = 4, 
                            long.df = ewave.ordinal,
                            initvect = PO.initvect3.4row)


## 5 row clusters
PO.initvect3.5row <- c(mu2=0.15,mu3=-0.15,mu4=0.15,
                       alpha1=3.5,alpha2=3.0,alpha3=2.5,alpha4=4.5) 

PO.fiverow3 <- rowclustering(formula = "Y~row", 
                            model = "POM", 
                            nclus.row = 5, 
                            long.df = ewave.ordinal,
                            initvect = PO.initvect3.5row)

PO.tworow3$EM.status$best.lli
PO.threerow3$EM.status$best.lli
PO.fourrow3$EM.status$best.lli
PO.fiverow3$EM.status$best.lli



#Proportional Odds Model : Column clustering : Y~column

## 2 column clusters
PO.initvect3.2col <- c(mu2=0.15,mu3=-0.15,mu4=0.15,
                       beta1=3.5) 
PO.twocol3 <- columnclustering(formula = "Y~column",
                              model = "POM",
                              nclus.column = 2,
                              long.df = ewave.ordinal,
                              initvect = PO.initvect3.2col)


## 3 column clusters
PO.initvect3.3col <- c(mu2=0.15,mu3=-0.15,mu4=0.15,
                       beta1=3.5,beta2=3.0) 
PO.threecol3 <- columnclustering(formula = "Y~column",
                                model = "POM",
                                nclus.column = 3,
                                long.df = ewave.ordinal,
                                initvect = PO.initvect3.3col)


## 4 column clusters
PO.initvect3.4col <- c(mu2=0.15,mu3=-0.15,mu4=0.15,
                       beta1=3.5,beta2=3.0,beta3=2.5) 
PO.fourcol3 <- columnclustering(formula = "Y~column",
                               model = "POM",
                               nclus.column = 4,
                               long.df = ewave.ordinal,
                               initvect = PO.initvect3.4col)


## 5 column clusters
PO.initvect3.5col <- c(mu2=0.15,mu3=-0.15,mu4=0.15,
                       beta1=3.5,beta2=3.0,beta3=2.5,beta4=4.5) 
PO.fivecol3 <- columnclustering(formula = "Y~column",
                               model = "POM",
                               nclus.column = 5,
                               long.df = ewave.ordinal,
                               initvect = PO.initvect3.5col)


## 6 column clusters
PO.initvect3.6col <- c(mu2=0.15,mu3=-0.15,mu4=0.15,
                       beta1=3.5,beta2=3.0,beta3=2.5,beta4=4.5,beta5=-2.0) 
PO.sixcol3 <- columnclustering(formula = "Y~column",
                              model = "POM",
                              nclus.column = 6,
                              long.df = ewave.ordinal,
                              initvect = PO.initvect3.6col)


## 7 column clusters
PO.initvect3.7col <- c(mu2=0.15,mu3=-0.15,mu4=0.15,
                      beta1=3.5,beta2=3.0,beta3=2.5,beta4=4.5,beta5=-2.0,beta6=-5.0) 

PO.sevencol3 <- columnclustering(formula = "Y~column",
                                model = "POM",
                                nclus.column = 7,
                                long.df = ewave.ordinal,
                                initvect = PO.initvect3.7col)
PO.twocol3$criteria$AIC
PO.threecol3$criteria$AIC
PO.fourcol3$criteria$AIC
PO.fivecol3$criteria$AIC
PO.sixcol3$criteria$AIC
PO.sevencol3$criteria$AIC

PO.twocol3$EM.status$best.lli
PO.threecol3$criteria$AIC
PO.fourcol3$criteria$AIC
PO.fivecol3$criteria$AIC
PO.sixcol3$criteria$AIC
PO.sevencol3$criteria$AIC
