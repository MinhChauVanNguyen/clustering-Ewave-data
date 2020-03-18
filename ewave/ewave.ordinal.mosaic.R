

############################ SPACED MOSAIC PLOT ################################

library(vcd)
library(grid)
setwd("~/Documents/STAT487/reformatted")

ewave <- read.csv("ewaveOrdinal.csv", header = TRUE)  # n=76, m=235
ewave.ordinal.mosaic <- ewave[complete.cases(ewave),]

ewave.mat <- as.matrix(ewave.ordinal.mosaic)
head(ewave.mat)
ewave.mat <- ewave.mat[,-1]
n <- nrow(ewave.mat)
m <- ncol(ewave.mat)
n;m
phi.vect <- c(0, 0.4, 0.6, 1)
R <- 5

Labels <- list(categ = c("1","2","3","4"),
               cluster = paste("R", 1:R, sep=""),
               row = paste("r", 1:n, sep=""),
               col = paste("c", 1:m, sep=""))

ClusterRowY <- array(NA,n)
for(i in 1:n){
  ClusterRowY[i] <- sample(1:R, 1, replace=TRUE)
}
rownames(ClusterRowY) <- Labels$row
ClusterRowY


spaced.mosaic.plot(y.mat = ewave.mat,
               v.phi.est = phi.vect,
                  R.best = R,
             ClusterRowY = ClusterRowY,
                  labels = Labels)



C <- 7
Labels.COL <- list(categ = c("1","2","3","4"),
               cluster = paste("C", 1:C, sep=""),
               row = paste("r", 1:n, sep=""),
               col = paste("c", 1:m, sep=""))

ClusterColY <- array(NA,m)
for(i in 1:m){
  ClusterColY[i] <- sample(1:C, 1, replace=TRUE)
}
rownames(ClusterColY) <- Labels.COL$col
ClusterColY
spaced.mosaic.plot2(y.mat = ewave.mat,
                   v.phi.est = phi.vect,
                   C.best = C,
                   ClusterColY = ClusterColY,
                   labels = Labels.COL)




############################  ELBOW PLOT FUNCTION  ###################################

elbow.func <- function(aic2, aic3, aic4, aic5, aic6, aic7, type){
  if (!is.character(type) || !is.vector(type) || length(type)!=1) 
    stop("type must be a string, 'ROW' or 'COL'.")
  if 
  (!(type %in% c("ROW","COL"))) 
    stop("type must be either 'ROW' or COL'.")
  if(type=="ROW"){
    x <- c(2,3,4,5)
    y <- c(aic2=aic2,aic3=aic3,aic4=aic4,aic5=aic5)
    elbow <- plot(x, y, xlim = range(x), ylim = range(y),
                  xlab = "Number of Row Clusters",
                  ylab = "AIC score",
                  main = "Elbow Plot for determining the \n number of Row clusters",
                  pch = 8, cex = 2, lwd = 2)
    lines(x[order(x)], y[order(x)], xlim=range(x), ylim=range(y), pch=16)
    symbols(x = 5, y=aic5, circles=0.08, add=T, inches=F, fg="red")
  }else{
    x <- c(2,3,4,5,6,7)
    y <- c(aic2=aic2,aic3=aic3,aic4=aic4,aic5=aic5,aic6=aic6,aic7=aic7)
    elbow <- plot(x, y, xlim = range(x), ylim = range(y),
       xlab = "Number of Column Clusters",
       ylab = "AIC score",
       main = "Elbow Plot for determining the \n number of Column clusters",
       pch = 8, cex = 2, lwd = 2)
    lines(x[order(x)], y[order(x)], xlim=range(x), ylim=range(y), pch=16)
    symbols(x = 7, y=aic7, circles=0.2, add=T, inches=F, fg="blue")
  }
  print(elbow)
}

elbow.func(Bin2.row,Bin3.row,Bin4.row,Bin5.row, type="ROW")
elbow.func(Bin2.col,Bin3.col,Bin4.col,Bin5.col, Bin6.col, Bin7.col, type="COL")


############################# ELBOW PLOT FOR OSM #####################################

two <- OSM.2rows.cluster2$criteria$AIC 
three <- OSM.3rows.cluster2$criteria$AIC 
four <- OSM.4rows.cluster2$criteria$AIC 
five <- OSM.5rows.cluster2$criteria$AIC 

two.col <- OSM2.2col$criteria$AIC
three.col <- OSM2.3col$criteria$AIC
four.col <- OSM2.4col$criteria$AIC
five.col <- OSM2.5col$criteria$AIC
six.col <- OSM2.6col$criteria$AIC
seven.col <- OSM2.7col$criteria$AIC

############################ ELBOW PLOT FOR POM ################################

PO2.row <- PO.tworow$criteria$AIC
PO3.row <- PO.threerow$criteria$AIC
PO4.row <- PO.fourrow$criteria$AIC
PO5.row <- PO.fiverow$criteria$AIC

PO2.col <- PO.twocol$criteria$AIC
PO3.col <- PO.threecol$criteria$AIC
PO4.col <- PO.fourcol$criteria$AIC
PO5.col <- PO.fivecol$criteria$AIC
PO6.col <- PO.sixcol$criteria$AIC
PO7.col <- PO.sevencol$criteria$AIC


############################ ELBOW PLOT FOR BINARY ################################

Bin2.row <- tworow.binary$criteria$AIC
Bin3.row <- threerow.binary$criteria$AIC
Bin4.row <- fourrow.binary$criteria$AIC
Bin5.row <- fiverow.binary$criteria$AIC

Bin2.col <- twocol.binary$criteria$AIC
Bin3.col <- threecol.binary$criteria$AIC
Bin4.col <- fourcol.binary$criteria$AIC
Bin5.col <- fivecol.binary$criteria$AIC
Bin6.col <- sixcol.binary$criteria$AIC
Bin7.col <- sevencol.binary$criteria$AIC

