#######################################################################
## FUNCTION TO GENERATE A SPACED MOSAIC PLOT FOR ROW CLUSTERING #######
#######################################################################
#Author: Daniel Fernandez. Victoria University of Wellington. Sep2013##
#######################################################################

#phi.est: vector. Estimated "score" parameters (phi) regarding the Ordinal Stereotype Model (according to the increasing constraint defined in Anderson, 1984)
#q: integer. Number of ordinal categories.
#n: integer. Number of rows (e.g. spp.)
#m: integer. Number of columns (e.g. sites)
#y.mat: matrix(n x m). Ordinal data.
#data.label: list. first dimension contains the categories labels (categ),
#                  second dimension contains the cluster labels (cluster),
#                  third dimension contains the row labels in the data (row),
#                  four dimension contains the columns labels in the data (col)
#ClusterRowY: vector. Contain the row cluster for each row in the data.


#Author: Daniel Fernandez. Victoria University of Wellington. Sep2013

spaced.mosaic.plot <- function(y.mat, v.phi.est, R.best, ClusterRowY, labels){
  n <- nrow(y.mat)
  m <- ncol(y.mat)
  q <- length(unique(c(y.mat))) #number of categories data

  ClusterRowY.sorted <- sort(ClusterRowY)  
  
  phi.space <- v.phi.est*5

  freq.matrix <- matrix(0,nrow=q,ncol=n)
  rownames(freq.matrix) <- seq(1,q,1)
  for (i in 1:n) freq.matrix[names(table(y.mat[i,])),i] <-  table(y.mat[i,])
  rownames(freq.matrix) <- labels$categ
  colnames(freq.matrix) <- labels$row
  
  freq.matrix <- freq.matrix[,rownames(ClusterRowY.sorted)]
  colnames(freq.matrix) <- unname(ClusterRowY.sorted)

  retmat <- matrix(NA,nrow=q,ncol=R.best)
  retmat <- rowsum(t(freq.matrix), group=unname(ClusterRowY.sorted))
  qname <- colnames(retmat)
  cname <- c()
  for (r in 1:R.best) 
    cname <- c(cname,paste("R",as.character(r),sep=""))
  retmat <- unname(retmat) 
  dimnames(retmat) <- list(Row=cname, Col=qname)

  phi.space.increm <- array(NA,q-1)
  for (l in 1:(q-1)) phi.space.increm[l] <- phi.space[l+1]-phi.space[l]
  R.space.increm <- rep(0.5,R.best-1)
  my.spacing <- list(unit(R.space.increm, "lines"),
                   unit(phi.space.increm, "lines"))

  retmat <- as.table(retmat)
  
  pdf(paste("MosaicPlot_R=", R.best,".pdf",sep=""), paper="a4r") 
  mosaic(retmat, main = "Row Clustering Results",
       tl_labels=c(TRUE, TRUE),        
       set_varnames = c(Row = "Row Cluster",Col = "Ordinal Categories"),
       gp_varnames = gpar(fontsize = 12, fontface = 2),
       gp_labels = gpar(fontsize = 10, fontface = 2),
       gp = gpar(fill = "lightsalmon"),
       #gp_text = gpar(text=retmat),
       margins = c(left = 0, right = 0, top = 5, bottom = 3),
       rot_labels = c(0, 0, 0, 0), pop=FALSE)
       labeling_cells(text = retmat, clip = FALSE, gp_text = gpar(fontsize = 12, fontface = 3))(retmat) 

  dev.off()

  pdf(paste("MosaicPlot_SPACING_R=", R.best,".pdf",sep=""),paper="a4r") 
  mosaic(retmat, main = "Row Clustering Results. Scaled Space (Fitted Scores)",
       tl_labels = c(TRUE, TRUE),
       set_varnames = c(Row = "Row Cluster",Col = "Ordinal Categories"),
       gp_varnames = gpar(fontsize = 12, fontface = 2),
       gp_labels = gpar(fontsize = 10, fontface = 2),
       gp = gpar(fill = "lightsalmon"),
       margins = c(left = 0, right = 0, top = 5, bottom = 3),
       rot_labels = c(0, 0, 0, 0),
       spacing = my.spacing, pop=FALSE)	

       color.fill <- c("yellow","red","green","orange","black","grey","mediumorchid1","blue")
       for (r in 1:R.best)
       {
  	 for (k in 1:q)
  	 { 
           seekViewport(paste("cell:Row=R",r,",Col=",labels$categ[k],sep=""))
    	   grid.rect(x=1, y=0.3,
           width=unit(phi.space.increm[k]-0.02,"lines"),height=unit(0.25,"lines"), 
	   just="left",gp=gpar(fill=color.fill[k], lty=1))
  	}
       }

       labeling_cells(text = retmat, clip = FALSE, gp_text = gpar(fontsize = 12, fontface = 3))(retmat)
    
  dev.off()

  sum.of.the.rows <- rowSums(freq.matrix, na.rm = FALSE, dims = 1)
  without.row.cluster <- as.table(t(matrix(sum.of.the.rows,nrow=q,ncol=1)))
  without.row.cluster <- unname(without.row.cluster) 
  dimnames(without.row.cluster) <- list(Row=" ", Col=qname)
  
  pdf(paste("MosaicPlot_withoutClustering.pdf",sep=""),paper="a4r")
  mosaic(without.row.cluster, main = "Results without Row Clustering/Spacing",
       tl_labels = c(TRUE, TRUE),
       set_varnames = c(Row = "Without Row Clustering",Col = "Ordinal Categories"),
       gp_varnames = gpar(fontsize = 12, fontface = 2),
       gp_labels = gpar(fontsize = 10, fontface = 2),
       gp = gpar(fill = "lightsalmon"),
       margins = c(left = 0, right = 0, top = 5, bottom = 3),
       rot_labels = c(0, 0, 0, 0), pop=FALSE)
       labeling_cells(text = without.row.cluster, clip = FALSE,
                      gp_text = gpar(fontsize = 12, fontface = 3))(without.row.cluster)

  dev.off()
  
  return(retmat)
}


#######################################################################
## FUNCTION TO GENERATE A SPACED MOSAIC PLOT FOR COLUMN CLUSTERING #######
#######################################################################

spaced.mosaic.plot2 <- function(y.mat, v.phi.est, C.best, ClusterColY, labels){
  n <- nrow(y.mat)
  m <- ncol(y.mat)
  q <- length(unique(c(y.mat))) #number of categories data
  
  ClusterColY.sorted <- sort(ClusterColY)  
  
  phi.space <- v.phi.est*5
  
  freq.matrix <- matrix(0, nrow=q, ncol=m)
  rownames(freq.matrix) <- seq(1,q,1)
  for (i in 1:m) 
    freq.matrix[names(table(y.mat[,i])),i] <-  table(y.mat[,i])
  rownames(freq.matrix) <- labels$categ
  colnames(freq.matrix) <- labels$col
  
  freq.matrix <- freq.matrix[,rownames(ClusterColY.sorted)]
  colnames(freq.matrix) <- unname(ClusterColY.sorted)
  
  retmat <- matrix(NA,nrow=q, ncol=C.best)
  retmat <- rowsum(t(freq.matrix), group=unname(ClusterColY.sorted))
  qname <- colnames(retmat)
  cname <- c()
  for (c in 1:C.best) 
    cname <- c(cname,paste("C", as.character(c), sep=""))
  retmat <- unname(retmat) 
  dimnames(retmat) <- list(Row=cname, Col=qname)
  
  phi.space.increm <- array(NA, q-1)
  for (l in 1:(q-1)) 
    phi.space.increm[l] <- phi.space[l+1]-phi.space[l]
  C.space.increm <- rep(0.5,C.best-1)
  my.spacing <- list(unit(C.space.increm, "lines"),
                     unit(phi.space.increm, "lines"))
  
  retmat <- as.table(retmat)
  
  pdf(paste("MosaicPlot_C=", C.best,".pdf",sep=""), paper="a4r") 
  mosaic(retmat, main = "Column Clustering Results",
         tl_labels=c(TRUE, TRUE),        
         set_varnames = c(Row = "Column Cluster", Col = "Ordinal Categories"),
         gp_varnames = gpar(fontsize = 12, fontface = 2),
         gp_labels = gpar(fontsize = 10, fontface = 2),
         gp = gpar(fill = "lightsalmon"),
         margins = c(left = 0, right = 0, top = 5, bottom = 3),
         rot_labels = c(0, 0, 0, 0), pop=FALSE)
  labeling_cells(text = retmat, clip = FALSE, gp_text = gpar(fontsize = 12, fontface = 3))(retmat) 
  
  dev.off()
  
  pdf(paste("MosaicPlot_SPACING_C=", C.best,".pdf", sep=""), paper="a4r") 
  mosaic(retmat, main = "Column Clustering Results. Scaled Space (Fitted Scores)",
         tl_labels = c(TRUE, TRUE),
         set_varnames = c(Row = "Column Cluster", Col = "Ordinal Categories"),
         gp_varnames = gpar(fontsize = 12, fontface = 2),
         gp_labels = gpar(fontsize = 10, fontface = 2),
         gp = gpar(fill = "lightsalmon"),
         margins = c(left = 0, right = 0, top = 5, bottom = 3),
         rot_labels = c(0, 0, 0, 0),
         spacing = my.spacing, pop=FALSE)	
  
  color.fill <- c("yellow","red","green","orange","black","grey","mediumorchid1","blue")
  
  for (c in 1:C.best){
    for (k in 1:q) 
    {
      seekViewport(paste("cell:Row=C", c, ",Col=", labels$categ[k], sep=""))
      grid.rect(x=1, y=0.3,
                width=unit(phi.space.increm[k]-0.02,"lines"), height=unit(0.25,"lines"), 
                just="left", gp=gpar(fill=color.fill[k], lty=1))
    }
  }
  
  labeling_cells(text = retmat, clip = FALSE, gp_text = gpar(fontsize = 12, fontface = 3))(retmat)
  
  dev.off()
  
  sum.of.the.columns <- rowSums(freq.matrix, na.rm = FALSE, dims = 1)
  without.column.cluster <- as.table(t(matrix(sum.of.the.columns, nrow=q, ncol=1)))
  without.column.cluster <- unname(without.column.cluster) 
  dimnames(without.column.cluster) <- list(Row=" ", Col=qname)
  
  pdf(paste("MosaicPlot_withoutCluster.pdf",sep=""), paper="a4r")
  mosaic(without.column.cluster, main = "Results without Clustering/Spacing",
         tl_labels = c(TRUE, TRUE),
         set_varnames = c(Row = "Without Clustering", Col = "Ordinal Categories"),
         gp_varnames = gpar(fontsize = 12, fontface = 2),
         gp_labels = gpar(fontsize = 10, fontface = 2),
         gp = gpar(fill = "lightsalmon"),
         margins = c(left = 0, right = 0, top = 5, bottom = 3),
         rot_labels = c(0, 0, 0, 0), pop=FALSE)
  labeling_cells(text = without.column.cluster, clip = FALSE,
                 gp_text = gpar(fontsize = 12, fontface = 3))(without.column.cluster)
  
  dev.off()
  
  return(retmat)
}

