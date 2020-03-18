library(ostereotype)
library(cluster)
library(clustord)
library(tidyr)
library(tidyverse)

# process the data
ewave.ordinal <- read.csv("ewaveOrdinal.csv", header = TRUE)  # n=76, m=235
binary.ewave <- read.csv("ewaveBinary.csv", header = TRUE)    # n=76, m=235

ewave.data <- function(ewave.data.input){
  ewave.data.input <- gather(data = ewave.data.input, 
                             key = COL, 
                             value = Y, 
                             E001:E235,
                             factor_key = TRUE)
  names(ewave.data.input)[1] <- "ROW"
  ewave.data.input$COL <- str_remove(ewave.data.input$COL, "E")
  ewave.data.input$Y <- factor(ewave.data.input$Y)
  # remove missing data
  ewave.data.input <- ewave.data.input[complete.cases(ewave.data.input),]
  return(ewave.data.input)
}

ewave.ordinal <- ewave.data(ewave.ordinal)
binary.ewave <- ewave.data(binary.ewave)

### ewave.ordinal
## Row clustering, Y~row+row*column
init.OSM.2r <- c(mu2=0.4,mu3=0.6,mu4=0.8,u2=1.8,u3=2,alpha1=0.4,beta=rnorm(234),gamma=rnorm(234))  
OSM.2rows.clust <- rowclustering(formula = "Y~row*column", 
                                   model = "OSM", nclus.row = 2, 
                                   long.df = ewave.ordinal,
                                   initvect = init.OSM.2r)

