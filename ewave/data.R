#############################################################################
# Read data from the three data sources
# Bequia - A variety of English
# ewave - Many varieties of English
# grambank - World languages
#############################################################################

#############################################################################
# List of files
datdir <- "../STAT487/reformatted/"
fnames <- list.files(datdir)
cbind(fnames)
##  [1,] "BQewavebmat.csv"      
##  [2,] "BQewavenmat.csv"      
##  [3,] "BQewaveymat.csv"      
##  [4,] "BQgrambankbmat.csv"   
##  [5,] "BQgrambanknmat.csv"   
##  [6,] "BQgrambankymat.csv"   
##  [7,] "BQspeakers.csv"       
##  [8,] "ewaveBinary.csv"      
##  [9,] "ewaveBinaryBQ.csv"    
## [10,] "ewaveFeatures.csv"    
## [11,] "ewaveFeaturesBQ.csv"  
## [12,] "ewaveGrambank.csv"    
## [13,] "ewaveOrdinal.csv"     
## [14,] "ewaveOrdinalBQ.csv"   
## [15,] "ewaveParameters.psv"  
## [16,] "ewaveParametersBQ.psv"
## [17,] "gbBinary.csv"         
## [18,] "gbBinaryEwave.csv"    
## [19,] "gbFeatures.csv"       
## [20,] "gbFeaturesEwave.csv"  
## [21,] "gbParameters.psv"     
## [22,] "gbParametersEwave.psv"

hh <- lapply(fnames,function(ff) {
  print(ff)
  rnames <- 1
  if(ff%in%c("BQspeakers.csv",
             "ewaveParameters.psv",
             "ewaveParametersBQ.psv",
             "gbFeatures.csv",
             "gbFeaturesEwave.csv",
             "gbParameters.psv",
             "gbParametersEwave.psv")) rnames <- NULL
  sep <- ","
  if(grepl("psv$",ff)) sep <- "|"
  datf <- read.table(paste0(datdir,ff),quote="\"",sep=sep,
                     header=TRUE,stringsAsFactors=FALSE,row.names=rnames,fill=TRUE)
  objname <- strsplit(ff,".",fixed=TRUE)[[1]][1]
  cmd <- paste0(objname," <<- datf")
  eval(parse(text=cmd))
  return(list(ff, names(datf), dim(datf), objname))
})
names(hh) <- fnames
t(sapply(hh, function(x) x[[3]]))
t(sapply(hh, function(x) x[[4]]))

# Bequia data: n=19 speakers, 36 ewave features (efeatures),
#                              4 grambank features (gfeatures)
## BQewavebmat.csv           19   36  - binary 0/1 speakers x efeatures
## BQewavenmat.csv           19   36  - number of opportunities: speakers x efeatures
## BQewaveymat.csv           19   36  - number of usages: speakers x efeatures         
## BQgrambankbmat.csv        19    4  - binary 0/1 speakers x gfeatures                
## BQgrambanknmat.csv        19    4  - number of opportunities: speakers x gfeatures  
## BQgrambankymat.csv        19    4  - number of usages: speakers x gfeatures         
## BQspeakers.csv            19    7  - speaker characteristics

# Ewave data: n=76 varieties of English, and 235 ewave features (efeatures)
#                                        and  26 grambank features (gfeatures)
## ewaveFeatures.csv         76  235  - full coded features ABCD: languages x efeatures
## ewaveFeaturesBQ.csv       76   36  - restriction to 36 efeatures coded in BQ
## ewaveBinary.csv           76  235  - binary 0/1: languages x efeatures
## ewaveBinaryBQ.csv         76   36  - restriction to 36 efeatures coded in BQ
## ewaveOrdinal.csv          76  235  - recoded features 1,2,3,4: languages x efeatures
## ewaveOrdinalBQ.csv        76   36  - restriction to 36 efeatures coded in BQ
## ewaveGrambank.csv         76   26  - binary 0/1: languages x gfeatures
## ewaveParameters.psv      235    8  - list of ewave feature descriptions
## ewaveParametersBQ.psv     36    8  - restriction to the 36 efeatures coded in BQ

# Grambank data: n=872 languages, 195 grambank features (gfeatures)
## gbBinary.csv             872  195  - binary 0/1 languages x gfeatures
## gbBinaryEwave.csv        872   26  - restriction to 26 gfeatures coded in ewave
## gbFeatures.csv        117739   10  - token values for the Grambank features
## gbFeaturesEwave.csv    17406   10  - restriction to the 26 gfeatures coded in ewave
## gbParameters.psv         195   32  - list of GB feature descriptions
## gbParametersEwave.psv     26   32  - restriction to 26 gfeatures coded in ewave

#############################################################################


#------settings------
n <- 100            #sample size                          
mu <- 5             #this is unknown in practice                         
beta <- 2.7         #this is unknown in practice
sigma <- 0.15       #this is unknown in practice
#--------------------

#------set the seed so that this example can be replicated------
set.seed(937)
#---------------------------------------------------------------

#------generate 1000 data sets and store betaHat------
betaHat <- numeric(1000)
for(i in 1:1000){
  #generate the binary covariate --> n Bernoulli trials
  x <- sample(x=c(0, 1), size=n, replace=TRUE, prob=c(0.5, 0.5))
  #generate the errors
  epsilon <- rnorm(n=n, mean=0, sd=sigma)
  #form the response variable      
  y <- mu + beta * x + epsilon 
  #the ith generated data set
  data_i <- data.frame(y=y, x=x)
  #fit the model
  mod <- lm(y~x, data=data_i)
  #store the estimate of beta
  betaHat[i] <- as.numeric(coef(mod)[2])     
  }    
#-----------------------------------------------------
#------E(betaHat) = beta?------
mean(betaHat)
# [1] 2.698609
#---------------------------
