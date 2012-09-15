### Andrew Trister
### Sage Bionetworks
### Seattle, WA
### 20120912

### Use the populated data to look at survival with POSTN in all GBMS in TCGA and REMBRANDT


require(ggplot2)
require(survival)
require(mixtools)
source("./src/ggkm.R")

#let's now put together just the GBM from the REMBRANDT 

#hist(rembrandtEset['10631_mt',which(rembrandtPat$Grade==" IV")])
temp <- sub("_eg","_mt",rownames(eset))
overlapGenes <- which(rownames(rembrandtEset) %in% temp)
allGBMeset <- cbind(eset,rembrandtEset[overlapGenes,which(rembrandtPat$Disease==" GBM")])


hist(allGBMeset['10631_eg',])

#and fit a model to it
POSTNclassAll <- normalmixEM(allGBMeset['10631_eg',],lambda=0.5,mu=c(6.5,11.5),sigma=1)
plot(POSTNclassAll,density=T)
summary(POSTNclassAll)

predictedClassesAll <- rep(1,length(POSTNclassAll$x))

for(i in 1:length(POSTNclassAll$x)) {
  if(POSTNclassAll$posterior[i,1] > POSTNclassAll$posterior[i,2])
    predictedClassesAll[i] <- 2
}




## let's make a survival object

gbmPatAll <- c()
gbmPatAll$survTime <- c(gbmPat$survTime,30.42*gbmPat.rembrandt$survTime[which(rembrandtPat$Disease==" GBM")])
gbmPatAll$surv <- c(gbmPat$surv,gbmPat.rembrandt$surv[which(rembrandtPat$Disease==" GBM")])

tmpSurvAll <- Surv(gbmPatAll$survTime,gbmPatAll$surv)


##let's look at survival based on class for POSTN
ggkm(survfit(tmpSurvAll ~ predictedClassesAll), ystratalabs = (c("POSTN High", "POSTN Low")), 
     timeby = 365,
     main = "GBM K-M Plot By POSTN Expression (REMBRANDT+TCGA)")

