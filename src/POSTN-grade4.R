### Andrew Trister
### Sage Bionetworks
### Seattle, WA
### 20120912

### Use the populated data to look at survival with POSTN in Grade 4 Gliomas (REMBRANDT)
##

require(ggplot2)
require(survival)
require(mixtools)
source("./src/ggkm.R")
##
## Now just the grade IV patients:

#and fit a model to it
hist(rembrandtEset['10631_mt',which(rembrandtPat$Grade==" IV")],main="Grade IV",xlab="POSTN Expression")
cutpoint <- quantile(rembrandtEset['10631_mt',which(rembrandtPat$Grade==" IV")])[3]+sd(rembrandtEset['10631_mt',which(rembrandtPat$Grade==" IV")])
cutpoint <- quantile(rembrandtEset['10631_mt',which(rembrandtPat$Grade==" IV")])[3]

POSTNclass.lowgrade$x <- POSTNclass.rembrandt$x[which(rembrandtPat$Grade==" IV" & rembrandtPat$Patient.ID!="E09448" )]
POSTNclass.lowgrade$posterior <- POSTNclass.rembrandt$posterior[which(rembrandtPat$Grade==" IV"& rembrandtPat$Patient.ID!="E09448" ),]


predictedClasses.lowgrade <- rep(1,length(POSTNclass.lowgrade$x))


#use this to repeat the previous calls
for(i in 1:length(POSTNclass.lowgrade$x)) {
  if(POSTNclass.lowgrade$posterior[i,1] > POSTNclass.lowgrade$posterior[i,2])
    predictedClasses.lowgrade[i] <- 2
}


#use this for cuts with cutpoint (one SD away from mean)
for(i in 1:length(POSTNclass.lowgrade$x)) {
  if(POSTNclass.lowgrade$x[i] < cutpoint)
    predictedClasses.lowgrade[i] <- 2
}


table(predictedClasses.lowgrade)

## let's make a survival object
lowgradePat <- c()
lowgradePat$surv <- gbmPat.rembrandt$surv[which(rembrandtPat$Grade==" IV"& rembrandtPat$Patient.ID!="E09448" )]
lowgradePat$survTime <- gbmPat.rembrandt$survTime[which(rembrandtPat$Grade==" IV"& rembrandtPat$Patient.ID!="E09448" )]
#lowgradePat$survTime <- 30.42* lowgradePat$survTime

tmpSurv.lowgrade <- Surv(lowgradePat$survTime,lowgradePat$surv)


##let's look at survival based on class for POSTN
ggkm(survfit(tmpSurv.lowgrade ~ predictedClasses.lowgrade), ystratalabs = (c("POSTN High", "POSTN Low")), 
     timeby = 12,
     main = "Grade IV K-M Plot By POSTN Expression (REMBRANDT)")