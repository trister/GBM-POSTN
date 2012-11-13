### Andrew Trister
### Sage Bionetworks
### Seattle, WA
### 20120912

### Use the populated data to look at survival with POSTN in Grade 2+3 Gliomas (REMBRANDT)
##

require(ggplot2)
require(survival)
require(mixtools)
source("./src/ggkm.R")

##
## Now just the grade II and III patients:

#and fit a model to it
hist(rembrandtEset['10631_mt',which(rembrandtPat$Disease!=" GBM" & rembrandtPat$Disease!=" MIXED")])

POSTNclass.lowgrade <- normalmixEM(rembrandtEset['10631_mt',
                                                 which(rembrandtPat$Disease!=" GBM" & rembrandtPat$Patient.ID!="E09448"
                                                       & rembrandtPat$Disease!=" MIXED")],
                                   lambda=0.5,mu=c(6.5,11.5),sigma=1)
plot(POSTNclass.lowgrade,density=T)
summary(POSTNclass.lowgrade)

predictedClasses.lowgrade <- rep(1,length(POSTNclass.lowgrade$x))

for(i in 1:length(POSTNclass.lowgrade$x)) {
  if(POSTNclass.lowgrade$posterior[i,1] > POSTNclass.lowgrade$posterior[i,2])
    predictedClasses.lowgrade[i] <- 2
}
table(predictedClasses.lowgrade)


#another option would be to use the predicted classes from the larger group, which might be better
#and just cull based on the disease

POSTNclass.lowgrade$x <- POSTNclass.rembrandt$x[which(rembrandtPat$Disease !=" GBM" & rembrandtPat$Patient.ID!="E09448"  & rembrandtPat$Disease != " MIXED")]
POSTNclass.lowgrade$posterior <- POSTNclass.rembrandt$posterior[which(rembrandtPat$Disease !=" GBM" & rembrandtPat$Patient.ID!="E09448" & rembrandtPat$Disease != " MIXED"),]



predictedClasses.lowgrade <- rep(1,length(POSTNclass.lowgrade$x))

for(i in 1:length(POSTNclass.lowgrade$x)) {
  if(POSTNclass.lowgrade$posterior[i,1] > POSTNclass.lowgrade$posterior[i,2])
    predictedClasses.lowgrade[i] <- 2
}
table(predictedClasses.lowgrade)

## let's make a survival object
lowgradePat <- c()
lowgradePat$surv <- gbmPat.rembrandt$surv[which(rembrandtPat$Disease!=" GBM" & rembrandtPat$Patient.ID!="E09448"  & rembrandtPat$Disease!=" MIXED" )]
lowgradePat$survTime <- gbmPat.rembrandt$survTime[which(rembrandtPat$Disease!=" GBM" & rembrandtPat$Patient.ID!="E09448"  & rembrandtPat$Disease!=" MIXED")]
#lowgradePat$survTime <- 30.42* lowgradePat$survTime

tmpSurv.lowgrade <- Surv(lowgradePat$survTime,lowgradePat$surv)


##let's look at survival based on class for POSTN
ggkm(survfit(tmpSurv.lowgrade ~ predictedClasses.lowgrade), ystratalabs = (c("POSTN High", "POSTN Low")), 
     timeby = 12,
     main = "Low Grade K-M Plot By POSTN Expression (REMBRANDT)")
