### Andrew Trister
### Sage Bionetworks
### Seattle, WA
### 20120912

### Use the populated data to look at survival with POSTN in REMBRANDT


require(ggplot2)
require(survival)
require(mixtools)
source("./src/ggkm.R")

hist(rembrandtEset['10631_mt',])

#and fit a model to it
POSTNclass.rembrandt <- normalmixEM(rembrandtEset['10631_mt',],lambda=0.5,mu=c(6.5,11.5),sigma=1)
plot(POSTNclass.rembrandt,density=T)
summary(POSTNclass.rembrandt)

predictedClasses.rembrandt <- rep(1,length(POSTNclass.rembrandt$x))

for(i in 1:length(POSTNclass.rembrandt$x)) {
  if(POSTNclass.rembrandt$posterior[i,1] > POSTNclass.rembrandt$posterior[i,2])
    predictedClasses.rembrandt[i] <- 2
}


## let's make a survival object
gbmPat.rembrandt <- c()
gbmPat.rembrandt$survTime <- as.numeric(rembrandtPat$Survival..months.)
gbmPat.rembrandt$surv <- ifelse(rembrandtPat$Survival..months.=="--", 0,1)
gbmPat.rembrandt$survTime[ which(gbmPat.rembrandt$surv==0)] <- as.numeric(rembrandtPat$Followup.Month)[ which(gbmPat.rembrandt$surv==0)]

tmpSurv.rembrandt <- Surv(gbmPat.rembrandt$survTime,gbmPat.rembrandt$surv)


##let's look at survival based on class for POSTN
ggkm(survfit(tmpSurv.rembrandt ~ predictedClasses.rembrandt), ystratalabs = (c("POSTN High", "POSTN Low")), 
     timeby = 12,
     main = "GBM K-M Plot By POSTN Expression (REMBRANDT)")


#boxplots:

plot(rembrandtEset['10631_mt',which(rembrandtPat$Grade!=" --" & rembrandtPat$Grade!=" I")]~factor(rembrandtPat$Grade[which(rembrandtPat$Grade!=" --" & rembrandtPat$Grade!=" I")]),
     xlab="Glioma Grade", ylab="Expression",main="Expression of POSTN")

aov.rembrandt <- aov(rembrandtEset['10631_mt',which(rembrandtPat$Grade!=" --" &rembrandtPat$Grade!=" I")]~factor(rembrandtPat$Grade[which(rembrandtPat$Grade!=" --" & rembrandtPat$Grade!=" I")]))
summary(aov.rembrandt)
print(model.tables(aov.rembrandt,"means"),digits=3)





