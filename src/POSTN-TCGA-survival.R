### Andrew Trister
### Sage Bionetworks
### Seattle, WA
### 20120912

### Use the populated data to look at survival with POSTN in TCGA


require(ggplot2)
require(survival)
require(mixtools)
source("./src/ggkm.R")



###
###
# Now let's look at POSTN (10631)
###
###

hist(eset['10631_eg',])

#and fit a model to it
POSTNclass <- normalmixEM(eset['10631_eg',],lambda=0.5,mu=c(6.5,11.5),sigma=1)
plot(POSTNclass,density=T)
summary(POSTNclass)

predictedClasses <- rep(1,length(POSTNclass$x))

for(i in 1:length(POSTNclass$x)) {
  if(POSTNclass$posterior[i,1] > POSTNclass$posterior[i,2])
    predictedClasses[i] <- 2
}


#use this for cuts with cutpoint (one SD away from mean)
#cutpoint <- quantile(eset['10631_eg',])[3]
#for(i in 1:length(POSTNclass$x)) {
#  if(POSTNclass$x[i] < cutpoint)
#    predictedClasses[i] <- 2
#}


## let's make a survival object
metadata$vital_status[metadata$vital_status=="[Not Available]"] <- NA
gbmPat <- c()
gbmPat$survTime <- as.numeric(metadata$days_to_death)
gbmPat$surv <- ifelse(metadata$vital_status=="DECEASED", 1,0)
gbmPat$survTime[ which(metadata$vital_status=="LIVING")] <- as.numeric(metadata$days_to_last_followup)[ which(metadata$vital_status=="LIVING")]

tmpSurv <- Surv(gbmPat$survTime,gbmPat$surv)


##let's look at survival based on class for POSTN
ggkm(survfit(tmpSurv ~ predictedClasses), ystratalabs = (c("POSTN High", "POSTN Low")), 
     timeby = 365,
     main = "GBM K-M Plot By POSTN Expression (TCGA)")

