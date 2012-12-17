### Andrew Trister
### Sage Bionetworks
### Seattle, WA
### 20121112

### Build an elasticnet model trained on POSTN class in TCGA to learn what other genes may be implicated
### and verify in REMBRANDT


require(synapseClient)
require(glmnet)
require(randomForest)
require(caret)
require(affy)
require(survival)
require(pROC)
require(pls)
require(gplots)
require(mixtools)
source("./src/ggkm.R")
library(org.Hs.eg.db)






#Now let's get an object to have geneIDs from Entrez
x <- org.Hs.egSYMBOL
# Get the gene symbol that are mapped to an entrez gene identifiers
mapped_genes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_genes])
if(length(xx) > 0) {
  # Get the SYMBOL for the first five genes
  xx[1:5]
  # Get the first one
  xx[[1]]
}




#fit a model to the distribution of POSTN in TCGA
POSTNclass <- normalmixEM(eset['10631_eg',],lambda=0.5,mu=c(6.5,11.5),sigma=1)
predictedClasses <- rep(1,length(POSTNclass$x))

for(i in 1:length(POSTNclass$x)) {
  if(POSTNclass$posterior[i,1] > POSTNclass$posterior[i,2])
    predictedClasses[i] <- 2
}


POSTNclass.rembrandt <- normalmixEM(rembrandtEset['10631_mt',],lambda=0.5,mu=c(6.5,11.5),sigma=1)
predictedClasses.rembrandt <- rep(1,length(POSTNclass.rembrandt$x))

for(i in 1:length(POSTNclass.rembrandt$x)) {
  if(POSTNclass.rembrandt$posterior[i,1] > POSTNclass.rembrandt$posterior[i,2])
    predictedClasses.rembrandt[i] <- 2
}


#cull the REMBRANDT U133Aplus2 to U133A genesets
temp <- unlist(lapply(rownames(rembrandtEset),function(x){
  return(gsub("_mt","_eg",x))
}))

rembrandtEset.culled <- rembrandtEset
rownames(rembrandtEset.culled) <- temp

temp <- intersect(temp,rownames(eset))

eset <- eset[temp,]
rembrandtEset.culled <- rembrandtEset.culled[temp,]

eset.culled <- eset[-grep("10631_eg",rownames(eset)),]
rembrandtEset.culled <- rembrandtEset.culled[-grep("10631_eg",rownames(rembrandtEset.culled)),]


# run the elastic net model
#  cv.fit <- cv.glmnet(x=t(eset.culled), y=eset["10631_eg",], nfolds=10, alpha=.1, family="gaussian")
#  plot(cv.fit)
#  fitEnet <- glmnet(x=t(eset.culled), y=eset["10631_eg",], family="gaussian", alpha=.1, lambda=cv.fit$lambda.min)
# 
# cv.fit <- cv.glmnet(x=t(eset), y=factor(predictedClasses), nfolds=10, alpha=.1, family="binomial")
# plot(cv.fit)
# fitEnet <- glmnet(x=t(eset), y=factor(predictedClasses), family="binomial", alpha=.1, lambda=cv.fit$lambda.min)


cv.fit <- cv.glmnet(x=t(eset.culled), y=factor(predictedClasses), nfolds=10, alpha=.1, family="binomial")
plot(cv.fit)
fitEnet <- glmnet(x=t(eset.culled), y=factor(predictedClasses), family="binomial", alpha=.1, lambda=cv.fit$lambda.min)














###############################################################################
#  Look now at all patients in REMBRANDT                                      #
###############################################################################


yhatEnet <- predict(fitEnet, t(rembrandtEset.culled), type="response", s="lambda.min")  #all of REMBRANDT

boxplot(yhatEnet ~ predictedClasses.rembrandt, ylab="3-year OS prediction (%)", xlab="predicted POSTN group", main="elastic net validation")
stripchart(yhatEnet ~ predictedClasses.rembrandt,pch=20, col="royalblue", vertical=TRUE, add=TRUE, cex=.6)

rocEnet <- roc(predictor=as.numeric(yhatEnet), response=as.numeric(predictedClasses.rembrandt-1),ci=TRUE)
plot.roc(rocEnet,col="red")

riskEnet <- as.vector(ifelse(yhatEnet >= median(yhatEnet), 1, 0))
names(riskEnet) <- rownames(yhatEnet)


plot(survfit(tmpSurv.rembrandt  ~ riskEnet), main="elastic net model", xlab="months",ylab="probability of OS",col= c("blue","magenta"),lwd=3)
survdiff(tmpSurv.rembrandt  ~ riskEnet, rho=0)



###############################################################################
#  Look now at only the Grade III                                             #
###############################################################################


lowgradePat <- c()
lowgradePat$surv <- gbmPat.rembrandt$surv[which(rembrandtPat$Grade==" III")]
lowgradePat$survTime <- gbmPat.rembrandt$survTime[which(rembrandtPat$Grade==" III")]


tmpSurv.lowgrade <- Surv(lowgradePat$survTime,lowgradePat$surv)


yhatEnet <- predict(fitEnet, t(rembrandtEset.culled[,which(rembrandtPat$Grade==" III")]), type="response", s="lambda.min")

boxplot(yhatEnet ~ predictedClasses.rembrandt[which(rembrandtPat$Grade==" III")], ylab="3-year OS prediction (%)", xlab="predicted POSTN group", main="elastic net validation")
stripchart(yhatEnet ~ predictedClasses.rembrandt[which(rembrandtPat$Grade==" III")],pch=20, col="royalblue", vertical=TRUE, add=TRUE, cex=.6)

rocEnet <- roc(predictor=as.numeric(yhatEnet), response=as.numeric(predictedClasses.rembrandt[which(rembrandtPat$Grade==" III")]),ci=TRUE)
plot.roc(rocEnet,col="red")

riskEnet <- as.vector(ifelse(yhatEnet >= median(yhatEnet), 1, 0))
names(riskEnet) <- rownames(yhatEnet)



ggkm(survfit(tmpSurv.lowgrade  ~ riskEnet),timeby=12,main="KM of REMBRANDT")
summary(survfit(tmpSurv.lowgrade  ~ riskEnet))

plot(survfit(tmpSurv.lowgrade ~ riskEnet), main="elastic net model", xlab="months",ylab="probability of OS",col= c("blue","magenta"),lwd=3)
survdiff(tmpSurv.lowgrade ~ riskEnet, rho=0)










###############################################################################
#  Look now at only the Grade II                                              #
###############################################################################


lowgradePat <- c()
lowgradePat$surv <- gbmPat.rembrandt$surv[which(rembrandtPat$Grade==" II")]
lowgradePat$survTime <- gbmPat.rembrandt$survTime[which(rembrandtPat$Grade==" II")]

tmpSurv.lowgrade <- Surv(lowgradePat$survTime,lowgradePat$surv)


yhatEnet <- predict(fitEnet, t(rembrandtEset.culled[,which(rembrandtPat$Grade==" II")]), type="response", s="lambda.min") 

boxplot(yhatEnet ~ predictedClasses.rembrandt[which(rembrandtPat$Grade==" II")], ylab="3-year OS prediction (%)", xlab="predicted POSTN group", main="elastic net validation")
stripchart(yhatEnet ~ predictedClasses.rembrandt[which(rembrandtPat$Grade==" II")],pch=20, col="royalblue", vertical=TRUE, add=TRUE, cex=.6)

rocEnet <- roc(predictor=as.numeric(yhatEnet), response=as.numeric(predictedClasses.rembrandt[which(rembrandtPat$Grade==" II")]),ci=TRUE)
plot.roc(rocEnet,col="red")

riskEnet <- as.vector(ifelse(yhatEnet >= 0.15, 1, 0))
names(riskEnet) <- rownames(yhatEnet)



ggkm(survfit(tmpSurv.lowgrade  ~ riskEnet),timeby=12,main="KM of REMBRANDT")
summary(survfit(tmpSurv.lowgrade  ~ riskEnet))


plot(survfit(tmpSurv.lowgrade ~ riskEnet), main="elastic net model", xlab="months",ylab="probability of OS",col= c("blue","magenta"),lwd=3)
survdiff(tmpSurv.lowgrade ~ riskEnet, rho=0)











heatmap(x=rembrandtEset.culled[which(abs(fitEnet$beta)>0),],
          ColSideColors=c("red","green","blue","yellow")[factor(rembrandtPat$Disease)],
          scale="row",col=redgreen(50),main="Elastic Net")



gene.names <- lapply(rownames(eset)[which(abs(fitEnet$beta)>0)],function(x){
  return(xx[strsplit(x,"_eg")[[1]]])
})

paste(unlist(gene.names),collapse=" ")



##
# Look at the Philips data

phillipsEset.culled <- phillipsEset[-grep("10631_mt",rownames(phillipsEset)),recurrent]
phillipsEset.culled <- cbind(phillipsEset.culled,phillipsEset[-grep("10631_mt",rownames(phillipsEset)),matchedPrimary])

yhatEnet <- predict(fitEnet, t(phillipsEset.culled), type="response", s="lambda.min")
boxplot(yhatEnet ~ c(rep(1,23),rep(0,23)), ylab="3-year OS prediction (%)", xlab="predicted POSTN group", main="elastic net validation")
stripchart(yhatEnet ~ c(rep(1,23),rep(0,23)),pch=20, col="royalblue", vertical=TRUE, add=TRUE, cex=.6)



riskEnet <- as.vector(ifelse(yhatEnet >= median(yhatEnet), 1, 0))
names(riskEnet) <- rownames(yhatEnet)


plot(survfit(tmpSurv.phillips  ~ riskEnet[1:23]), main="elastic net model", xlab="months",ylab="probability of OS",col= c("blue","magenta"),lwd=3)
survdiff(tmpSurv.rembrandt  ~ riskEnet, rho=0)






### Build a more robust classifer using bootstrapping



# first rescale the data (chemores) to have the same mean and variance than the LUAD
# Justin's function to rescale the VS to get the same mean/var than the TS
normalize_to_X <- function(mean.x, sd.x, Y){
  m.y <- rowMeans(Y)
  sd.y <- apply(Y, 1, sd)
  Y.adj <- (Y - m.y) * sd.x / sd.y  + mean.x 
  Y.adj[sd.y == 0] <- mean.x[sd.y==0]
  Y.adj
}

tmp <- apply(eset.culled,1,sd)
rembrandtEset.scaled <- normalize_to_X(rowMeans(eset.culled),tmp,rembrandtEset.culled)

#get rid of the least variant probes
tmp1 <- which(tmp>quantile(tmp,probs=0.2))
eset.culled2 <- eset.culled[tmp1,]
rembrandtEset.scaled <- rembrandtEset.scaled[tmp1,]
rm(tmp,tmp1)


require(glmnet)


# alpha = 0.1 optimize lambda nfolds=10
# bootstrap 100 times
N <- 100
fit <- c()
selected <- rep(0,length(rownames(eset.culled2)))
features <- c()
yhat_REMBRANDT <- c()
models <- 0
i <- 0
for(i in 1:N) {
  
  j <- sample(which(predictedClasses>0),replace=TRUE)
  cv.fit <- cv.glmnet(x=t(eset.culled2[,j]), y=factor(predictedClasses[j]), nfolds=10, alpha=.1, family="binomial")
  fit <- glmnet(x=t(eset.culled2[,j]),y=factor(predictedClasses[j]),family="binomial",alpha=.1,lambda=cv.fit$lambda.1se)
  features <- c(features,fit$beta)
  print(i)
}

  
#   if(length(which(abs(as.numeric(fit$beta))> 10^-4))>10) {
#     i=i+1
#     print(i)
#     features <- c(features,length(which(abs(as.numeric(fit$beta))> 10^-5)))
#     selected <- selected+(abs(as.numeric(fit$beta))>10^-5+0)
#     yhat_REMBRANDT <- cbind(yhat_REMBRANDT,predict(fit, t(rembrandtEset.scaled),type="response"))
#     models <- dim(yhat_REMBRANDT)[2]
#   } 
# }

 
# gene.names <- lapply(rownames(eset.culled2)[which(selected>90)],function(x){
#   return(xx[strsplit(x,"_eg")[[1]]])
# })
# 
# paste(unlist(gene.names),collapse=" ")


# weighted aggregation
#function from In Sock to help with the determination of the best features to select
weightAggregation<-function(resultsModel){
  
  ResultBS<-c()
  for(k in 1:length(resultsModel)){  
    #     a<-abs(resultsModel[[k]][-1])
    #     A<-sort(a,decreasing = T,index.return=T)
    ResultBS<-cbind(ResultBS,rank(abs(as.numeric(resultsModel[[k]][-1])),ties.method="min")/length(resultsModel[[k]]))
  }
  rownames(ResultBS)<-rownames(resultsModel[[1]])[-1]  
  reference <- apply(ResultBS,1,sum) 
  return(reference)
}


selected <- weightAggregation(features)
plot(sort(selected))

#use the weighted aggregation to grab the top 100
select_features <- sort(selected,decreasing=TRUE)[1:100]


heatmap(x=rembrandtEset.scaled[names(sort(select_features)[90:100]),],
        ColSideColors=c("red","green","blue","yellow")[factor(rembrandtPat$Disease)],
        scale="none",col=redgreen(50),main="Elastic Net")

gene.names <- lapply(names(sort(select_features)[90:100]),function(x){
  return(xx[strsplit(x,"_eg")[[1]]])
})

paste(unlist(gene.names),collapse=" ")


###############################################################################
#  Look now at only the Grade III                                             #
###############################################################################


lowgradePat <- c()
lowgradePat$surv <- gbmPat.rembrandt$surv[which(rembrandtPat$Grade==" III")]
lowgradePat$survTime <- gbmPat.rembrandt$survTime[which(rembrandtPat$Grade==" III")]


tmpSurv.lowgrade <- Surv(lowgradePat$survTime,lowgradePat$surv)


yhatEnet <- predict(fitEnet, t(rembrandtEset.culled[,which(rembrandtPat$Grade==" III")]), type="response", s="lambda.min")

boxplot(yhatEnet ~ predictedClasses.rembrandt[which(rembrandtPat$Grade==" III")], ylab="3-year OS prediction (%)", xlab="predicted POSTN group", main="elastic net validation")
stripchart(yhatEnet ~ predictedClasses.rembrandt[which(rembrandtPat$Grade==" III")],pch=20, col="royalblue", vertical=TRUE, add=TRUE, cex=.6)

rocEnet <- roc(predictor=as.numeric(yhatEnet), response=as.numeric(predictedClasses.rembrandt[which(rembrandtPat$Grade==" III")]),ci=TRUE)
plot.roc(rocEnet,col="red")

riskEnet <- as.vector(ifelse(yhatEnet >= median(yhatEnet), 1, 0))
names(riskEnet) <- rownames(yhatEnet)



ggkm(survfit(tmpSurv.lowgrade  ~ riskEnet),timeby=12,main="KM of REMBRANDT")
summary(survfit(tmpSurv.lowgrade  ~ riskEnet))

plot(survfit(tmpSurv.lowgrade ~ riskEnet), main="elastic net model", xlab="months",ylab="probability of OS",col= c("blue","magenta"),lwd=3)
survdiff(tmpSurv.lowgrade ~ riskEnet, rho=0)










###############################################################################
#  Look now at only the Grade II                                              #
###############################################################################


lowgradePat <- c()
lowgradePat$surv <- gbmPat.rembrandt$surv[which(rembrandtPat$Grade==" II")]
lowgradePat$survTime <- gbmPat.rembrandt$survTime[which(rembrandtPat$Grade==" II")]

tmpSurv.lowgrade <- Surv(lowgradePat$survTime,lowgradePat$surv)


yhatEnet <- predict(fitEnet, t(rembrandtEset.culled[,which(rembrandtPat$Grade==" II")]), type="response", s="lambda.min") 

boxplot(yhatEnet ~ predictedClasses.rembrandt[which(rembrandtPat$Grade==" II")], ylab="3-year OS prediction (%)", xlab="predicted POSTN group", main="elastic net validation")
stripchart(yhatEnet ~ predictedClasses.rembrandt[which(rembrandtPat$Grade==" II")],pch=20, col="royalblue", vertical=TRUE, add=TRUE, cex=.6)

rocEnet <- roc(predictor=as.numeric(yhatEnet), response=as.numeric(predictedClasses.rembrandt[which(rembrandtPat$Grade==" II")]),ci=TRUE)
plot.roc(rocEnet,col="red")

riskEnet <- as.vector(ifelse(yhatEnet >= 0.15, 1, 0))
names(riskEnet) <- rownames(yhatEnet)



ggkm(survfit(tmpSurv.lowgrade  ~ riskEnet),timeby=12,main="KM of REMBRANDT")
summary(survfit(tmpSurv.lowgrade  ~ riskEnet))


plot(survfit(tmpSurv.lowgrade ~ riskEnet), main="elastic net model", xlab="months",ylab="probability of OS",col= c("blue","magenta"),lwd=3)
survdiff(tmpSurv.lowgrade ~ riskEnet, rho=0)




