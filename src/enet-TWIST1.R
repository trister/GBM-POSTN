### Andrew Trister
### Sage Bionetworks
### Seattle, WA
### 20121119

### Build an elasticnet model trained on TWIST1 class in TCGA to learn what other genes may be implicated
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




#call median value of TWIST1 a cutoff

temp <- median(eset.culled['7291_eg',])
TWISTclass <- eset.culled['7291_eg',]

for(i in 1:length(TWISTclass)) {
  if(TWISTclass[i] < temp) {
    TWISTclass[i] <- 0
  } else TWISTclass[i] <- 1
}


temp <- median(rembrandtEset.scaled['7291_eg',])
TWISTclass.rembrandt <- rembrandtEset.scaled['7291_eg',]

for(i in 1:length(TWISTclass.rembrandt)) {
  if(TWISTclass.rembrandt[i] < temp) {
    TWISTclass.rembrandt[i] <- 0
  } else TWISTclass.rembrandt[i] <- 1
}


# 
# #cull the REMBRANDT U133Aplus2 to U133A genesets
# temp <- unlist(lapply(rownames(rembrandtEset),function(x){
#   return(gsub("_mt","_eg",x))
# }))
# 
# rembrandtEset.culled <- rembrandtEset
# rownames(rembrandtEset.culled) <- temp
# 
# temp <- intersect(temp,rownames(eset))
# 
# eset <- eset[temp,]
# rembrandtEset.culled <- rembrandtEset.culled[temp,]
# 
# eset.culled <- eset[-grep("10631_eg",rownames(eset)),]
# rembrandtEset.culled <- rembrandtEset.culled[-grep("10631_eg",rownames(rembrandtEset.culled)),]


# run the elastic net model
#  cv.fit <- cv.glmnet(x=t(eset.culled), y=eset["10631_eg",], nfolds=10, alpha=.1, family="gaussian")
#  plot(cv.fit)
#  fitEnet <- glmnet(x=t(eset.culled), y=eset["10631_eg",], family="gaussian", alpha=.1, lambda=cv.fit$lambda.min)
# 
# cv.fit <- cv.glmnet(x=t(eset), y=factor(predictedClasses), nfolds=10, alpha=.1, family="binomial")
# plot(cv.fit)
# fitEnet <- glmnet(x=t(eset), y=factor(predictedClasses), family="binomial", alpha=.1, lambda=cv.fit$lambda.min)


cv.fit <- cv.glmnet(x=t(eset.culled), y=factor(TWISTclass), nfolds=10, alpha=.1, family="binomial")
plot(cv.fit)
fitEnet <- glmnet(x=t(eset.culled), y=factor(TWISTclass), family="binomial", alpha=.1, lambda=cv.fit$lambda.min)




###############################################################################
#  Look now at all patients in REMBRANDT                                      #
###############################################################################


yhatEnet <- predict(fitEnet, t(rembrandtEset.scaled), type="response", s="lambda.min")  #all of REMBRANDT

boxplot(yhatEnet ~ TWISTclass.rembrandt, ylab="3-year OS prediction (%)", xlab="predicted TWIST group", main="elastic net validation")
stripchart(yhatEnet ~ TWISTclass.rembrandt,pch=20, col="royalblue", vertical=TRUE, add=TRUE, cex=.6)

rocEnet <- roc(predictor=as.numeric(yhatEnet), response=as.numeric(TWISTclass.rembrandt),ci=TRUE)
plot.roc(rocEnet,col="red")

riskEnet <- as.vector(ifelse(yhatEnet >= median(yhatEnet), 1, 0))
names(riskEnet) <- rownames(yhatEnet)


plot(survfit(tmpSurv.rembrandt  ~ riskEnet), main="elastic net model", xlab="months",ylab="probability of OS",col= c("blue","magenta"),lwd=3)
survdiff(tmpSurv.rembrandt  ~ riskEnet, rho=0)

ggkm(survfit(tmpSurv.rembrandt ~ riskEnet),timeby=12,main="KM of REMBRANDT",
     ystratalabs=c("high TWIST","low TWIST"))
#summary(survfit(tmpSurv.lowgrade  ~ riskEnet))








###############################################################################
#  Look now at only the Grade IV                                              #
###############################################################################


lowgradePat <- c()
lowgradePat$surv <- gbmPat.rembrandt$surv[which(rembrandtPat$Grade==" IV")]
lowgradePat$survTime <- gbmPat.rembrandt$survTime[which(rembrandtPat$Grade==" IV")]


tmpSurv.lowgrade <- Surv(lowgradePat$survTime,lowgradePat$surv)


yhatEnet <- predict(fitEnet, t(rembrandtEset.scaled[,which(rembrandtPat$Grade==" IV")]), type="response", s="lambda.min")

boxplot(yhatEnet ~ predictedClasses.rembrandt[which(rembrandtPat$Grade==" IV")], ylab="3-year OS prediction (%)", xlab="predicted TWIST group", main="elastic net validation")
stripchart(yhatEnet ~ predictedClasses.rembrandt[which(rembrandtPat$Grade==" IV")],pch=20, col="royalblue", vertical=TRUE, add=TRUE, cex=.6)

rocEnet <- roc(predictor=as.numeric(yhatEnet), response=as.numeric(predictedClasses.rembrandt[which(rembrandtPat$Grade==" III")]),ci=TRUE)
plot.roc(rocEnet,col="red")

riskEnet <- as.vector(ifelse(yhatEnet >= median(yhatEnet), 1, 0))
names(riskEnet) <- rownames(yhatEnet)



ggkm(survfit(tmpSurv.lowgrade  ~ riskEnet),timeby=12,main="KM of REMBRANDT")
summary(survfit(tmpSurv.lowgrade  ~ riskEnet))

plot(survfit(tmpSurv.lowgrade ~ riskEnet), main="elastic net model", xlab="months",ylab="probability of OS",col= c("blue","magenta"),lwd=3)
survdiff(tmpSurv.lowgrade ~ riskEnet, rho=0)







###############################################################################
#  Look now at only the Grade III                                             #
###############################################################################


lowgradePat <- c()
lowgradePat$surv <- gbmPat.rembrandt$surv[which(rembrandtPat$Grade==" III")]
lowgradePat$survTime <- gbmPat.rembrandt$survTime[which(rembrandtPat$Grade==" III")]


tmpSurv.lowgrade <- Surv(lowgradePat$survTime,lowgradePat$surv)


yhatEnet <- predict(fitEnet, t(rembrandtEset.scaled[,which(rembrandtPat$Grade==" III")]), type="response", s="lambda.min")

boxplot(yhatEnet ~ predictedClasses.rembrandt[which(rembrandtPat$Grade==" III")], ylab="3-year OS prediction (%)", xlab="predicted TWIST group", main="elastic net validation")
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


yhatEnet <- predict(fitEnet, t(rembrandtEset.scaled[,which(rembrandtPat$Grade==" II")]), type="response", s="lambda.min") 

boxplot(yhatEnet ~ predictedClasses.rembrandt[which(rembrandtPat$Grade==" II")], ylab="3-year OS prediction (%)", xlab="predicted TWIST group", main="elastic net validation")
stripchart(yhatEnet ~ predictedClasses.rembrandt[which(rembrandtPat$Grade==" II")],pch=20, col="royalblue", vertical=TRUE, add=TRUE, cex=.6)

rocEnet <- roc(predictor=as.numeric(yhatEnet), response=as.numeric(predictedClasses.rembrandt[which(rembrandtPat$Grade==" II")]),ci=TRUE)
plot.roc(rocEnet,col="red")

riskEnet <- as.vector(ifelse(yhatEnet >= median(yhatEnet), 1, 0))
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








#  build a more robust classifier( bootstrapping)
df <- data.frame(value=TWISTclass-1)
rownames(df) <- colnames(eset.culled)

tmp <- c()
try.cv.fit <- c()
tryfit <- c()
for(i in c(1:100)) {
  print(i)
  N <- sample(colnames(eset.culled),559,replace=TRUE)
  
#   try.cv.fit <- cv.glmnet(x=t(eset.culled[,N]), y=eset["10631_eg",N], nfolds=10, alpha=.1, family="gaussian")
#   tryfit <- as.numeric(glmnet(x=t(eset.culled[,N]), y=eset["10631_eg",N], family="gaussian", alpha=0.1,lambda=try.cv.fit$lambda.min)$beta)
  
    try.cv.fit <- cv.glmnet(x=t(eset.culled[,N]), y=factor(df[N,]), nfolds=10, alpha=0.1, family="binomial")
    tryfit <- as.numeric(glmnet(x=t(eset.culled[,N]), y=factor(df[N,]), family="binomial", alpha=0.1, lambda=try.cv.fit$lambda.min)$beta)
  tmp <- cbind(tmp,tryfit)
}

temp <- apply(tmp,1,function(x){length(which(abs(x)>0))})
select_features <- which(temp>=quantile(temp,probs=.99))

w <- as.data.frame(t(eset[select_features,]))
w$predictedClass <- predictedClasses-1

# run the logit model on the top selected features
bootEnetfit <- glm(w$predictedClass~. ,data = w, family = binomial())
rm(w)

yhatbootEnet <- predict(object=bootEnetfit, newdata=as.data.frame(t(rembrandtEset.culled)), type="response")

boxplot(yhatbootEnet ~ (TWISTclass.rembrandt), ylab="prediction of TWIST expression", xlab="TWIST", main="bootstrap elastic net - molecular features")
stripchart(yhatbootEnet ~ (TWISTclass.rembrandt),pch=20, col="royalblue", vertical=TRUE, add=TRUE, cex=.6)




gene.names <- lapply(rownames(eset)[select_features],function(x){
  return(xx[strsplit(x,"_eg")[[1]]])
})

paste(unlist(gene.names),collapse=" ")







