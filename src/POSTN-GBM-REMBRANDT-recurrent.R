### Andrew Trister
### Sage Bionetworks
### Seattle, WA
### 20130905

### Use the populated data to look at survival with POSTN in only recurrent GBMs in REMBRANDT


require(ggplot2)
require(survival)
require(mixtools)
require(pROC)
source("./src/ggkm.R")

#let's now put together just the GBM from the REMBRANDT 


# first rescale the data  to have the same mean and variance 
# Justin's function to rescale the VS to get the same mean/var than the TS
normalize_to_X <- function(mean.x, sd.x, Y){
  m.y <- rowMeans(Y)
  sd.y <- apply(Y, 1, sd)
  Y.adj <- (Y - m.y) * sd.x / sd.y  + mean.x 
  Y.adj[sd.y == 0] <- mean.x[sd.y==0]
  Y.adj
}

## DEFINE generatePcPlot() FUNCTION
# Create an SVD object, make a dataframe out of the first two factors
# and plot it using ggplot
generatePcPlot <- function(fullMatrix){
  require(ggplot2)
  require(corpcor)
  ## Make a ggplot-friendly dataframe from the singular value decomposition
  svdObj <- fast.svd(fullMatrix)
  pcDF <- data.frame(svdObj$v[ , 1:2])
  colnames(pcDF) <- c('PrinComp1', 'PrinComp2')
  ## Plot it
  pcPlot <- ggplot(pcDF, aes(PrinComp1, PrinComp2)) +
    geom_point(aes(colour = factor(studyIndicator), size = 20)) +
    scale_size(guide = 'none')
}




#hist(rembrandtEset['10631_mt',which(rembrandtPat$Grade==" IV")])
temp <- sub("_eg","_mt",rownames(eset))
overlapGenes <- which(rownames(rembrandtEset) %in% temp)

tmp <- apply(eset,1,sd)
tmp1 <- normalize_to_X(rowMeans(eset),tmp,rembrandtEset[overlapGenes,])


# studyIndicator <- c(rep('TCGA', ncol(eset)),
#                     rep('REMBRANDT',length(intersect(which(rembrandtPat$Disease==" GBM"),which(rembrandtPat$Prior.Therapy.Surgery.Tumor.Histology!=" --")))))
# 


rembrandtEset.recurrent <- rembrandtEset[overlapGenes,intersect(which(rembrandtPat$Disease==" GBM"),which(rembrandtPat$Prior.Therapy.Surgery.Tumor.Histology!=" --"))]



# boxplot(rembrandtEset.recurrent['10631_mt',])

hist(rembrandtEset.recurrent['10631_mt',])

#and fit a model to it
POSTNclassRembrandt.recurrent <- normalmixEM(rembrandtEset.recurrent['10631_mt',],lambda=0.5,mu=c(6.5,11.5),sigma=1)
plot(POSTNclassRembrandt.recurrent,density=T)
summary(POSTNclassRembrandt.recurrent)

predictedClassesRembrandt.recurrent <- rep(1,length(POSTNclassRembrandt.recurrent$x))

for(i in 1:length(POSTNclassRembrandt.recurrent$x)) {
  if(POSTNclassRembrandt.recurrent$posterior[i,1] > POSTNclassRembrandt.recurrent$posterior[i,2])
    predictedClassesRembrandt.recurrent[i] <- 2
}




## let's make a survival object

gbmPatRembrandt.recurrent <- c()
gbmPatRembrandt.recurrent$survTime <- 30.42*gbmPat.rembrandt$survTime[intersect(which(rembrandtPat$Disease==" GBM"),which(rembrandtPat$Prior.Therapy.Surgery.Tumor.Histology!=" --"))]
gbmPatRembrandt.recurrent$surv <- gbmPat.rembrandt$surv[intersect(which(rembrandtPat$Disease==" GBM"),which(rembrandtPat$Prior.Therapy.Surgery.Tumor.Histology!=" --"))]

tmpSurvRembrandt.recurrent <- Surv(gbmPatRembrandt.recurrent$survTime,gbmPatRembrandt.recurrent$surv)


##let's look at survival based on class for POSTN
ggkm(survfit(tmpSurvRembrandt.recurrent ~ predictedClassesRembrandt.recurrent), ystratalabs = (c("POSTN High", "POSTN Low")), 
     timeby = 365,
     main = "GBM K-M Plot By POSTN Expression (REMBRANDT+TCGA)")










### Analyze!###

rembrandtEset.recurrent.culled <- rembrandtEset.recurrent[-grep("10631_mt",rownames(rembrandtEset.recurrent)),]


tmp <- apply(rembrandtEset.recurrent.culled,1,sd)

#get rid of the least variant probes
tmp1 <- which(tmp>quantile(tmp,probs=0.2))
rembrandtEset.recurrent.culled <- rembrandtEset.recurrent.culled[tmp1,]
rm(tmp,tmp1)





cv.fit <- cv.glmnet(x=t(rembrandtEset.recurrent.culled), y=factor(predictedClassesRembrandt.recurrent), nfolds=10, alpha=.1, family="binomial")
plot(cv.fit)
fitEnet <- glmnet(x=t(rembrandtEset.recurrent.culled), y=factor(predictedClassesRembrandt.recurrent), family="binomial", alpha=.1, lambda=cv.fit$lambda.min)

# 
# cv.fit <- cv.glmnet(x=t(allGBMeset.culled), y=allGBMeset['10631_eg',], nfolds=10,alpha=0.1)
# plot(cv.fit)
# fitEnet <- glmnet(x=t(allGBMeset.culled), y=allGBMeset['10631_eg',], alpha=.1, lambda=cv.fit$lambda.min)


gene.names.recurrent <- lapply(rownames(rembrandtEset.recurrent.culled)[which(abs(fitEnet$beta)>0)],function(x){
  return(xx[strsplit(x,"_mt")[[1]]])
})

paste(unlist(gene.names.recurrent),collapse=" ")




#  build a more robust classifier( N= 500 bootstraps)
M <- 500
tmp <- c()
try.cv.fit <- c()
tryfit <- c()
for(i in c(1:M)) {
  print(i)
  N <- sample(1:length(colnames(rembrandtEset.recurrent.culled)),replace=TRUE)
  try.cv.fit <- cv.glmnet(x=t(rembrandtEset.recurrent.culled[,N]), y=factor(predictedClassesRembrandt.recurrent[N]), nfolds=10, alpha=.1, family="binomial")
  tryfit <- as.numeric(glmnet(x=t(rembrandtEset.recurrent.culled[,N]), y=factor(predictedClassesRembrandt.recurrent[N]), family="binomial", alpha=.1, lambda=try.cv.fit$lambda.min)$beta)
  tmp <- cbind(tmp,tryfit)
}

tmp <- apply(tmp,1,function(x){length(which(abs(x)>0))})
bootsave.rembrandt.recurrent <- tmp

select_features <- which(tmp>=quantile(tmp,probs=.999))

select_features <- which(tmp>=quantile(tmp,probs=.99))

# select_features <- which(tmp>=quantile(tmp,probs=.995))

# w <- as.data.frame(t(allGBMeset.culled[select_features,]))



# temp <- apply(allGBMeset.primary.culled,1,sd)
# tempEset <- normalize_to_X(rowMeans(allGBMeset.primary.culled),temp,rembrandtEset.scaled)
# 
# rembrandtEset.primary.boot <- tempEset[select_features,]

# gene.names.boot <- lapply(rownames(eset)[which(tmp==500)],function(x){
#   return(xx[strsplit(x,"_eg")[[1]]])
# })

gene.names.rembrandt.recurrent <- lapply(rownames(rembrandtEset.recurrent.culled)[select_features],function(x){
  return(xx[strsplit(x,"_mt")[[1]]])
})

paste(unlist(gene.names.rembrandt.recurrent[order(tmp[select_features],decreasing=T)]),collapse=" ")
tmp[select_features][order(tmp[select_features],decreasing=T)]





# 
#make a bootEnet
cv.fit <- cv.glmnet(x=t(rembrandtEset.recurrent.culled[select_features,]), y=factor(predictedClassesRembrandt.recurrent), nfolds=10, alpha=.1, family="binomial")
plot(cv.fit)
bootFitEnet <- glmnet(x=t(rembrandtEset.recurrent.culled[select_features,]), y=factor(predictedClassesRembrandt.recurrent), family="binomial", alpha=.1, lambda=cv.fit$lambda.min)

###### 
# 
# 
# # boostFitEnet <- glmnet(x=t(allGBMeset.culled[select_features,]), y=factor(predictedClassesAll), family="binomial", alpha=.1)
# # boostEnetfit <- glm(w$POSTN~. ,data = w, family = "binomial")
# rm(w)
# 
# 
# 
# #determine the correlation of the selected genes
# bootFitEnet$beta
# 
# 


# yhatboostEnet <- predict(object=boostEnetfit, newdata=as.data.frame(t(z)), type="response")
# 
# boxplot(yhatboostEnet ~ zhuClin$os3yr, ylab="prediction of 3-year OS probability (%)", xlab="3-year OS", main="elastic net - molecular + clinical features")
# stripchart(yhatboostEnet ~ zhuClin$os3yr,pch=20, col="royalblue", vertical=TRUE, add=TRUE, cex=.6)

# 
# 
# 
# 
# 
# gene.names <- lapply(rownames(eset)[which(abs(fitEnet$beta)>0)],function(x){
#   return(xx[strsplit(x,"_eg")[[1]]])
# })
# 
# paste(unlist(gene.names),collapse=" ")
# 
# 
# 
# 
# design <- model.matrix(~0+factor(predictedClasses))
# #design <- model.matrix(~age.group)
# colnames(design) = c("low", "high")
# 
# fit <- lmFit(eset.culled2, design)
# fit <- eBayes(fit)
# topTable(fit, coef=2)
# 
# contrast.matrix <- makeContrasts(low-high, levels=design)
# fit2 <- contrasts.fit(fit,contrast.matrix)
# fit2 <- eBayes(fit2)
# topTable(fit2,adjust="BH")#,lfc=0.5)
# results <- decideTests(fit2,p.value=0.05,lfc=0.5)
# table(results)
# 
# gene.names.selected <- lapply(rownames(eset.culled2)[which(results!=0)],function(x){
#   return(xx[strsplit(x,"_eg")[[1]]])
# })
# 
# paste(unlist(gene.names.selected),collapse=" ")
# 
# paste(unlist(intersect(gene.names,gene.names.selected)),collapse=" ")
# 


######
## See if this can predict the recurrence in Phillips
#####

tempEset <- cbind(phillipsEset.scaled[,matchedPrimary],phillipsEset.scaled[,recurrent])
temp <- apply(rembrandtEset.recurrent.culled,1,sd)
tempEset <- normalize_to_X(rowMeans(rembrandtEset.recurrent.culled),temp,tempEset)
tempClasses <- c(rep(1,23),rep(2,23))
phillipsEset.paired.boot <- tempEset[select_features,]

yhatEnet <- predict(bootFitEnet, t(phillipsEset.paired.boot), type="response",s="lambda.min")  #all of REMBRANDT


boxplot(yhatEnet ~ tempClasses, ylab="Predicted Recurrence", xlab="Recurrence or Primary", main="elastic net validation")
stripchart(yhatEnet ~ tempClasses,pch=20, col="royalblue", vertical=TRUE, add=TRUE, cex=.6)

rocEnet <- roc(predictor=as.numeric(yhatEnet), response=as.numeric(tempClasses-1),ci=TRUE)
plot.roc(rocEnet,col="red")
# 
# riskEnet <- as.vector(ifelse(yhatEnet >= median(yhatEnet), 1, 0))
# names(riskEnet) <- rownames(yhatEnet)
# 
# 
# plot(survfit(tmpSurv.rembrandt  ~ riskEnet), main="elastic net model", xlab="months",ylab="probability of OS",col= c("blue","magenta"),lwd=3)
# 
# survdiff(tmpSurv.rembrandt  ~ riskEnet, rho=0)




###############################################################################
#  Look now at all patients in REMBRANDT                                      #
###############################################################################


yhatEnet <- predict(fitEnet, t(rembrandtEset.scaled), type="response", s="lambda.min")  #all of REMBRANDT
yhatEnet <- predict(bootFitEnet, t(rembrandtEset.boot), type="response",s="lambda.min")  #all of REMBRANDT


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


yhatEnet <- predict(fitEnet, t(rembrandtEset.scaled[,which(rembrandtPat$Grade==" III")]), type="response", s="lambda.min")
yhatEnet <- predict(bootFitEnet, t(rembrandtEset.boot[,which(rembrandtPat$Grade==" III")]), type="response",s="lambda.min")  #all of REMBRANDT


boxplot(yhatEnet ~ predictedClasses.rembrandt[which(rembrandtPat$Grade==" III")], ylab="Predicted group", xlab="POSTN group", main="elastic net validation")
stripchart(yhatEnet ~ predictedClasses.rembrandt[which(rembrandtPat$Grade==" III")],pch=20, col="royalblue", vertical=TRUE, add=TRUE, cex=.6)

rocEnet <- roc(predictor=as.numeric(yhatEnet), response=as.numeric(predictedClasses.rembrandt[which(rembrandtPat$Grade==" III")]),ci=TRUE)
plot.roc(rocEnet,col="red")

riskEnet <- as.vector(ifelse(yhatEnet >= median(yhatEnet), 1, 0))
names(riskEnet) <- rownames(yhatEnet)



#ggkm(survfit(tmpSurv.lowgrade  ~ riskEnet),timeby=12,main="KM of REMBRANDT")
# summary(survfit(tmpSurv.lowgrade  ~ riskEnet))

plot(survfit(tmpSurv.lowgrade ~ riskEnet), main="elastic net model Grade III", xlab="months",ylab="probability of OS",col= c("blue","magenta"),lwd=3)
survdiff(tmpSurv.lowgrade ~ riskEnet, rho=0)










###############################################################################
#  Look now at only the Grade II                                              #
###############################################################################


lowgradePat <- c()
lowgradePat$surv <- gbmPat.rembrandt$surv[which(rembrandtPat$Grade==" II")]
lowgradePat$survTime <- gbmPat.rembrandt$survTime[which(rembrandtPat$Grade==" II")]

tmpSurv.lowgrade <- Surv(lowgradePat$survTime,lowgradePat$surv)


yhatEnet <- predict(fitEnet, t(rembrandtEset.scaled[,which(rembrandtPat$Grade==" II")]), type="response", s="lambda.min") 
yhatEnet <- predict(bootFitEnet, t(rembrandtEset.boot[,which(rembrandtPat$Grade==" II")]), type="response",s="lambda.min") 


boxplot(yhatEnet ~ predictedClasses.rembrandt[which(rembrandtPat$Grade==" II")], ylab="Predicted group", xlab="POSTN group", main="elastic net validation")
stripchart(yhatEnet ~ predictedClasses.rembrandt[which(rembrandtPat$Grade==" II")],pch=20, col="royalblue", vertical=TRUE, add=TRUE, cex=.6)

rocEnet <- roc(predictor=as.numeric(yhatEnet), response=as.numeric(predictedClasses.rembrandt[which(rembrandtPat$Grade==" II")]),ci=TRUE)
plot.roc(rocEnet,col="red")

riskEnet <- as.vector(ifelse(yhatEnet >= 0.15, 1, 0))
names(riskEnet) <- rownames(yhatEnet)



#ggkm(survfit(tmpSurv.lowgrade  ~ riskEnet),timeby=12,main="KM of REMBRANDT")
# summary(survfit(tmpSurv.lowgrade  ~ riskEnet))


plot(survfit(tmpSurv.lowgrade ~ riskEnet), main="elastic net model Grade II", xlab="months",ylab="probability of OS",col= c("blue","magenta"),lwd=3)
survdiff(tmpSurv.lowgrade ~ riskEnet, rho=0)
#survdiff(tmpSurv.lowgrade ~ rembrandtEset["10631_mt",which(rembrandtPat$Grade==" II")]>8.75, rho=0)





###############################################################################
#  Look now at only the Grade II and III                                      #
###############################################################################


lowgradePat <- c()
lowgradePat$surv <- gbmPat.rembrandt$surv[which(rembrandtPat$Grade!=" IV")]
lowgradePat$survTime <- gbmPat.rembrandt$survTime[which(rembrandtPat$Grade!=" IV")]


tmpSurv.lowgrade <- Surv(lowgradePat$survTime,lowgradePat$surv)


yhatEnet <- predict(fitEnet, t(rembrandtEset.scaled[,which(rembrandtPat$Grade!=" IV")]), type="response", s="lambda.min")
yhatEnet <- predict(bootFitEnet, t(rembrandtEset.boot[,which(rembrandtPat$Grade!=" IV")]), type="response",s="lambda.min") 


boxplot(yhatEnet ~ predictedClasses.rembrandt[which(rembrandtPat$Grade!=" IV")], ylab="Predicted group", xlab="POSTN group", main="elastic net validation")
stripchart(yhatEnet ~ predictedClasses.rembrandt[which(rembrandtPat$Grade!=" IV")],pch=20, col="royalblue", vertical=TRUE, add=TRUE, cex=.6)

rocEnet <- roc(predictor=as.numeric(yhatEnet), response=as.numeric(predictedClasses.rembrandt[which(rembrandtPat$Grade!=" IV")]),ci=TRUE)
plot.roc(rocEnet,col="red")

riskEnet <- as.vector(ifelse(yhatEnet >= median(yhatEnet), 1, 0))

# riskEnet <- as.vector(ifelse(yhatEnet >= 0.02, 1, 0))
names(riskEnet) <- rownames(yhatEnet)



ggkm(survfit(tmpSurv.lowgrade  ~ riskEnet),timeby=12,main="KM of REMBRANDT")
# summary(survfit(tmpSurv.lowgrade  ~ riskEnet))

plot(survfit(tmpSurv.lowgrade ~ riskEnet), main="elastic net model Grade II and III", xlab="months",ylab="probability of OS",col= c("blue","magenta"),lwd=3)
survdiff(tmpSurv.lowgrade ~ riskEnet, rho=0)










heatmap(x=rembrandtEset.culled[which(abs(fitEnet$beta)>0),],
        ColSideColors=c("red","green","blue","yellow")[factor(rembrandtPat$Disease)],
        scale="row",col=redgreen(50),main="Elastic Net")



gene.names <- lapply(rownames(eset)[which(abs(fitEnet$beta)>0)],function(x){
  return(xx[strsplit(x,"_eg")[[1]]])
})

paste(unlist(gene.names),collapse=" ")

# 
# 
# ##
# # Look at the Phillips data
# 
# phillipsEset.culled <- phillipsEset[-grep("10631_mt",rownames(phillipsEset)),recurrent]
# phillipsEset.culled <- cbind(phillipsEset.culled,phillipsEset[-grep("10631_mt",rownames(phillipsEset)),matchedPrimary])
# 
# yhatEnet <- predict(fitEnet, t(phillipsEset.culled), type="response", s="lambda.min")
# # boxplot(yhatEnet ~ c(rep(1,23),rep(0,23)), ylab="3-year OS prediction (%)", xlab="predicted POSTN group", main="elastic net validation")
# # stripchart(yhatEnet ~ c(rep(1,23),rep(0,23)),pch=20, col="royalblue", vertical=TRUE, add=TRUE, cex=.6)
# 
# 
# 
# riskEnet <- as.vector(ifelse(yhatEnet >= median(yhatEnet), 1, 0))
# names(riskEnet) <- rownames(yhatEnet)
# 
# 
# plot(survfit(tmpSurv.phillipsAll  ~ riskEnet[1:23]), main="elastic net model", xlab="months",ylab="probability of OS",col= c("blue","magenta"),lwd=3)
# survdiff(tmpSurv.phillipsAll  ~ riskEnet, rho=0)
# 








###############################################################################
#  Look now at all patients in Phillips                                       #
###############################################################################


temp <- apply(allGBMeset.culled,1,sd)
tempEset <- normalize_to_X(rowMeans(allGBMeset.culled),temp,phillipsEset.scaled[-grep("10631",rownames(phillipsEset.scaled)),])
templist <- grep("primary",phillipsMetadata$X.Sample_characteristics_ch1)
phillipsEset.boot <- tempEset[select_features,templist]


# yhatEnet <- predict(fitEnet, t(rembrandtEset.scaled), type="response", s="lambda.min")  #all of REMBRANDT
yhatEnet <- predict(bootFitEnet, t(phillipsEset.boot), type="response",s="lambda.min")  #all of REMBRANDT


# boxplot(yhatEnet ~ predictedClasses.rembrandt, ylab="3-year OS prediction (%)", xlab="predicted POSTN group", main="elastic net validation")
# stripchart(yhatEnet ~ predictedClasses.rembrandt,pch=20, col="royalblue", vertical=TRUE, add=TRUE, cex=.6)

# rocEnet <- roc(predictor=as.numeric(yhatEnet), response=as.numeric(predictedClasses.rembrandt-1),ci=TRUE)
# plot.roc(rocEnet,col="red")

# riskEnet <- as.vector(ifelse(yhatEnet >= median(yhatEnet), 1, 0))
riskEnet <- as.vector(ifelse(yhatEnet >= quantile(yhatEnet)[2], 1, 0))

names(riskEnet) <- rownames(yhatEnet)


plot(survfit(tmpSurv.phillipsAll  ~ riskEnet), main="elastic net model", xlab="weeks",ylab="probability of OS",col= c("blue","magenta"),lwd=3)

ggkm(survfit(tmpSurv.phillipsAll  ~ riskEnet),timeby=52,main="KM of Phillips")

survdiff(tmpSurv.phillipsAll  ~ riskEnet, rho=0)





### Build a more robust classifer using bootstrapping

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
    ResultBS<-cbind(ResultBS,rank(abs(as.numeric(resultsModel[[k]])),ties.method="min")/length(resultsModel[[k]]))
  }
  rownames(ResultBS)<-rownames(resultsModel[[1]])  
  reference <- apply(ResultBS,1,sum) 
  return(reference)
}




selected <- weightAggregation(features)

plot(sort(selected))
abline(v=length(selected)-100)

#use the weighted aggregation to grab the top 100
select_features <- sort(selected,decreasing=TRUE)[1:100]

#find the correlations between the selected features and POSTN
signFeatures <- unlist(lapply(names(select_features),function(y){
  return(lm(eset["10631_eg",]~eset[y,])$coefficients[2])
}))


temp <- names(sort(select_features)[100:1])
temp <- names(select_features)[which(signFeatures<=0)]

gene.names <- lapply(temp,function(x){
  return(xx[strsplit(x,"_eg")[[1]]])
})

paste(unlist(gene.names),collapse=" ")


heatmap(x=rembrandtEset.scaled[names(sort(select_features)[90:100]),],
        ColSideColors=c("red","green","blue","yellow")[factor(rembrandtPat$Disease)],
        scale="none",col=redgreen(50),main="Elastic Net")


###############################################################################
#  Look now at only the Grade III                                             #
###############################################################################


lowgradePat <- c()
lowgradePat$surv <- gbmPat.rembrandt$surv[which(rembrandtPat$Grade==" III")]
lowgradePat$survTime <- gbmPat.rembrandt$survTime[which(rembrandtPat$Grade==" III")]


tmpSurv.lowgrade <- Surv(lowgradePat$survTime,lowgradePat$surv)


yhatEnet <- predict(fitEnet, t(rembrandtEset.scaled[,which(rembrandtPat$Grade==" III")]), type="response", s="lambda.min")

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





