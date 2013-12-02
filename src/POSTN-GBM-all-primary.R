### Andrew Trister
### Sage Bionetworks
### Seattle, WA
### 20130905

### Use the populated data to look at survival with POSTN in all PRIMARY GBMS in TCGA and REMBRANDT


require(ggplot2)
require(survival)
require(mixtools)
require(glmnet)
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


studyIndicator <- c(rep('TCGA', ncol(eset)),
                    rep('REMBRANDT',length(intersect(which(rembrandtPat$Disease==" GBM"),which(rembrandtPat$Prior.Therapy.Surgery.Tumor.Histology==" --")))))



allGBMeset.primary <- cbind(eset,rembrandtEset[overlapGenes,intersect(which(rembrandtPat$Disease==" GBM"),which(rembrandtPat$Prior.Therapy.Surgery.Tumor.Histology==" --"))])


# Plot the raw expression data
rawPcPlot <- generatePcPlot(allGBMeset.primary) + 
  opts(title = 'Raw Unnormalized Matrix\n')
rawPcPlot


boxplot(rembrandtEset['10631_mt',intersect(which(rembrandtPat$Disease==" GBM"),which(rembrandtPat$Prior.Therapy.Surgery.Tumor.Histology==" --"))],
        rembrandtEset['10631_mt',intersect(which(rembrandtPat$Disease==" GBM"),which(rembrandtPat$Prior.Therapy.Surgery.Tumor.Histology!=" --"))])

allGBMeset.primary <- cbind(eset,tmp1[,intersect(which(rembrandtPat$Disease==" GBM"),which(rembrandtPat$Prior.Therapy.Surgery.Tumor.Histology==" --"))])


# Plot the normalized expression data
normalizedPcPlot <- generatePcPlot(allGBMeset.primary) + 
  opts(title = 'Normalized Matrix\n')
normalizedPcPlot



hist(allGBMeset.primary['10631_eg',])

#and fit a model to it
POSTNclassAll.primary <- normalmixEM(allGBMeset.primary['10631_eg',],lambda=0.5,mu=c(6.5,11.5),sigma=1)
plot(POSTNclassAll.primary,density=T)
summary(POSTNclassAll.primary)

predictedClassesAll.primary <- rep(1,length(POSTNclassAll.primary$x))

for(i in 1:length(POSTNclassAll.primary$x)) {
  if(POSTNclassAll.primary$posterior[i,1] > POSTNclassAll.primary$posterior[i,2])
    predictedClassesAll.primary[i] <- 2
}




## let's make a survival object

gbmPatAll.primary <- c()
gbmPatAll.primary$survTime <- c(gbmPat$survTime,30.42*gbmPat.rembrandt$survTime[intersect(which(rembrandtPat$Disease==" GBM"),which(rembrandtPat$Prior.Therapy.Surgery.Tumor.Histology==" --"))])
gbmPatAll.primary$surv <- c(gbmPat$surv,gbmPat.rembrandt$surv[intersect(which(rembrandtPat$Disease==" GBM"),which(rembrandtPat$Prior.Therapy.Surgery.Tumor.Histology==" --"))])

tmpSurvAll.primary <- Surv(gbmPatAll.primary$survTime,gbmPatAll.primary$surv)


##let's look at survival based on class for POSTN
ggkm(survfit(tmpSurvAll.primary ~ predictedClassesAll.primary), ystratalabs = (c("POSTN High", "POSTN Low")), 
     timeby = 365,
     main = "GBM K-M Plot By POSTN Expression (REMBRANDT+TCGA)")










### Analyze!###

allGBMeset.primary.culled <- allGBMeset.primary[-grep("10631_eg",rownames(allGBMeset.primary)),]


tmp <- apply(allGBMeset.primary.culled,1,sd)

#get rid of the least variant probes
tmp1 <- which(tmp>quantile(tmp,probs=0.2))
allGBMeset.primary.culled <- allGBMeset.primary.culled[tmp1,]
rm(tmp,tmp1)





cv.fit <- cv.glmnet(x=t(allGBMeset.primary.culled), y=factor(predictedClassesAll.primary), nfolds=10, alpha=.1, family="binomial")
plot(cv.fit)
fitEnet <- glmnet(x=t(allGBMeset.primary.culled), y=factor(predictedClassesAll.primary), family="binomial", alpha=.1, lambda=cv.fit$lambda.min)

# 
# cv.fit <- cv.glmnet(x=t(allGBMeset.culled), y=allGBMeset['10631_eg',], nfolds=10,alpha=0.1)
# plot(cv.fit)
# fitEnet <- glmnet(x=t(allGBMeset.culled), y=allGBMeset['10631_eg',], alpha=.1, lambda=cv.fit$lambda.min)


gene.names <- lapply(rownames(allGBMeset.primary.culled)[which(abs(fitEnet$beta)>0)],function(x){
  return(xx[strsplit(x,"_eg")[[1]]])
})

paste(unlist(gene.names),collapse=" ")




#  build a more robust classifier( N= 500 bootstraps)
M <- 500
tmp <- c()
try.cv.fit <- c()
tryfit <- c()
for(i in c(1:M)) {
  print(i)
  N <- sample(1:length(colnames(allGBMeset.primary.culled)),replace=TRUE)
  try.cv.fit <- cv.glmnet(x=t(allGBMeset.primary.culled[,N]), y=factor(predictedClassesAll.primary[N]), nfolds=10, alpha=.1, family="binomial")
  tryfit <- as.numeric(glmnet(x=t(allGBMeset.primary.culled[,N]), y=factor(predictedClassesAll.primary[N]), family="binomial", alpha=.1, lambda=try.cv.fit$lambda.min)$beta)
  tmp <- cbind(tmp,tryfit)
}

tmp <- apply(tmp,1,function(x){length(which(abs(x)>0))})
bootsave.primary <- tmp
select_features <- which(tmp>=quantile(tmp,probs=.999))

select_features <- which(tmp>=quantile(tmp,probs=.99))

select_features <- which(tmp>=quantile(tmp,probs=.995))

# w <- as.data.frame(t(allGBMeset.culled[select_features,]))



temp <- apply(allGBMeset.primary.culled,1,sd)
tempEset <- normalize_to_X(rowMeans(allGBMeset.primary.culled),temp,rembrandtEset.scaled)

rembrandtEset.primary.boot <- tempEset[select_features,]

# gene.names.boot <- lapply(rownames(eset)[which(tmp==500)],function(x){
#   return(xx[strsplit(x,"_eg")[[1]]])
# })

gene.names.primary.boot <- lapply(rownames(allGBMeset.primary.culled)[select_features],function(x){
  return(xx[strsplit(x,"_eg")[[1]]])
})

paste(unlist(gene.names.primary.boot[order(tmp[select_features],decreasing=T)]),collapse=" ")
tmp[select_features][order(tmp[select_features],decreasing=T)]






#may a bootEnet
cv.fit <- cv.glmnet(x=t(allGBMeset.culled[select_features,]), y=factor(predictedClassesAll), nfolds=10, alpha=.1, family="binomial")
plot(cv.fit)
bootFitEnet <- glmnet(x=t(allGBMeset.culled[select_features,]), y=factor(predictedClassesAll), family="binomial", alpha=.1, lambda=cv.fit$lambda.min)

# 
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

## this is for grade III and IVs
templist <- grep("primary",phillipsMetadata$X.Sample_characteristics_ch1)
phillipsEset.boot <- tempEset[select_features,templist]


## now do only the grade IVs
# templist <- grep("primary",phillipsMetadata$X.Sample_characteristics_ch1)
# tempGrade <- unlist(lapply(phillipsMetadata$X.Sample_characteristics_ch1, function(x){
#   strsplit(x,";")[[1]][1]
# }))
# templist <- templist[!templist %in%which(tempGrade=="WHO grade: III")]
# # intersect(templist,which(tempGrade=="WHO grade: III"))




#now start the analysis
phillipsEset.boot <- tempEset[select_features,templist]
cutpoint <- quantile(phillipsEset.scaled['10631_mt',templist])[3] #median cutpoint

tempPredict <- rep(1,length(templist))

#use this for cuts with cutpoint 
for(i in 1:length(templist)) {
  if(phillipsEset.scaled["10631_mt",templist[i]] < cutpoint)
    tempPredict[i] <- 2
}


table(tempPredict)



phillipsSurvAll <- unlist(lapply(phillipsMetadata[templist,'X.Sample_characteristics_ch1'],function(x){
  temp <- strsplit(x,";")[[1]][5]
  temp <- strsplit(temp,": ")[[1]][2]
  if(temp=="no") 
    return(1)
  else return(0)
}))

phillipsSurvTimeAll <- unlist(lapply(phillipsMetadata[templist,'X.Sample_characteristics_ch1'],function(x){
  temp <- strsplit(x,";")[[1]][4]
  return(strsplit(temp,": ")[[1]][2])
}))
phillipsSurvTimeAll <- as.numeric(phillipsSurvTimeAll)

phillipsGradeAll <- unlist(lapply(phillipsMetadata[templist,'X.Sample_characteristics_ch1'],function(x){
  return(strsplit(x,";")[[1]][1])
}))

tmpSurv.phillipsAll <- Surv(phillipsSurvTimeAll,phillipsSurvAll)




# yhatEnet <- predict(fitEnet, t(rembrandtEset.scaled), type="response", s="lambda.min")  #all of REMBRANDT
yhatEnet <- predict(bootFitEnet, t(phillipsEset.boot), type="response",s="lambda.min")  #all of REMBRANDT
# yhatEnet <- predict(fitEnet, t(phillipsEset.boot), type="response",s="lambda.min")  #all of REMBRANDT

boxplot(yhatEnet~phillipsGradeAll=="WHO grade: III")
t.test(yhatEnet~phillipsGradeAll=="WHO grade: III")


# boxplot(yhatEnet ~ predictedClasses.rembrandt, ylab="3-year OS prediction (%)", xlab="predicted POSTN group", main="elastic net validation")
# stripchart(yhatEnet ~ predictedClasses.rembrandt,pch=20, col="royalblue", vertical=TRUE, add=TRUE, cex=.6)

# rocEnet <- roc(predictor=as.numeric(yhatEnet), response=as.numeric(phillipsGradeAll=="WHO grade: III"),ci=TRUE)
rocEnet <- roc(predictor=as.numeric(yhatEnet), response=as.numeric(tempPredict-1),ci=TRUE)
plot.roc(rocEnet,col="red")

#  riskEnet <- as.vector(ifelse(yhatEnet >= median(yhatEnet), 1, 0))
riskEnet <- as.vector(ifelse(yhatEnet >= quantile(yhatEnet)[2], 1, 0))

names(riskEnet) <- rownames(yhatEnet)


plot(survfit(tmpSurv.phillipsAll  ~ riskEnet), main="elastic net model", xlab="weeks",ylab="probability of OS",col= c("blue","magenta"),lwd=3)

ggkm(survfit(tmpSurv.phillipsAll  ~ riskEnet),timeby=52,main="KM of Phillips")

survdiff(tmpSurv.phillipsAll  ~ riskEnet, rho=0)












###############################################################################
#  Look now at all patients in Murat                                          #
###############################################################################


temp <- apply(allGBMeset.culled,1,sd)
tempEset <- normalize_to_X(rowMeans(allGBMeset.culled),temp,phillipsEset.scaled[-grep("10631",rownames(phillipsEset.scaled)),])

## this is for grade III and IVs
templist <- grep("primary",phillipsMetadata$X.Sample_characteristics_ch1)
phillipsEset.boot <- tempEset[select_features,templist]


## now do only the grade IVs
# templist <- grep("primary",phillipsMetadata$X.Sample_characteristics_ch1)
# tempGrade <- unlist(lapply(phillipsMetadata$X.Sample_characteristics_ch1, function(x){
#   strsplit(x,";")[[1]][1]
# }))
# templist <- templist[!templist %in%which(tempGrade=="WHO grade: III")]
# # intersect(templist,which(tempGrade=="WHO grade: III"))




#now start the analysis
phillipsEset.boot <- tempEset[select_features,templist]
cutpoint <- quantile(phillipsEset.scaled['10631_mt',templist])[3] #median cutpoint

tempPredict <- rep(1,length(templist))

#use this for cuts with cutpoint 
for(i in 1:length(templist)) {
  if(phillipsEset.scaled["10631_mt",templist[i]] < cutpoint)
    tempPredict[i] <- 2
}


table(tempPredict)



phillipsSurvAll <- unlist(lapply(phillipsMetadata[templist,'X.Sample_characteristics_ch1'],function(x){
  temp <- strsplit(x,";")[[1]][5]
  temp <- strsplit(temp,": ")[[1]][2]
  if(temp=="no") 
    return(1)
  else return(0)
}))

phillipsSurvTimeAll <- unlist(lapply(phillipsMetadata[templist,'X.Sample_characteristics_ch1'],function(x){
  temp <- strsplit(x,";")[[1]][4]
  return(strsplit(temp,": ")[[1]][2])
}))
phillipsSurvTimeAll <- as.numeric(phillipsSurvTimeAll)

phillipsGradeAll <- unlist(lapply(phillipsMetadata[templist,'X.Sample_characteristics_ch1'],function(x){
  return(strsplit(x,";")[[1]][1])
}))

tmpSurv.phillipsAll <- Surv(phillipsSurvTimeAll,phillipsSurvAll)




# yhatEnet <- predict(fitEnet, t(rembrandtEset.scaled), type="response", s="lambda.min")  #all of REMBRANDT
yhatEnet <- predict(bootFitEnet, t(phillipsEset.boot), type="response",s="lambda.min")  #all of REMBRANDT
# yhatEnet <- predict(fitEnet, t(phillipsEset.boot), type="response",s="lambda.min")  #all of REMBRANDT

boxplot(yhatEnet~phillipsGradeAll=="WHO grade: III")
t.test(yhatEnet~phillipsGradeAll=="WHO grade: III")


# boxplot(yhatEnet ~ predictedClasses.rembrandt, ylab="3-year OS prediction (%)", xlab="predicted POSTN group", main="elastic net validation")
# stripchart(yhatEnet ~ predictedClasses.rembrandt,pch=20, col="royalblue", vertical=TRUE, add=TRUE, cex=.6)

# rocEnet <- roc(predictor=as.numeric(yhatEnet), response=as.numeric(phillipsGradeAll=="WHO grade: III"),ci=TRUE)
rocEnet <- roc(predictor=as.numeric(yhatEnet), response=as.numeric(tempPredict-1),ci=TRUE)
plot.roc(rocEnet,col="red")

#  riskEnet <- as.vector(ifelse(yhatEnet >= median(yhatEnet), 1, 0))
riskEnet <- as.vector(ifelse(yhatEnet >= quantile(yhatEnet)[2], 1, 0))

names(riskEnet) <- rownames(yhatEnet)


plot(survfit(tmpSurv.phillipsAll  ~ riskEnet), main="elastic net model", xlab="weeks",ylab="probability of OS",col= c("blue","magenta"),lwd=3)

ggkm(survfit(tmpSurv.phillipsAll  ~ riskEnet),timeby=52,main="KM of Phillips")

survdiff(tmpSurv.phillipsAll  ~ riskEnet, rho=0)





