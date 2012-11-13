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
cv.fit <- cv.glmnet(x=t(eset.culled), y=eset["10631_eg",], nfolds=10, alpha=.1, family="gaussian")
plot(cv.fit)
fitEnet <- glmnet(x=t(eset.culled), y=eset["10631_eg",], family="gaussian", alpha=.1, lambda=cv.fit$lambda.min)


# cv.fit <- cv.glmnet(x=t(eset.culled), y=factor(predictedClasses), nfolds=10, alpha=.1, family="binomial")
# plot(cv.fit)
# fitEnet <- glmnet(x=t(eset.culled), y=factor(predictedClasses), family="binomial", alpha=.1, lambda=cv.fit$lambda.min)

yhatEnet <- predict(fitEnet, t(rembrandtEset.culled), type="response", s="lambda.min")
boxplot(yhatEnet ~ predictedClasses.rembrandt, ylab="3-year OS prediction (%)", xlab="predicted POSTN group", main="elastic net validation")
stripchart(yhatEnet ~ predictedClasses.rembrandt,pch=20, col="royalblue", vertical=TRUE, add=TRUE, cex=.6)

rocEnet <- roc(predictor=as.numeric(yhatEnet), response=as.numeric(predictedClasses.rembrandt-1),ci=TRUE)
plot.roc(rocEnet,col="red")

riskEnet <- as.vector(ifelse(yhatEnet >= median(yhatEnet), 1, 0))
names(riskEnet) <- rownames(yhatEnet)


plot(survfit(tmpSurv.rembrandt  ~ riskEnet), main="elastic net model", xlab="months",ylab="probability of OS",col= c("blue","magenta"),lwd=3)
survdiff(tmpSurv.rembrandt  ~ riskEnet, rho=0)



ggkm(survfit(tmpSurv.rembrandt  ~ riskEnet),timeby=12,main="KM of REMBRANDT")
ggkm(survfit(tmpSurv.rembrandt[which(rembrandtPat$Disease==" GBM")]~riskEnet[which(rembrandtPat$Disease==" GBM")]),timeby=12,
     main="KM of REMBRANDT GBM")

ggkm(survfit(tmpSurv.rembrandt[which(rembrandtPat$Disease==" ASTROCYTOMA")]~riskEnet[which(rembrandtPat$Disease==" ASTROCYTOMA")]),timeby=12,
     ystratalabs = c("Unfavorable","Favorable"),main="KM of REMBRANDT ASTROCYTOMA")

ggkm(survfit(tmpSurv.rembrandt[which(rembrandtPat$Disease==" OLIGODENDROGLIOMA")]~riskEnet[which(rembrandtPat$Disease==" OLIGODENDROGLIOMA")]),timeby=12,
     ystratalabs = c("Unfavorable","Favorable"),main="KM of REMBRANDT OLIGODENDROGLIOMA")


heatmap(x=rembrandtEset.culled[which(abs(fitEnet$beta)>0),],
          ColSideColors=c("red","green","blue","yellow")[factor(rembrandtPat$Disease)],
          scale="row",col=redgreen(50),main="Elastic Net")



gene.names <- lapply(rownames(eset)[which(abs(fitEnet$beta)>0)],function(x){
  return(xx[strsplit(x,"_eg")[[1]]])
})

paste(unlist(gene.names),collapse=" ")













#  build a more robust classifier( bootstrapping)
df <- data.frame(value=predictedClasses-1)
rownames(df) <- colnames(eset)

tmp <- c()
try.cv.fit <- c()
tryfit <- c()
for(i in c(1:100)) {
  print(i)
  N <- sample(colnames(eset.culled),559,replace=TRUE)
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
boxplot(yhatbootEnet ~ (predictedClasses.rembrandt-1), ylab="prediction of POSTN expression", xlab="POSTN", main="bootstrap elastic net - molecular features")
stripchart(yhatbootEnet ~ (predictedClasses.rembrandt-1),pch=20, col="royalblue", vertical=TRUE, add=TRUE, cex=.6)




gene.names <- lapply(rownames(eset)[select_features],function(x){
  return(xx[strsplit(x,"_eg")[[1]]])
})

paste(unlist(gene.names),collapse=" ")







