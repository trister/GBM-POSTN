### Andrew Trister
### Sage Bionetworks
### Seattle, WA
### 20130210

### Use the elasticnet model built in previous steps to verify genes in the matched set


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




# first rescale the data  to have the same mean and variance 
# Justin's function to rescale the VS to get the same mean/var than the TS
normalize_to_X <- function(mean.x, sd.x, Y){
  m.y <- rowMeans(Y)
  sd.y <- apply(Y, 1, sd)
  Y.adj <- (Y - m.y) * sd.x / sd.y  + mean.x 
  Y.adj[sd.y == 0] <- mean.x[sd.y==0]
  Y.adj
}

phillipsEset.culled <- phillipsEset[-grep("10631_mt",rownames(phillipsEset)),]

tmp <- apply(eset.culled,1,sd)
phillipsEset.scaled <- normalize_to_X(rowMeans(eset.culled),tmp,phillipsEset.culled)



#get rid of the least variant probes
tmp1 <- which(tmp>quantile(tmp,probs=0.2))
eset.culled2 <- eset.culled[tmp1,]
phillipsEset.scaled <- phillipsEset.scaled[tmp1,]
rm(tmp,tmp1)



cv.fit <- cv.glmnet(x=t(eset.culled2), y=factor(predictedClasses), nfolds=10, alpha=.1, family="binomial")
plot(cv.fit)
fitEnet <- glmnet(x=t(eset.culled2), y=factor(predictedClasses), family="binomial", alpha=.1, lambda=cv.fit$lambda.min)







###############################################################################
#  Look now at all samples in Phillips                                #
###############################################################################


yhatEnet <- predict(fitEnet, t(phillipsEset.scaled), type="response", s="lambda.min")  #all of Phillips


#boxplot(yhatEnet~grepl("recurrent",phillipsMetadata$X.Sample_characteristics_ch1))
boxplot(yhatEnet[recurrent],yhatEnet[which(colnames(phillipsEset) %in% matchedPrimary)]) 
t.test(yhatEnet[recurrent],yhatEnet[which(colnames(phillipsEset) %in% matchedPrimary)]) 






###############################################################################
#  Look now at the matched samples in Phillips                                #
###############################################################################


#  yhatEnet <- predict(fitEnet, t(phillipsEset.scaled[,c(recurrent,which(colnames(phillipsEset) %in% matchedPrimary))]), type="response", s="lambda.min")  #all of Phillips

yhatEnet <- predict(fitEnet, t(phillipsEset.scaled), type="response", s="lambda.min")  #all of Phillips


#boxplot(yhatEnet~grepl("recurrent",phillipsMetadata$X.Sample_characteristics_ch1))
boxplot(yhatEnet[recurrent],yhatEnet[which(colnames(phillipsEset) %in% matchedPrimary)]) 
t.test(yhatEnet[recurrent],yhatEnet[which(colnames(phillipsEset) %in% matchedPrimary)]) 




###############################################################################
#  What about only the grade IIIs?                                            #
###############################################################################

tempPrimary <- intersect(which(rownames(phillipsMetadata) %in% matchedPrimary),grep("WHO grade: III",phillipsMetadata$X.Sample_characteristics_ch1))
tempPrimary <- grep("WHO grade: III;",phillipsMetadata$X.Sample_characteristics_ch1[rownames(phillipsMetadata) %in% matchedPrimary])
tempRecurrent <- unlist(lapply(phillipsMetadata$X.Sample_characteristics_ch1[tempPrimary],function(x){
  strsplit(x,"Matching sample: ")[[1]][2]
}))
tempRecurrent <- which(rownames(phillipsMetadata) %in% tempRecurrent)

boxplot(phillipsEset['10631_mt',tempRecurrent],phillipsEset['10631_mt',tempPrimary],main="POSTN in Matched Primary and Recurrent Glioma (Phillips)")
axis(1,at=c(1,2),labels=c("Recurrent","Primary"))

t.test(phillipsEset['10631_mt',tempRecurrent],phillipsEset['10631_mt',tempPrimary])




yhatEnet <- predict(fitEnet, t(phillipsEset.scaled), type="response", s="lambda.min")  #all of Phillips
# yhatEnet <- predict(fitEnet, t(phillipsEset.scaled[,c(tempRecurrent,tempPrimary)]), type="response", s="lambda.min")  #all of Phillips


#boxplot(yhatEnet~grepl("recurrent",phillipsMetadata$X.Sample_characteristics_ch1))
boxplot(yhatEnet[tempRecurrent],yhatEnet[tempPrimary]) 
t.test(yhatEnet[tempRecurrent],yhatEnet[tempPrimary]) 




###############################################################################
#  What about only the grade IVs?                                            #
###############################################################################

tempPrimary <- intersect(which(rownames(phillipsMetadata) %in% matchedPrimary),grep("WHO grade: IV",phillipsMetadata$X.Sample_characteristics_ch1))
tempPrimary <- grep("WHO grade: IV",phillipsMetadata$X.Sample_characteristics_ch1[rownames(phillipsMetadata) %in% matchedPrimary])
tempRecurrent <- unlist(lapply(phillipsMetadata$X.Sample_characteristics_ch1[tempPrimary],function(x){
  strsplit(x,"Matching sample: ")[[1]][2]
}))
tempRecurrent <- which(rownames(phillipsMetadata) %in% tempRecurrent)

boxplot(phillipsEset['10631_mt',tempRecurrent],phillipsEset['10631_mt',tempPrimary],main="POSTN in Matched Primary and Recurrent Glioma (Phillips)")
axis(1,at=c(1,2),labels=c("Recurrent","Primary"))

t.test(phillipsEset['10631_mt',tempRecurrent],phillipsEset['10631_mt',tempPrimary])




yhatEnet <- predict(fitEnet, t(phillipsEset.scaled), type="response", s="lambda.min")  #all of Phillips
# yhatEnet <- predict(fitEnet, t(phillipsEset.scaled[,c(tempRecurrent,tempPrimary)]), type="response", s="lambda.min")  #all of Phillips


#boxplot(yhatEnet~grepl("recurrent",phillipsMetadata$X.Sample_characteristics_ch1))
boxplot(yhatEnet[tempRecurrent],yhatEnet[tempPrimary]) 
t.test(yhatEnet[tempRecurrent],yhatEnet[tempPrimary]) 

