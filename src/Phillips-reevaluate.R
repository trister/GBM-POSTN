### Andrew Trister
### Sage Bionetworks
### Seattle, WA
### 20130213

### This will process the data for matched samples from Phillips separately. and then normalize the distributions
### GEO GSE4271


require(synapseClient)
require(affy)
require(snm)
require(mGenomics)
require(ggplot2)
require(samr)
require(glmnet)
require(randomForest)
require(caret)
require(survival)

 synapseLogin()

phillipsCEL <- loadEntity("syn1444929")
CELdir <- file.path(phillipsCEL$cacheDir, phillipsCEL$files)
cdfs <- sapply(CELdir,whatcdf)

#get the metadata

phillipsMetadata <- read.delim2("./data/GSE4271-GPL96_series_matrix.txt",header=T,stringsAsFactors=F)

rownames(phillipsMetadata) <- phillipsMetadata$X.Sample_title




#only do the HG-U133A
files <- names(which(cdfs=="HG-U133A"))


#make the subset of files for the patients with recurrent tumors
recurrent <- grep("recurrent",phillipsMetadata$X.Sample_characteristics_ch1)
files_recurrent <- phillipsMetadata$X.Sample_geo_accession[grep("recurrent",phillipsMetadata$X.Sample_characteristics_ch1)]
files_recurrent <- unlist(lapply(files_recurrent,function(x){
  grep(x,files)
}))

abatch <- ReadAffy(filenames=file.path(files[files_recurrent]))

rawExpr <- log2(pm(abatch))
#rec <- ifelse(1:100 %in% recurrent, "recurrent", "primary")

sr <- snm(rawExpr, int.var=data.frame(array=factor(1:ncol(rawExpr))))
srEval <- fs(sr$norm.dat)
plot(srEval$d)
xyplot(srEval$v[,2] ~ srEval$v[,1])
plot(srEval$v[,1])


myDir <- tempfile()
dir.create(myDir)
for( i in files[files_recurrent] ){
  system(paste("ln -s ", i, " ", myDir, sep=""), ignore.stdout=T, ignore.stderr=T)
}

normDat <- affyWorkflow.SNM(rawDataDir=myDir, int.var=data.frame(array=factor(1:ncol(rawExpr))))

recurrentEset <- exprs(normDat$hgu133a$eset)










#make the subset of files for the patients with matched tumors
files_primary <- unlist(lapply(phillipsMetadata$X.Sample_characteristics_ch1[recurrent],function(x){
    strsplit(x,"Matching sample: ")[[1]][2]
  }))

files_primary <- phillipsMetadata[files_primary,"X.Sample_geo_accession"]
files_primary <- unlist(lapply(files_primary,function(x){
  grep(x,files)
}))

abatch <- ReadAffy(filenames=file.path(files[files_primary]))

rawExpr <- log2(pm(abatch))
#rec <- ifelse(1:100 %in% recurrent, "recurrent", "primary")

sr <- snm(rawExpr, int.var=data.frame(array=factor(1:ncol(rawExpr))))
srEval <- fs(sr$norm.dat)
plot(srEval$d)
xyplot(srEval$v[,2] ~ srEval$v[,1])
plot(srEval$v[,1])


myDir <- tempfile()
dir.create(myDir)
for( i in files[files_primary] ){
  system(paste("ln -s ", i, " ", myDir, sep=""), ignore.stdout=T, ignore.stderr=T)
}

normDat <- affyWorkflow.SNM(rawDataDir=myDir, int.var=data.frame(array=factor(1:ncol(rawExpr))))

primaryEset <- exprs(normDat$hgu133a$eset)



# first rescale the data  to have the same mean and variance 
# Justin's function to rescale the VS to get the same mean/var than the TS
normalize_to_X <- function(mean.x, sd.x, Y){
  m.y <- rowMeans(Y)
  sd.y <- apply(Y, 1, sd)
  Y.adj <- (Y - m.y) * sd.x / sd.y  + mean.x 
  Y.adj[sd.y == 0] <- mean.x[sd.y==0]
  Y.adj
}

tmp <- apply(primaryEset,1,sd)
recurrentEset.scaled <- normalize_to_X(rowMeans(primaryEset),tmp,recurrentEset)


#get rid of the least variant probes
tmp1 <- which(tmp>quantile(tmp,probs=0.2))
primaryEset.scaled <- primaryEset[tmp1,]
recurrentEset.scaled <- recurrentEset.scaled[tmp1,]
phillipsEset.scaled <- phillipsEset[tmp1,]
rm(tmp,tmp1)

colnames(primaryEset.scaled) <- unlist(lapply(colnames(primaryEset.scaled),function(x){
  strsplit(x,".CEL")[[1]][1]
}))

colnames(recurrentEset.scaled) <- unlist(lapply(colnames(recurrentEset.scaled),function(x){
  strsplit(x,".CEL")[[1]][1]
}))



differenceMatched <- recurrentEset.scaled['10631_mt',]/primaryEset.scaled['10631_mt',]
plot(differenceMatched[order(differenceMatched)],main="Relative POSTN expression in recurrent Glioma")











###
# Build an elasticnet model of recurrent versus primary
###

# tempEset <- cbind(primaryEset.scaled,recurrentEset.scaled)
# tempEset <- cbind(primaryEset,recurrentEset)
tempEset <- cbind(phillipsEset[,matchedPrimary],phillipsEset[,recurrent])
#tempClasses <- rep(1:2,each=23)
tempClasses <- c(-1:-23,1:23)

# samfit <- SAM(tempEset,y=tempClasses,resp.type="Two class paired",
#               geneid=rownames(primaryEset),fdr.output=0.001,logged2=TRUE)


 samfit <- SAM(tempEset,y=tempClasses,resp.type="Two class paired",
               geneid=rownames(primaryEset),logged2=TRUE)



genes.interest <- c(samfit$siggenes.table$genes.up[,2],samfit$siggenes.table$genes.lo[,2])
 genes.interest <- samfit$siggenes.table$genes.up[,2]

gene.names <- lapply(genes.interest,function(x){
  return(xx[strsplit(x,"_mt")[[1]]])
})

paste(unlist(gene.names),collapse=" ")



tempClasses <- rep(1:2,each=23)
phillipsFit <- cv.glmnet(x=t(tempEset[genes.interest,]), y=factor(tempClasses), nfolds=10, alpha=.1, family="binomial")
plot(phillipsFit)
phillipsFit <- glmnet(x=t(tempEset[genes.interest,]), y=factor(tempClasses), family="binomial", alpha=.1, lambda=phillipsFit$lambda.min)


temp <- unlist(lapply(genes.interest,function(x){
  sub("_mt","_eg",x)
}))

# yhatEnet <- predict(phillipsFit, t(phillipsEset[genes.interest,]), type="response", s="lambda.min")  #all of Phillips
yhatEnet <- predict(phillipsFit, t(rembrandtEset.culled[temp,]), type="response", s="lambda.min")

riskEnet <- as.vector(ifelse(yhatEnet >= median(yhatEnet), 1, 0))
names(riskEnet) <- rownames(yhatEnet)


plot(survfit(tmpSurv.rembrandt ~ riskEnet), main = "GBM K-M Plot By POSTN Expression (TCGA)")
survdiff(tmpSurv.rembrandt~riskEnet)


#boxplots:

plot(rembrandtEset['10631_mt',which(rembrandtPat$Grade!=" --" & rembrandtPat$Grade!=" I")]~factor(rembrandtPat$Grade[which(rembrandtPat$Grade!=" --" & rembrandtPat$Grade!=" I")]),
     xlab="Glioma Grade", ylab="Expression",main="Expression of POSTN")

aov.rembrandt <- aov(rembrandtEset['10631_mt',which(rembrandtPat$Grade!=" --" &rembrandtPat$Grade!=" I")]~factor(rembrandtPat$Grade[which(rembrandtPat$Grade!=" --" & rembrandtPat$Grade!=" I")]))
summary(aov.rembrandt)
print(model.tables(aov.rembrandt,"means"),digits=3)
