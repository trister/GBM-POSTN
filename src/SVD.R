## Andrew Trister
## Sage Bionetworks
## Seattle, WA
## 20121119

## FIRST: GENERATING UN-NORMALIZED, UN-BACKGROUND CORRECTED DATA
## We'll use the 'rma()' function with the normalize and background arguments set to false.

## REQUIRED LIBRARIES
require(Biobase)
require(affy)
require(snm)
require(ggplot2)
require(synapseClient)


# Justin Guinney's function to rescale the validation data to get the same mean/var than the training set
normalize_to_X <- function(mean.x, sd.x, Y){
  m.y <- rowMeans(Y)
  sd.y <- apply(Y, 1, sd)
  Y.adj <- (Y - m.y) * sd.x / sd.y  + mean.x 
  Y.adj[sd.y == 0] <- mean.x[sd.y==0]
  Y.adj
}


#cull the REMBRANDT U133Aplus2 to U133A genesets
temp <- unlist(lapply(rownames(rembrandtEset),function(x){
  return(gsub("_mt","_eg",x))
}))

rembrandtEset.culled <- rembrandtEset
rownames(rembrandtEset.culled) <- temp

temp <- intersect(temp,rownames(eset))

eset.culled <- eset[temp,]
rembrandtEset.culled <- rembrandtEset.culled[temp,]


fullRawMat <- cbind(eset.culled,rembrandtEset.culled)

studyIndicator <- c(rep('TCGA', ncol(eset.culled)),
                    rep('REMBRANDT', ncol(rembrandtEset.culled)))

#Plot the principal components
svdObj <- svd(fullRawMat)
plot(svdObj$v[,1],svdObj$v[,2],
     col=c("royalblue","orange","aquamarine4","brown2")[as.factor(studyIndicator)],
     bg=c("royalblue","orange","aquamarine4","brown2")[as.factor(studyIndicator)],
     pch=c(8,24,23,19)[as.factor(studyIndicator)], font=2,
     cex=0.9, xlab="Principal component 1", ylab="Principal component 2")


#plot the scaled SNM SVD




rembrandtEset.scaled <- normalize_to_X(rowMeans(eset.culled), apply(eset.culled, 1, sd), rembrandtEset.culled)

# describe the latent structure
s <- svd(cbind(rembrandtEset.scaled,eset.culled))
plot(s$v[,1],s$v[,2],
     col=c("royalblue","orange","aquamarine4","brown2")[as.factor(studyIndicator)],
     bg=c("royalblue","orange","aquamarine4","brown2")[as.factor(studyIndicator)],
     pch=c(8,24,23,19)[as.factor(studyIndicator)], font=2,
     cex=0.9, xlab="Principal component 1", ylab="Principal component 2")


