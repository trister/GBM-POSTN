### Andrew Trister
### Sage Bionetworks
### Seattle, WA
### 20130809

### Look at the limma analysis for each of the different datasets.
###

require("limma")


### First start with just the TCGA patients


design <- model.matrix(~0+factor(predictedClasses))
# #design <- model.matrix(~age.group)
colnames(design) = c("low", "high")

fit <- lmFit(eset.culled2, design)
fit <- eBayes(fit)
topTable(fit, coef=2)
 
contrast.matrix <- makeContrasts(low-high, levels=design)
fit2 <- contrasts.fit(fit,contrast.matrix)
fit2 <- eBayes(fit2)
topTable(fit2,adjust="BH")#,lfc=0.5)
results <- decideTests(fit2,p.value=0.001,lfc=1)
table(results)
 
gene.names.selected <- lapply(rownames(eset.culled2)[which(results!=0)],function(x){
  return(xx[strsplit(x,"_eg")[[1]]])
})

paste(unlist(gene.names.selected),collapse=" ")

paste(unlist(intersect(gene.names,gene.names.selected)),collapse=" ")






### Now on the combination of all GBM patients


design <- model.matrix(~0+factor(predictedClassesAll))
# #design <- model.matrix(~age.group)
colnames(design) = c("low", "high")

fit <- lmFit(allGBMeset.culled, design)
fit <- eBayes(fit)
topTable(fit, coef=2)

contrast.matrix <- makeContrasts(low-high, levels=design)
fit2 <- contrasts.fit(fit,contrast.matrix)
fit2 <- eBayes(fit2)
topTable(fit2,adjust="BH")#,lfc=0.5)
results <- decideTests(fit2,p.value=0.001,lfc=1)
table(results)

gene.names.selected <- lapply(rownames(eset.culled2)[which(results!=0)],function(x){
  return(xx[strsplit(x,"_eg")[[1]]])
})

paste(unlist(gene.names.selected),collapse=" ")

paste(unlist(intersect(gene.names,gene.names.selected)),collapse=" ")

