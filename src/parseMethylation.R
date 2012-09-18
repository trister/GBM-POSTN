### Andrew Trister
### Sage Bionetworks
### Seattle, WA
### 20120917

### Prepare all of the methylation data to upload to Synapse
###

require(synapseClient)

synapseLogin()

#first make a new project
#myProject <- Project(list(name="TCGA GBM Methylation"))
#myProject <- createEntity(myProject)

#myStudy <- Study(list(name="Level 3 Methylation", parentId=propertyValue(myProject,"id")))
#myStudy <- createEntity(myStudy)

#myData <- ExpressionData(list(name="Combined 27 and 450", parentId = properties(myStudy)$id))
#myData <- createEntity(myData)

#myData <- loadEntity('syn1334844')


# start by running a loop over the Methylation27 files
mainDir <- "~/DNA_Methylation/JHU_USC__HumanMethylation27/Level_3/"
allFiles <- list.files(mainDir)

combination <- lapply(allFiles,function(x){
  return(read.delim(paste(mainDir,x,sep="/"), stringsAsFactors=F))
})

methylationTemp <- c()
columns <- c()
for(i in 1:length(combination)){
  methylationTemp <- cbind(methylationTemp, as.numeric(combination[[i]]$beta.value))
  columns <- c(columns,combination[[i]]$barcode[1])
}


colnames(methylationTemp) <- columns
rownames(methylationTemp) <- combination[[1]]$probe.name


methylation <- methylationTemp

#now repeat for the Methylation450 files
mainDir <- "~/DNA_Methylation/JHU_USC__HumanMethylation450/Level_3/"
allFiles <- list.files(mainDir)

combination <- lapply(allFiles,function(x){
  return(read.delim(paste(mainDir,x,sep="/"), stringsAsFactors=F))
})

intersectProbes <- intersect(rownames(methylation), combination[[1]]$probe.name)
methylation <- methylation[intersectProbes,]

methylationTemp <- c()
columns <- c()

for(i in 1:length(combination)){
  methylationTemp <- cbind(methylationTemp, as.numeric(combination[[i]]$beta.value))
  columns <- c(columns,combination[[i]]$barcode[1])
}


colnames(methylationTemp) <- substring(columns,1,12)
rownames(methylationTemp) <- combination[[1]]$probe.name

intersect(names(methylationTemp),names(methylation))

methylation <- cbind(methylation,methylationTemp[intersectProbes,])


#now push the entire thing up to Synapse
myData <- addObject(myData, methylation)
myData <- storeEntity(myData)

#myData <- loadEntity('syn1334844')