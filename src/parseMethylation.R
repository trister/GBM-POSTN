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

#myData <- ExpressionData(list(name="New combined 27 and 450", parentId = properties(myStudy)$id))
#myData <- createEntity(myData)


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

#find the subjects that were repeated
columnsTemp <- substring(columns,1,12)

valTemp <- lapply(unique(columnsTemp[which(duplicated(columnsTemp))]),function(x){
  temp <- methylationTemp[,grep(x, columns)]
  return(apply(temp,1,function(y){mean(y,na.rm=T)}))
})

for(i in 1:length(unique(columnsTemp[which(duplicated(columnsTemp))]))) {
  temp <- grep(unique(columnsTemp[which(duplicated(columnsTemp))])[i],columns)
  methylationTemp <- methylationTemp[,-temp]
  columns <- columns[-temp]
  methylationTemp <- cbind(methylationTemp,valTemp[[i]])
  columns <- c(columns,unique(columnsTemp[which(duplicated(columnsTemp))])[i])
  
}

colnames(methylationTemp) <- substring(columns,1,12)
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



#find the subjects that were repeated
columnsTemp <- substring(columns,1,12)

valTemp <- lapply(unique(columnsTemp[which(duplicated(columnsTemp))]),function(x){
  temp <- methylationTemp[,grep(x, columns)]
  return(apply(temp,1,function(y){mean(y,na.rm=T)}))
})

for(i in 1:length(unique(columnsTemp[which(duplicated(columnsTemp))]))) {
  temp <- grep(unique(columnsTemp[which(duplicated(columnsTemp))])[i],columns)
  methylationTemp <- methylationTemp[,-temp]
  columns <- columns[-temp]
  methylationTemp <- cbind(methylationTemp,valTemp[[i]])
  columns <- c(columns,unique(columnsTemp[which(duplicated(columnsTemp))])[i])
  
}


colnames(methylationTemp) <- substring(columns,1,12)
rownames(methylationTemp) <- combination[[1]]$probe.name

#now check to make sure that there are no repeats in the original set from 27k
temp <- intersect(colnames(methylationTemp),colnames(methylation))
methylationTemp <- methylationTemp[,!colnames(methylationTemp) %in% temp]



methylation <- cbind(methylation,methylationTemp[intersectProbes,])


#now push the entire thing up to Synapse
#myData <- getEntity('syn1334844')

myData <- ExpressionData(list(name="New combined 27 and 450", parentId = myData$properties$parentId))
myData <- createEntity(myData)


myData <- addObject(myData, methylation)
myData <- storeEntity(myData) #syn1352976



