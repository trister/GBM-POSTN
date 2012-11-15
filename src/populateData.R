### Andrew Trister
### Sage Bionetworks
### Seattle, WA
### 20120912

### Load all of the required TCGA-GBM data from Synapse
###


require(synapseClient)
require(affy)


####
#### Start with TCGA Coherent Sets

#'syn313583', 'syn372761' - these are the SNM normalized data
dataReturn <- loadEntity('syn313583')
eset <- exprs(dataReturn$objects$eset)

miRNAdata <- loadEntity('syn313603') #SNM normalized miRNA data
miRNAdata <- exprs(miRNAdata$objects$eset)

#proteindata <- loadEntity('syn361966')
#proteindata <- exprs(proteindata$objects$eset)

metadataLoad <- loadEntity('syn673127') #all of the coherent metadata
metadata <- metadataLoad$objects$metadata #extract the R object
metadataAll <- metadata

#exonLoad <- loadEntity('syn313638') #Using Exon expression
#huexdata <- exprs(exonLoad$objects$eset)

mutationLoad <- loadEntity('syn411426') #this gets the level3 calls for mutations from the Freeze 1.1
mutationData <- read.delim(paste(mutationLoad$cacheDir,mutationLoad$files,sep="/"))

#methylationLoad <- loadEntity('syn412284') #frozen 1.1 mutation calls
#methylationData <- read.delim(paste(methylationLoad$cacheDir,methylationLoad$files,sep="/"))))


#load the preprocessed level 3 calls (both 27 and 450)
#methylationLoad <- loadEntity('syn1334844')
methylationLoad <- loadEntity('syn1352976')
methylationData <- methylationLoad$objects$methylation



# get the intersection of the HGU133A files that also have metadata
rows133a <- unlist(lapply(colnames(eset),function(x){ 
  return(grep(x,metadata[,1]))
}))

metadata <- metadata[rows133a,]
rownames(metadata) <- metadata[,1]

inCommon <- intersect(colnames(eset), rownames(metadata))
eset <- eset[,inCommon]
metadata <- metadata[inCommon,]








####
#### and now repeat the entire thing on REMBRANDT across tumor types

rembrandtDataReturn <- loadEntity('syn376920') #REMBRANDT CEL SNM Normalized
rembrandtEset <- rembrandtDataReturn$objects$rembrandt.matrix.matched
rembrandtPat <- rembrandtDataReturn$objects$phenotype.data.matched

rembrandtPat[which(rembrandtPat$Patient.ID=="E09448"),"Grade"] <- " --"



###
### and now repeat the entire thing on recurrent tumors from Heidi Phillips

phillipsDataReturn <- loadEntity('syn360484') #GEO entry GSE4271
phillipsEset <- exprs(phillipsDataReturn$objects$eset)

phillipsMetadata <- read.delim2("./data/GSE4271-GPL96_series_matrix.txt",header=T,stringsAsFactors=F)

rownames(phillipsMetadata) <- phillipsMetadata$X.Sample_title
colnames(phillipsEset) <- rownames(phillipsMetadata)



