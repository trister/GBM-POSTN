### Andrew Trister
### Sage Bionetworks
### Seattle, WA
### 20121018

### This will process the data for the Phillips to look at changes in matched samples
### GEO GSE4271


require(synapseClient)
require(affy)
require(snm)
require(mGenomics)
require(ggplot2)

# synapseLogin()

phillipsCEL <- loadEntity("syn1444929")
CELdir <- file.path(phillipsCEL$cacheDir, phillipsCEL$files)
cdfs <- sapply(CELdir,whatcdf)

#get the metadata

phillipsMetadata <- read.delim2("./data/GSE4271-GPL96_series_matrix.txt",header=T,stringsAsFactors=F)

rownames(phillipsMetadata) <- phillipsMetadata$X.Sample_title




#only do the HG-U133A
files <- names(which(cdfs=="HG-U133A"))
abatch <- ReadAffy(filenames=file.path(files))

rawExpr <- log2(pm(abatch))
rec <- ifelse(1:100 %in% recurrent, "recurrent", "primary")

sr <- snm(rawExpr, int.var=data.frame(array=factor(1:ncol(rawExpr))))
srEval <- fs(sr$norm.dat)
plot(srEval$d)
xyplot(srEval$v[,2] ~ srEval$v[,1], groups=rec)
plot(srEval$v[,1])


myDir <- tempfile()
dir.create(myDir)
for( i in files ){
  system(paste("ln -s ", i, " ", myDir, sep=""), ignore.stdout=T, ignore.stderr=T)
}

normDat <- affyWorkflow.SNM(rawDataDir=myDir, int.var=data.frame(array=factor(1:ncol(rawExpr))))

postn <- exprs(normDat$hgu133a$eset)["10631_mt",]

recNames <- phillipsMetadata[grep("recurrent",phillipsMetadata$X.Sample_characteristics_ch1), "X.Sample_supplementary_file"]
primNames <- unlist(lapply(phillipsMetadata$X.Sample_characteristics_ch1[match(recNames, phillipsMetadata$"X.Sample_supplementary_file")],function(x){
  strsplit(x,"Matching sample: ")[[1]][2]
}))
primNames <- phillipsMetadata[primNames, "X.Sample_supplementary_file"]

recNames <- sub(".gz", "", basename(recNames), fixed=T)
primNames <- sub(".gz", "", basename(primNames), fixed=T)

postnFC <- postn[recNames] - postn[primNames]

plot(postn[primNames],postn[recNames])
abline(0,1)
abline(-1,1)
abline(1,1)

grp <- ifelse(names(postn) %in% recNames, "rec","prim")
boxplot(postn~grp)
















# ####
# # Old stuff
# 
# int <- intensity(abatch)
# gene2row.cdf <- getGene2Row(NULL,annotation(abatch))
# 
# pms <- int[unlist(gene2row.cdf),]
# pms[pms<=1] <- 1
# data <- log2(pms)
# 
# colnames(data) <- rownames(phillipsMetadata)
# 
# 
# recurrent <- phillipsMetadata[grep("recurrent",phillipsMetadata$X.Sample_characteristics_ch1), "X.Sample_supplementary_file"]
# 
# matchedPrimary <- unlist(lapply(phillipsMetadata$X.Sample_characteristics_ch1[recurrent],function(x){
#   strsplit(x,"Matching sample: ")[[1]][2]
# }))
# 
# 
# 
# group <- rep(c("primary","recurrence"),each=length(recurrent))
# 
# df <- data.frame(id=seq_along(group),group,data[10631,matchedPrimary],data[10631,recurrent])
# dfm <- melt(df,id.var=c("id","group"))
# 
# ggplot(dfm,aes(variable,2^value,group=id)) + geom_path(alpha=0.5) +coord_trans(y="log2")
# 
# 
# data.culled <- data[,c(recurrent,matchedPrimary)]
# 
# 
# 
# 
# 
# phillipsMetadata$sampleType <- sapply(phillipsMetadata$X.Sample_characteristics_ch1,function(x) {
#   if(grepl("primary",x)) return(1)
#   else return(0)
# })
# 
# 
# int.var <- data.frame(array=factor(1:ncol(data)))
# 
# scanDate <- abatch@protocolData@data$ScanDate
# scanDate <- sapply(strsplit(scanDate, " "),"[[",1)
# scanDate <- as.Date(scanDate,format="%m/%d/%y")
# 
# batch <- as.factor(scanDate)
# adj.var <- model.matrix(~batch)
# bio.var <- model.matrix(~phillipsMetadata$sampleType)
# 
# snm.fit <- snm(data,bio.var=bio.var,adj.var=adj.var,int.var=int.var,rm.adj=T,diagnose=F)
# 
# svdObj <- fast.svd(snm.fit$norm.dat)
# pcDF <- data.frame(svdObj$v[,1:2])
# colnames(pcDF) <- c("PrinComp1", "PrinComp2")
# 
# pcPlot <- ggplot(pcDF,aes(PrinComp1,PrinComp2)) + 
#   geom_point(aes(color=batch),size=20) + scale_size(guide="none")
# 
# pcPlot