### Andrew Trister
### Sage Bionetworks
### Seattle, WA
### 20120912

### Use the populated data to look at survival with POSTN in Phillips GEO GSE4271

require(ggplot2)
require(survival)
require(mixtools)
source("./src/ggkm.R")

#build a survival object

phillipsSurv <- unlist(lapply(phillipsMetadata[matchedPrimary,'X.Sample_characteristics_ch1'],function(x){
  temp <- strsplit(x,";")[[1]][5]
  temp <- strsplit(temp,": ")[[1]][2]
  if(temp=="no") 
    return(1)
  else return(0)
}))


phillipsSurvTime <- unlist(lapply(phillipsMetadata[matchedPrimary,'X.Sample_characteristics_ch1'],function(x){
  temp <- strsplit(x,";")[[1]][4]
  return(strsplit(temp,": ")[[1]][2])
}))
phillipsSurvTime <- as.numeric(phillipsSurvTime)


tmpSurv.phillips <- Surv(phillipsSurvTime,phillipsSurv)

temp <- quantile(differenceMatched)
phillipsGroup <- unlist(lapply(differenceMatched,function(x){
  if(x<0.8) 1
  else if(x>1.2) 2
  else 3
}))

##let's look at survival based on class for POSTN
ggkm(survfit(tmpSurv.phillips ~ phillipsGroup), ystratalabs = (c("POSTN Lower", "POSTN Higher", "POSTN Equal")), 
     timeby = 52,
     main = "Change in POSTN Expression with Recurrence (Phillips)")






