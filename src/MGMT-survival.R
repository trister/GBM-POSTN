### Andrew Trister
### Sage Bionetworks
### Seattle, WA
### 20120917

### Make a call of whether a patient is MGMT methylated based on the Bady Paper (Acta Neurologica 2012)
###

require(survival)
require(ggplot2)
source("./src/ggkm.R")


#first grab the two probes that are necessary for the analysis
culledMethylation <- methylationData[c('cg12434587','cg12981137'),unique(colnames(methylationData))]

MGMT <- apply(culledMethylation,2,function(x){
  M1 <- log2(x[1]/(1-x[1]))
  M2 <- log2(x[2]/(1-x[2]))
  alpha <- 4.3215+0.5271*M1+0.9265*M2
  if (exp(alpha)/(1+exp(alpha)) > 0.358)
    return(1)
  else return(0)
})

culledMethylation <- rbind(culledMethylation,MGMT)

#now find the subset of the patients that have methylation calls
MGMTculled <- culledMethylation[,which(colnames(culledMethylation) %in% metadata$bcr_patient_barcode)]
metaMeth <- metadata[which(unique(as.character(metadata$bcr_patient_barcode)) %in% colnames(culledMethylation)),]

#now make a survival object
metaMeth$vital_status[metaMeth$vital_status=="[Not Available]"] <- NA
gbmPatMeth <- c()
gbmPatMeth$survTime <- as.numeric(metaMeth$days_to_death)
gbmPatMeth$surv <- ifelse(metaMeth$vital_status=="DECEASED", 1,0)
gbmPatMeth$survTime[ which(metaMeth$vital_status=="LIVING")] <- as.numeric(metaMeth$days_to_last_followup)[ which(metadata$vital_status=="LIVING")]

tmpSurv <- Surv(gbmPatMeth$survTime,gbmPatMeth$surv)


##let's look at survival based on class for POSTN
ggkm(survfit(tmpSurv ~ MGMTculled['MGMT',]), ystratalabs = (c("MGMT Normal", "MGMT Methylated")), 
     timeby = 365,
     main = "GBM K-M Plot By MGMT Methylation (TCGA)")
