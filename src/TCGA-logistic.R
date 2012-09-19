### Andrew Trister
### Sage Bionetworks
### Seattle, WA
### 20120914

### Populate methylation and mutation and perform a logistic regression.


require(ggplot2)
require(survival)
require(mixtools)


####
### Logistic regression of POSTN high and low versus other features in TCGA
####

table(predictedClasses)
metadata$karnofsky_performance_score[metadata$karnofsky_performance_score=="[Pending]"] <- NA
metadata$karnofsky_performance_score <- as.numeric(metadata$karnofsky_performance_score)



#put in the MGMT status
metadata$MGMT <- rep(NA,length(metadata$bcr_patient_barcode))

for (i in 1:length(metadata$bcr_patient_barcode)) {
  if(length(grep(metadata$bcr_patient_barcode[i],colnames(MGMTculled)))!=0) {
    metadata$MGMT[i] <- MGMTculled['MGMT',grep(metadata$bcr_patient_barcode[i],colnames(MGMTculled))]
  }
}

#Find the IDH1 mutations
length(unique(mutationData$Tumor_Sample_Barcode))



IDH1 <- substring(mutationData$Tumor_Sample_Barcode[grep("IDH1", mutationData$Hugo_Symbol)],1,15)
PTEN <- substring(mutationData$Tumor_Sample_Barcode[grep("PTEN", mutationData$Hugo_Symbol)],1,15)
P53 <- substring(mutationData$Tumor_Sample_Barcode[grep("P53", mutationData$Hugo_Symbol)],1,15)
PI3K <- substring(mutationData$Tumor_Sample_Barcode[grep("PI3K", mutationData$Hugo_Symbol)],1,15)
RB <- substring(mutationData$Tumor_Sample_Barcode[grep("RB", mutationData$Hugo_Symbol)],1,15)
AKT <- substring(mutationData$Tumor_Sample_Barcode[grep("AKT", mutationData$Hugo_Symbol)],1,15)
RTK <- substring(mutationData$Tumor_Sample_Barcode[grep("RTK", mutationData$Hugo_Symbol)],1,15)



temp <- substring(colnames(eset),1,15)
#metadata$IDH1 <- rep(0,length(metadata$bcr_patient_barcode))
metadata$RTK <- metadata$PTEN <-  metadata$P53 <- metadata$PI3K <- metadata$RB <- metadata$AKT <- metadata$IDH1 <- rep(NA,length(metadata$bcr_patient_barcode))

allMutant <- which(temp %in% substring(mutationData$Tumor_Sample_Barcode,1,15))

metadata$IDH1[allMutant] <- metadata$PTEN[allMutant] <- metadata$P53[allMutant] <- metadata$PI3K[allMutant] <- 0
metadata$RB[allMutant] <- metadata$AKT[allMutant] <- metadata$RTK[allMutant] <- 0

metadata$IDH1[which(temp %in% IDH1)] <- 1
metadata$PTEN[which(temp %in% PTEN)] <- 1
metadata$P53[which(temp %in% P53)] <- 1
metadata$PI3K[which(temp %in% PI3K)] <- 1
metadata$RB[which(temp %in% RB)] <- 1
metadata$AKT[which(temp %in% AKT)] <- 1
metadata$RTK[which(temp %in% RTK)] <- 1

as.factor(metadata$AKT)

  
    
#logistic regression time
mylogit <- glm(predictedClasses~metadata$karnofsky_performance_score+metadata$age_at_initial_pathologic_diagnosis+
  as.factor(metadata$gender)+as.factor(metadata$IDH1)+as.factor(metadata$P53)+
  as.factor(metadata$PTEN)+as.factor(metadata$AKT)+as.factor(metadata$RTK)+
  as.factor(metadata$RB) +as.factor(metadata$MGMT),na.action=na.exclude)

summary(mylogit)


exp(coef(mylogit))
exp(cbind(OR=coef(mylogit),confint(mylogit)))



