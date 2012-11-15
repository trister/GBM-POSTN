### Andrew Trister
### Sage Bionetworks
### Seattle, WA
### 20120914

### Master script calls for POSTN investigation in gliomas
###


require(synapseClient)

setwd("~/GBM-POSTN/")
synapseLogin()

#get all of the data structures necessary from Synapse
source("./src/populateData.R") 

#now perform a survival analysis on the TCGA based on expression of POSTN
source("./src/POSTN-TCGA-survival.R")
#check the MGMT status
source("./src/MGMT-survival.R")
#and a logistic regression on the same data
source("./src/TCGA-logistic.R")

#and do the survival for all of REMBRANDT
source("./src/POSTN-REMBRANDT-survival.R")



#and now break down each of the grades in REMBRANDT

source("./src/POSTN-grade2and3.R")
source("./src/POSTN-grade2.R")
source("./src/POSTN-grade3.R")
source("./src/POSTN-grade4.R")

#and now all GBMS from REMBRANDT and TCGA
source("./src/POSTN-GBM-all.R")

#and now look at the GEO set from Heidi Phillips GSE4721
source("./src/Phillips-matched-POSTN.R") #matched analysis
source("./src/POSTN-Phillips-survival.R") #survival


#use elastic net to train on TCGA and validate in REMBRANDT
source("./src/enet-POSTN.R")

