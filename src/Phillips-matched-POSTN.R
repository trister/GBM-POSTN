### Andrew Trister
### Sage Bionetworks
### Seattle, WA
### 20120912

### Use the populated data to look at Phillips data for matched samples at recurrence
### for differences in POSTN expression


require(ggplot2)
require(survival)
require(mixtools)
require(reshape2)
source("./src/ggkm.R")


hist(phillipsEset['10631_mt',])

boxplot(phillipsEset['10631_mt',grep("recurrent",phillipsMetadata$X.Sample_characteristics_ch1)],
        phillipsEset['10631_mt',-grep("recurrent",phillipsMetadata$X.Sample_characteristics_ch1)],
        main="POSTN in Recurrent vs. Primary Glioma (Phillips)")

axis(1,at=c(1,2),labels=c("Recurrent","Primary"))

t.test(phillipsEset['10631_mt',grep("recurrent",phillipsMetadata$X.Sample_characteristics_ch1)],
       phillipsEset['10631_mt',-grep("recurrent",phillipsMetadata$X.Sample_characteristics_ch1)],
       main="POSTN in Recurrent vs. Primary Glioma (Phillips)")



boxplot(phillipsEset['10631_mt',intersect(grep("recurrent",phillipsMetadata$X.Sample_characteristics_ch1),grep("WHO grade: IV",phillipsMetadata$X.Sample_characteristics_ch1))],
        phillipsEset['10631_mt',-c(grep("recurrent",phillipsMetadata$X.Sample_characteristics_ch1),grep("WHO grade: III",phillipsMetadata$X.Sample_characteristics_ch1))],
        main="POSTN in Recurrent vs. Primary GBM (Phillips)")
axis(1,at=c(1,2),labels=c("Recurrent","Primary"))

t.test(phillipsEset['10631_mt',intersect(grep("recurrent",phillipsMetadata$X.Sample_characteristics_ch1),grep("WHO grade: IV",phillipsMetadata$X.Sample_characteristics_ch1))],
       phillipsEset['10631_mt',-c(grep("recurrent",phillipsMetadata$X.Sample_characteristics_ch1),grep("WHO grade: III",phillipsMetadata$X.Sample_characteristics_ch1))],
       main="POSTN in Recurrent vs. Primary GBM (Phillips)")







#look at the matched samples
recurrent <- grep("recurrent",phillipsMetadata$X.Sample_characteristics_ch1)

matchedPrimary <- unlist(lapply(phillipsMetadata$X.Sample_characteristics_ch1[recurrent],function(x){
  strsplit(x,"Matching sample: ")[[1]][2]
}))

boxplot(phillipsEset['10631_mt',recurrent],phillipsEset['10631_mt',matchedPrimary],main="POSTN in Matched Primary and Recurrent Gliomas (Phillips)")
axis(1,at=c(1,2),labels=c("Recurrent","Primary"))

t.test(phillipsEset['10631_mt',recurrent],phillipsEset['10631_mt',matchedPrimary])



differenceMatched <- phillipsEset['10631_mt',recurrent]/phillipsEset['10631_mt',matchedPrimary]

plot(differenceMatched[order(differenceMatched)],main="Relative POSTN expression in recurrent Glioma")




recurrentGBM <- intersect(recurrent,grep("WHO grade: IV",phillipsMetadata$X.Sample_characteristics_ch1))
matchedPrimaryGBM <- intersect(matchedPrimary,phillipsMetadata$X.Sample_title[grep("WHO grade: IV",phillipsMetadata$X.Sample_characteristics_ch1)])

recurrentGBM <- unlist(lapply(phillipsMetadata[matchedPrimaryGBM,'X.Sample_characteristics_ch1'],function(x){
  strsplit(x,"Matching sample: ")[[1]][2]
}))

#boxplot(phillipsEset['10631_mt',recurrentGBM],phillipsEset['10631_mt',matchedPrimaryGBM],main="POSTN in Matched Primary GIII and Recurrent GBM (Phillips)")
boxplot(phillipsEset['10631_mt',recurrentGBM],phillipsEset['10631_mt',matchedPrimaryGBM],main="POSTN in Matched Primary and Recurrent GBM (Phillips)")
axis(1,at=c(1,2),labels=c("Recurrent","Primary"))

t.test(phillipsEset['10631_mt',recurrentGBM],phillipsEset['10631_mt',matchedPrimaryGBM])


# now make a nice ggplot
group <- rep(c("primary","recurrence"),each=length(recurrent))
df <- data.frame(id=seq_along(group),group,phillipsEset['10631_mt',matchedPrimary],phillipsEset['10631_mt',recurrent])
dfm <- melt(df,id.var=c("id","group"))

ggplot(dfm,aes(variable,2^value,group=id)) + geom_path(alpha=0.5) +coord_trans(y="log2")
