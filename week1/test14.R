library(devtools)
library(Biobase)
library(dendextend)
library(rafalib)

rm(list=ls())

con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
fdata = fData(mp)

#par(mfrow=c(2,1))

# By default calculates the distance between rows
dist1 = dist(t(edata))
hclust1=hclust(dist1)
myplclust(hclust1,labels = pdata$sample.id, lab.col = as.numeric(pdata$study),hang=0.1)

## ------------------------------------------------------------------------
low_genes=rowMeans(edata) <100
filt_data=filter(as.data.frame(edata),!low_genes)
dist2=dist(t(filt_data))
hclust2=hclust(dist2)
myplclust(hclust2,labels = pdata$sample.id, lab.col = as.numeric(pdata$study),hang=0.1)

## ------------------------------------------------------------------------
dist3=dist(t(log2(edata+1)))
hclust3=hclust(dist3)
myplclust(hclust3,labels = pdata$sample.id, lab.col = as.numeric(pdata$study) ,hang=0.1)

## ----session_info--------------------------------------------------------
devtools::session_info()