library(devtools)
library(Biobase)
library(dendextend)

BiocManager::install("rafalib")
library(rafalib)

rm(list=ls())

con = url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
fdata = fData(mp)
study_id = as.integer(pdata$study == "Montgomery")
set.seed(1235)
log2_edata = log2(edata + 1)
k1 = kmeans(t(log2_edata), centers=2)
matplot(t(k1$centers), col=1:2, type="l",lwd=3)
table(k1$cluster)
# Assume cluster 1 is study_id_1 and get overlap
# Map cluster to 0/1 same as study_id
kmeans_clusters = as.integer(k1$cluster == 1)
frac_correct = sum(kmeans_clusters == study_id) / length(study_id)
frac_correct
# frac_correct = 84%
# cuttree
d1 = dist(t(log2_edata))
h1 = hclust(d1)
c1 = cutree(h1, k = 2)
table(c1)
# Assume cluster 1 is study 1
hclust_clusters = as.integer(c1 == 1)
frac_correct = sum(hclust_clusters == study_id) / length(study_id)
frac_correct
# frac_correct = 45%
