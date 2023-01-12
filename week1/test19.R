library(devtools)
library(Biobase)
library(dendextend)

BiocManager::install("rafalib")
library(rafalib)
  
rm(list=ls())

tropical=  c('darkorange', 'dodgerblue', 'hotpink', 'limegreen', 'yellow')
palette(tropical)

```
{r question9, echo=FALSE}  
```
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
fdata = fData(mp)

# Hierarchical Cluster1
d1 = dist(t(edata))
h1 = hclust(d1)
plot(h1)
myplclust(h1, lab.col = as.numeric(pdata$study))


# Cluster2
d2 = dist(t(edata[rowMeans(edata) >= 100,]))
h2 = hclust(d2)
plot(h2)
myplclust(h2, lab.col = as.numeric(pdata$study))
# Cluster3
d3 = dist(t(log2(edata + 1)))
h3 = hclust(d3)
plot(h3)
myplclust(h3, lab.col = as.numeric(pdata$study))

## ----session_info--------------------------------------------------------
devtools::session_info()