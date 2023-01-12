
## ----load_hidden, echo=FALSE, results="hide", warning=FALSE--------------
suppressPackageStartupMessages({
  library(devtools)
  library(Biobase)
})

## ----load----------------------------------------------------------------
library(devtools)
library(Biobase)

par(pch = 19)

con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
fdata = fData(mp)

edata = log2(edata + 1)
edata_centered = edata - rowMeans(edata)
svd1 = svd(edata_centered)

set.seed(333)
k1 = kmeans(t(edata), centers=2)
cor.test(svd1$v[,1], k1$cluster)
# matplot(t(k1$centers), col=1:2, type="l",lwd=3)
# table(k1$cluster)
