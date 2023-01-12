## ----load_hidden, echo=FALSE, results="hide", warning=FALSE--------------
suppressPackageStartupMessages({
  library(devtools)
  library(Biobase)
})

## ----load----------------------------------------------------------------
library(devtools)
library(Biobase)

par(pch = 19)
par(mfrow=c(1,3))

con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
fdata = fData(mp)
ls()

## ------------------------------------------------------------------------
svd1 = svd(edata)
#plot(svd1$d^2/sum(svd1$d^2),ylab="Percent Variance Explained",col=2)
svd1[1]

## ------------------------------------------------------------------------
edata = log2(edata + 1)
svd2 = svd(edata)
#plot(svd2$d^2/sum(svd2$d^2),ylab="Percent Variance Explained",col=2)

## ------------------------------------------------------------------------
edata_centered = edata - rowMeans(edata)
svd3 = svd(edata_centered)
#plot(svd3$d^2/sum(svd3$d^2),ylab="Percent Variance Explained",col=2)
