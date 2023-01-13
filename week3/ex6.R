library(devtools)
library(Biobase)
library(limma)
library(edge)
library(genefilter)

con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
fdata = fData(mp)

# 
# con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
# load(file=con)
# close(con)
# mp = montpick.eset
# pdata=pData(mp)
# edata=as.data.frame(exprs(mp))
# fdata = fData(mp)


edata = log2(as.matrix(edata) + 1)

## ------------------------------------------------------------------------
tstats_obj = genefilter::rowttests(edata,pdata$population)

## ------------------------------------------------------------------------
fstats_obj = rowFtests(edata,as.factor(pdata$population))

par(mfrow=c(1,2))
hist(tstats_obj$p.value, col=2)
hist(fstats_obj$p.value, col=2)

par(mfrow=c(1,2))
hist(tstats_obj$statistic, col=2)
hist(fstats_obj$statistic, col=2)

