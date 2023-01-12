rm(list=ls())
library(devtools)
library(Biobase)
library(snpStats)
library(broom)
library(MASS)
library(DESeq2)
library(qvalue)

con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
edata = edata[rowMeans(edata) > 100,]
fdata = fData(mp)

de = DESeqDataSetFromMatrix(edata, pdata, ~study)
glm_de = DESeq(de)
result_de= results(glm_de)
#hist(result_nb$stat)

# using limma test the differences
edata = log2(as.matrix(edata) + 1)
mod = model.matrix(~ as.factor(pdata$study))
fit_limma = lmFit(edata,mod)
ebayes_limma = eBayes(fit_limma)
top = topTable(ebayes_limma,number=dim(edata)[1], sort.by="none")
# correlation in the statistics between two analyses
cor(result_de$stat, top$t)

# make an MA-plot
y = cbind(result_de$stat, top$t)
limma::plotMA(y)


fp_bonf_de = p.adjust(result_de$pvalue,method="bonferroni")
hist(fp_bonf_de,col=3)
quantile(fp_bonf_de)

fp_bonf_limma = p.adjust(top$P.Value,method="bonferroni")
hist(fp_bonf_limma,col=3)
quantile(fp_bonf_limma)

sum(fp_bonf_de < 0.05)
sum(fp_bonf_limma < 0.05)

