library(devtools)
library(Biobase)
library(GenomicRanges)
library(SummarizedExperiment)

rm(list=ls())

data(sample.ExpressionSet, package = "Biobase")
se = makeSummarizedExperimentFromExpressionSet(sample.ExpressionSet)

assay(se)
colData(se)
rowRanges(se)
rowRanges(se)

## ----session_info--------------------------------------------------------
devtools::session_info()
