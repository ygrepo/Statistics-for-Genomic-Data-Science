## ----load----------------------------------------------------------------
library(devtools)
library(Biobase)
library(DESeq2)


rm(list=ls())

con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bottomly_eset.RData")
load(file=con)
close(con)
bot = bottomly.eset
pdata_bot=pData(bot)

con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
pdata_bm=pData(bm)
# Examine covariates
## BodyMap
table(pdata_bm$tissue.type,pdata_bm$gender)
table(pdata_bm$tissue.type,pdata_bm$age)
table(pdata_bm$tissue.type,pdata_bm$race)
## Bottomly
table(pdata_bot$strain,pdata_bot$experiment.number)
table(pdata_bot$strain,pdata_bot$lane.number)

## ----session_info--------------------------------------------------------
devtools::session_info()

