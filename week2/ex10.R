## ----global_palette, results = 'asis'------------------------------------
rm(list=ls())
tropical=  c('darkorange', 'dodgerblue', 'hotpink', 'limegreen', 'yellow')
palette(tropical)

## ----global_options,warning=FALSE,message=FALSE--------------------------
## see ch. 10 Hooks of Xie's knitr book
library(knitr)
knit_hooks$set(setPch = function(before, options, envir) {
  if(before) par(pch = 19)
})
opts_chunk$set(setPch = TRUE)

## ----global_plot,warning=FALSE, message=FALSE----------------------------
knitr::opts_chunk$set(fig.width=5, fig.height=5, size="footnotesize",
                      warning=FALSE, message=FALSE)
knitr::knit_hooks$set(small.mar = function(before, options, envir) {
  if (before) graphics::par(mar = c(5,5,1.5,1))
})

## ----load_hidden, echo=FALSE, results="hide", warning=FALSE--------------
suppressPackageStartupMessages({
  library(devtools)
  library(Biobase)
  library(sva)
  library(bladderbatch)
  library(snpStats)
})

# BiocManager::install("sva")


## ----load----------------------------------------------------------------
library(devtools)
library(Biobase)
library(sva)
library(bladderbatch)
library(snpStats)

con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
edata = exprs(bm)
pdata_bm=pData(bm)

set.seed(33353)
pdata_bm = na.omit(pdata_bm)
edata = edata[,rownames(pdata_bm), drop=FALSE]
edata = log2(edata + 1)
edata = edata[rowMeans(edata) > 1, ]

mod = model.matrix(~pdata_bm$age,data=pdata_bm)
mod0 = model.matrix(~1, data=pdata_bm)

sva1 = sva(edata,mod,mod0,n.sv=2)

cor(sva1$sv, pdata_bm$age)


# correlation between surrogate for batch and race
cor(sva1$sv, as.numeric(pdata_bm$race))
```
```{r}
# correlation between surrogate for batch and gender
cor(sva1$sv, as.numeric(pdata_bm$gender))
```
