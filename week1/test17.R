library(devtools)
library(Biobase)
library(dendextend)

rm(list=ls())

con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
edata = exprs(bm)

row_sums = rowSums(edata)
edata2 = edata[order(-row_sums),]
index = 1:500
heatmap(edata2[index,],Rowv=NA,Colv=NA)

index2 = which(rank(row_sums) < 500 )
heatmap(edata[index2,],Colv=NA)

index3 = which(rank(-row_sums) < 500 )
heatmap(edata[index3,],Rowv=NA,Colv=NA)
        
heatmap(edata2[index3,],Rowv=NA,Colv=NA)
