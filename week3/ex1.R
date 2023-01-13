
rm(list=ls())
tropical=  c('darkorange', 'dodgerblue', 'hotpink', 'limegreen', 'yellow')
palette(tropical)


library(devtools)
library(Biobase)
library(snpStats)
library(broom)
library(MASS)
library(DESeq2)




data(for.exercise)
use <- seq(1, ncol(snps.10), 10)
sub.10 <- snps.10[,use]
snpdata = sub.10@.Data
status = subject.support$cc

snp = as.numeric(snpdata[,3])
snp[snp==0] = NA

lm1 = lm(status ~ snp)
tidy(lm1)

glm2 = stats::glm(status ~ snp,family="binomial")
tidy(glm2)

par(mfrow=c(1,2))
plot(status ~ snp,pch=19)
abline(lm1,col="darkgrey",lwd=5)
plot(glm2$residuals)

