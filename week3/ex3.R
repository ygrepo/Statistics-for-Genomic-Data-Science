rm(list=ls())
tropical=  c('darkorange', 'dodgerblue', 'hotpink', 'limegreen', 'yellow')
palette(tropical)
library(devtools)
library(Biobase)
library(sva)
library(bladderbatch)

library(snpStats)
library(broom)
data(for.exercise)
use <- seq(1, ncol(snps.10), 10)
sub.10 <- snps.10[,use]
snpdata = sub.10@.Data
status = subject.support$cc

# fit a logistic regression model
snp10 = as.numeric(snpdata[,10])
snp10[snp10==0] = NA
glm10 = glm(status ~ snp10, family="binomial")
tidy(glm10)

snp10_dom = (snp10 == 2)
glm10_dom = glm(status ~ snp10_dom, family="binomial")
tidy(glm10_dom)

par(mfrow=c(1,2))
boxplot(glm10$residuals, status, add=T, border=1:2)
boxplot(glm10_dom$residuals, status, add=T, border=1:2)