library(snpStats)
library(broom)
data(for.exercise)
use <- seq(1, ncol(snps.10), 10)
sub.10 <- snps.10[,use]
snpdata = sub.10@.Data
status = subject.support$cc


# fit an additive logistic regression model to each SNP
results = rep(NA, dim(snpdata)[2])
for (i in 1:ncol(snpdata)){
  snpdata_i = as.numeric(snpdata[,i])
  snpdata_i[snpdata_i == 0] = NA
  glm_i = glm(status ~ snpdata_i, family = "binomial")
  #results[i] = tidy(glm_i)$estimate[2]
  results[i] = tidy(glm_i)$statistic[2]
}


# square the coefficients
results_coeff_square =  results^2


# correlation with the results from using snp.rhs.tests and chi.squared
glm_all = snp.rhs.tests(status ~ 1, snp.data = sub.10)
cor(results_coeff_square, chi.squared(glm_all))