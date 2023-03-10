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
  library(UsingR)
})

## ----load----------------------------------------------------------------
  library(devtools)
  library(Biobase)
  library(UsingR)

## ----install_packages, eval=FALSE----------------------------------------
#  install.packages(c("devtools","UsingR"))
#  source("http://www.bioconductor.org/biocLite.R")
#  biocLite(c("Biobase"))

## ----galton,fig.height=4,fig.width=8-------------------------------------
data(galton)
par(mfrow=c(1,2))
hist(galton$child,col="blue",breaks=100)
hist(galton$parent,col="blue",breaks=100)

## ---- dependson="galton",fig.height=4,fig.width=4------------------------
hist(galton$child,col="blue",breaks=100)

## ---- dependson="galton",fig.height=4,fig.width=4------------------------
hist(galton$child,col="blue",breaks=100)
meanChild <- mean(galton$child)
lines(rep(meanChild,100),seq(0,150,length=100),col="red",lwd=5)

## ---- dependson="galton",fig.height=4,fig.width=4------------------------
plot(galton$parent,galton$child,pch=19,col="blue")

## ---- dependson="galton",fig.height=4,fig.width=4------------------------
plot(galton$parent,galton$child,pch=19,col="blue")
near65 <- galton[abs(galton$parent - 65)<1, ]
points(near65$parent,near65$child,pch=19,col="red")
lines(seq(64,66,length=100),rep(mean(near65$child),100),col="red",lwd=4)

## ---- dependson="galton",fig.height=4,fig.width=4------------------------
plot(galton$parent,galton$child,pch=19,col="blue")
near71 <- galton[abs(galton$parent - 71)<1, ]
points(near71$parent,near71$child,pch=19,col="red")
lines(seq(70,72,length=100),rep(mean(near71$child),100),col="red",lwd=4)

## ---- dependson="lm1",fig.height=4,fig.width=4---------------------------
plot(galton$parent,galton$child,pch=19,col="blue")
lm1 <- lm(galton$child ~ galton$parent)
lines(galton$parent,lm1$fitted,col="red",lwd=3)

## ---- dependson="lm1",fig.height=4,fig.width=4---------------------------
plot(galton$parent,galton$child,pch=19,col="blue")
lines(galton$parent,lm1$fitted,col="red",lwd=3)

## ---- dependson="lm1",fig.height=4,fig.width=8---------------------------
par(mfrow=c(1,2))
plot(galton$parent,galton$child,pch=19,col="blue")
lines(galton$parent,lm1$fitted,col="red",lwd=3)
plot(galton$parent,lm1$residuals,col="blue",pch=19)
abline(c(0,0),col="red",lwd=3)

## ----session_info--------------------------------------------------------
devtools::session_info()


