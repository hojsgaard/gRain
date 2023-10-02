### R code from vignette source 'grain-intro.Rnw'

###################################################
### code chunk number 1: grain-intro.Rnw:15-22
###################################################
library(knitr)
dir.create("figures")
opts_chunk$set(fig.height=2.5,
               fig.path='figures/grain-',
               warning=FALSE, message=FALSE
)
options("prompt"="> ","width"=85)


###################################################
### code chunk number 2: grain-intro.Rnw:42-45
###################################################
require(gRain)
prettyVersion <- packageDescription("gRain")$Version
prettyDate <- format(Sys.Date())


