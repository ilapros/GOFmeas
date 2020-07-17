## ---- echo = FALSE------------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "## ")
# par(bg = "white")

## -----------------------------------------------------------------------------
library(GOFmeas)
library(lmom)
set.seed(154)
# generate a region with (mostly) GEV distribution
# this is a list
stats <- list(quagev(runif(30,0,1),c(289, 45 ,-0.22)),
              quagev(runif(45,0,1),c(189, 40 ,-0.15)),
              quagev(runif(25,0,1),c(122, 10 ,-0.24)),
              quagev(runif(43,0,1),c( 59,  8 ,-0.18)),
              quagev(runif(32,0,1),c( 62, 10 ,-0.21)),
              quagev(runif(28,0,1),c( 91,  9 ,-0.25)),
              quaglo(runif(27,0,1),c(202, 25 ,-0.17)))
set.seed(15544) # needed to ensure the results are always the same
tt <- GOFmeasures(stats)
class(tt) # the GOFmeas class
# print and summary available for the GOFmeas class
tt
summary(tt)

## -----------------------------------------------------------------------------
ob1 <- GOFmeasures(stats); ob1
ob2 <- GOFmeasures(stats); ob2

## -----------------------------------------------------------------------------
head(ob1$mcmom)
GOFmeasures(stats, mcmom = tt$mcmom)

## -----------------------------------------------------------------------------
obj3 <- GOFmeasures(stats, Nsim = 1200)
obj3
dim(obj3$mcmom)

## ----figs,fig.width = 5, fig.height = 5,dev='png',fig.path="f"----------------
GOFmeasures(stations = stats, mcmom = obj3$mcmom, plot = TRUE, conf.lev = 0.9)
# change the confidence level used to testing
GOFmeasures(stations = stats, mcmom = obj3$mcmom, plot = TRUE, conf.lev = 0.5)

# for lower levels the ellipse shrinks: less distributions accepted 

## -----------------------------------------------------------------------------
# only the HW measure
GOFmeasures(stations = stats, type="HW_GOF", mcmom = obj3$mcmom)
# only the KP measure
GOFmeasures(stations = stats, type="KP_GOF", mcmom = obj3$mcmom)

## -----------------------------------------------------------------------------
ll <- do.call(rbind,lapply(stats,samlmu)) ## L-moments matrix 
ss <- unlist(lapply(stats,length))  ## sample size vector
GOFmeasures(lmom=ll, n.amax = ss)

