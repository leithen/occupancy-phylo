##************************************************************************
## Importance sampling code for Perry
##************************************************************************
rm(list=ls())
setwd('~/Dropbox/ccb_banding/analyses/occupancy-phylo/for_perry')

library(dclone)
library(R2jags)
library(scales)

source('functionsFaster.R') ## These are the functions you sent to us,
                            ## with a few modifications and additional
                            ## functions.  Completely new, or updated
                            ## functions begin on line 411

## compile C code and load resultant .so object
fn.1 <- 'phyloOccLL.o'
fn.2 <- 'phyloOccLL.so'
fn.3 <- 'phyloOccLL.c'
if(file.exists(fn.1)) file.remove(fn.1)
if(file.exists(fn.2)) file.remove(fn.2)
system(sprintf('R CMD SHLIB %s', fn.3))
dyn.load(fn.2)

## load the data which was simulated with lambda=0.5 and 32 species
## (lambda=0.5 is important, because as will be shown later, the model
## with lambda is the less likely model)
load('~/Dropbox/ccb_banding/analyses/occupancy-phylo/saved/dclone/POM_mods_lambda_0.5_nsp_32_V8.rdata',
     verbose=TRUE)
## load('POM_mods_lambda_0.5_nsp_32_V8.rdata', verbose=T)

####################################################################
## Calculate LL with small number of draws
####################################################################

## If this is repeated a few times the resulting p-values vary quite 
## a bit

set.seed(1)
n.batch <- 1

mod.lam <- pom.LL(mod.lam,
                  n.batches=n.batch,
                  batch.size=10000,
                  save.all.LLs=T)
round(attr(mod.lam,  'stat.sum'), 4)

mod.null <- pom.LL(mod.null,
                   n.batches=n.batch,
                   batch.size=10000,
                   save.all.LLs=T)

## Total number of draws
length(attr(mod.lam,  'lik.vals'))
length(attr(mod.null, 'lik.vals'))

## SE.LL reports the standard error, calculated in line 480 in
## functionsFaster.R
round(attr(mod.lam,  'stat.sum'), 4)
round(attr(mod.null, 'stat.sum'), 4)

## Histograms of LLs for both models. Vertical line represents mean,
## as calculated on line 473 in functionsFaster.R.  Blue
## (respectively, red) corresponds to draws from models with
## (respectively without) lambda.
hist(attr(mod.lam, 'lik.vals'), col=alpha('blue', 0.5))
abline(v=attr(mod.lam, 'stat.sum')['LL'], col='blue')
hist(attr(mod.null, 'lik.vals'), add=T, col=alpha('red', 0.5))
abline(v=attr(mod.null, 'stat.sum')['LL'], col='red')

anova.pom(mod.lam, mod.null)


######################################################################
### Repeat with larger number of draws (e.g., greater n.batch)
######################################################################

## Same behavior even with more draws

n.batch <- 200

#set.seed(1)

mod.lam <- pom.LL(mod.lam,
                  n.batches=n.batch,
                  batch.size=10000,
                  save.all.LLs=F,
                  n.cores=4)

mod.null <- pom.LL(mod.null,
                   n.batches=n.batch,
                   batch.size=10000,
                   save.all.LLs=F,
                   n.cores=4)

## Total number of draws
length(attr(mod.lam, 'lik.vals'))
length(attr(mod.null, 'lik.vals'))

## SE.LL reports the standard error, calculated in line 480 in
## functionsFaster.R
round(attr(mod.lam, 'stat.sum'), 4)
round(attr(mod.null, 'stat.sum'), 4)

## Histograms of LLs for both models. Vertical line represents mean,
## as calculated on line 478 in functionsFaster.R
# hist(attr(mod.lam, 'lik.vals'), col=alpha('blue', 0.5), breaks=100)
# abline(v=attr(mod.lam, 'stat.sum')['LL'], col='blue')
# hist(attr(mod.null, 'lik.vals'), add=T, col=alpha('red', 0.5), breaks=50)
# abline(v=attr(mod.null, 'stat.sum')['LL'], col='red')

anova.pom(mod.lam, mod.null)



## Runs with 8 million draws
# First run
# > anova.pom(mod.lam, mod.null)
#          df      AIC    logLik  dAIC   Test    L.Ratio   p-value
# mod.lam  39 7555.745 -3738.872 1.943                NA        NA
# mod.null 38 7553.802 -3738.901 0.000 1 vs 2 0.05780759 0.8099953
# Second run
# > anova.pom(mod.lam, mod.null)
#          df      AIC    logLik  dAIC   Test    L.Ratio   p-value
# mod.lam  39 7555.743 -3738.871 1.921                NA        NA
# mod.null 38 7553.822 -3738.911 0.000 1 vs 2 0.07870103 0.7790655
# Third run
# > anova.pom(mod.lam, mod.null)
#          df      AIC    logLik  dAIC   Test    L.Ratio   p-value
# mod.lam  39 7555.740 -3738.870 1.927                NA        NA
# mod.null 38 7553.813 -3738.907 0.000 1 vs 2 0.07334779 0.7865233




