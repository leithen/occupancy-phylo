## ************************************************************
## run analysis on simulated data using data cloning
## ************************************************************
setwd('~/Dropbox/ccb_banding')
rm(list=ls())
source('analyses/occupancy-phylo/src/initialize.R')
## ************************************************************

## ************************************************************
## create data
set.seed(1)
nsp <- 16
tree <- balanced.tree(nsp)
## pom.sim(tree=tree1)
dd.model <- make.data(nsp=nsp,
                      tree=tree,
                      ## beta.trait=0,
                      lambda.trait=0,
                      lambda=0.5,
                      sigma.psi.beta=1,
                      imperfect.detection=TRUE)

## load model file
source(file.path('analyses/occupancy-phylo/src/models',
                 'IvesHelmus_p_mII_cts_trait_lambda_rand.R'))

analyse <- function(dd) {

  z.init <- dd$z.init
  my.inits <- function() {
    list(Z=z.init)
  }
  ## dd$Z <- dd$z.init

  dd <- list(data=dd,
             inits=my.inits,
             params=get.params())
  
  run.R2jags.model(dd)
}

res <- analyse(dd.model)

cols <- c('mean', '2.5%', '97.5%', 'Rhat', 'n.eff')
res$BUGSoutput$summary[c('lambda',
                         'sigma.psi.beta'), cols]
## ************************************************************
