## ************************************************************
## run analysis on simulated data
## ************************************************************
setwd('~/Dropbox/ccb_banding')
rm(list=ls())
source('analyses/occupancy-phylo/src/initialize.R')
## ************************************************************

## Set tree so as to re-run multiple times but control the tree
set.seed(1)
num.sp <- 32
tree <- balanced.tree(num.sp)
## pom.sim(tree=tree1)
dd.model <- make.data(nsp=num.sp,
                      tree=tree,
                      beta.trait=0,
                      lambda.trait=0,
                      lambda=1,
                      sigma.psi.beta=5,
                      sigma.psi.site=1,
                      sigma.psi.year=1,
                      sigma.p.0=1,
                      mu.p.0=0,
                      mu.psi.0=0.5,
                      sigma.psi.0=2,
                      imperfect.detection=TRUE)

phylosig(dd.model$inputs$tree,
         dd.model$inputs$psi.beta.sp, method="lambda", test=T)

## load model file
source('analyses/occupancy-phylo/src/models/jags_lambda.R')
##  source(file.path('analyses/occupancy-phylo/src/models',
##                  'IvesHelmus_p_mII_cts_trait_lambda_rand.R'))

analyse <- function(dd, ...) {

  z.init <- dd$Z
  my.inits <- function() {
    list(Z=z.init)
  }

  dd <- list(data=dd,
             inits=my.inits,
             params=get.params())
  
  run.R2jags.model(dd, ...)
}

## Run first model in which lambda is free to vary
res <- analyse(dd.model, ni=6000, nb=3000) 
res
logit.lambda <- res$BUGSoutput$sims.matrix[,'logit.lambda']
hist(expit(logit.lambda))
## traceplot(res)

## These parameters are saved to input into the second model to help
## with model mixing.
MU.lambda <- mean(logit.lambda)
TAU.lambda <- 1/sd(logit.lambda)^2

dd.model$MU.lambda <- MU.lambda 
dd.model$TAU.lambda <- TAU.lambda

source('~/Dropbox/ccb_banding/analyses/occupancy-phylo/src/models/jags_lambda_select.R')
res2 <- analyse(dd.model, ni=6000, nb=3000) ###Run second model with indicator to switch lambda on and off
res2

bf.pre <- res2$BUGSoutput$summary['IND.lambda','mean'] #Pull out indicator posterior

bf <- bf.pre/(1-bf.pre) ### Calculate Bayes Factor (How many times more likely is the full model than the null?)
bf

pseudoLR <- 2*log(bf) ### Convert Bayes Factor to something that looks like a likelihood ratio
pseudoLR

pchisq(pseudoLR, df=1, lower.tail=F) ### Do a completely silly 'LR' test. This doesn't actually make any sense, but helps me visualize how 'significant' the parameter is.

expit(MU.lambda)
expit(res$BUGSoutput$summary['logit.lambda',c('mean','2.5%','97.5%')])



#traceplot(res2)


## ************************************************************
##############################################################
### Create a single model to compare increased accuracy of slope estimates 
### when lambda is estimated
##############################################################
## ************************************************************
setwd('~/Dropbox/ccb_banding')
rm(list=ls())
source('analyses/occupancy-phylo/src/initialize.R')
## ************************************************************


pom.coefs<-function(dd.simulated, res, parameter.group, plot=T, color.gradient=c("red", "blue", "orange"), ...) {

tree<-dd.simulated$inputs$tree

model.summary<-res$BUGSoutput$summary[,c('mean', '2.5%', '97.5%','Rhat') ]

psi.beta.sp.ID<-grep(parameter.group, row.names(model.summary))

coefs<-model.summary[psi.beta.sp.ID,]

coefs

######################################################
###Add Phylogeny
######################################################
require(ape)

#Load tree


row.names(coefs)<-as.character(tree$tip.label)


tip.labels.in.order<-tree$tip.label[tree$edge[,2][tree$edge[,2] %in% 1:dim(coefs)[1]]] #Have to reorder that slopes final plots in same order as tree

coefs.final <-coefs[match( tip.labels.in.order, row.names(coefs)),]

if (plot==TRUE) {
############################################################
###With Tip Labels
############################################################
par(mar=c(6,0,2,0))
#layout(matrix(c(1,2,3), nrow=1), widths=c(1,1,.2))
layout(matrix(c(1,2), nrow=1), widths=c(1, 1))

colz.tree<-colorize(coefs.final[,1], colors= color.gradient)

plot(tree, show.tip.label=T, cex=.7, label.offset=2)
tiplabels(pch=21, bg=colz.tree)

colz.plot<-colorize(coefs.final[,1], colors=color.gradient)

plot(y=1:dim(coefs.final)[1], x= coefs.final[,1], bty='n', yaxt='n', ylab="", pch=21, bg= colz.plot, ...)
arrows(y0=1:dim(coefs.final)[1], x0= coefs.final[,2], x1= coefs.final[,3], length=0)
abline(v=0, lty=2)
}

return(coefs.final)
}


load('~/Dropbox/ccb_banding/analyses/occupancy-phylo/saved/RJMCMC/Det25/POM_mods_lambda_1_nsp_32_RJMCMC_V3.rdata', verbose=T)

res.select
res.full

dd.model<-dd.simulated

analyse <- function(dd, ...) {

  z.init <- dd$Z
  my.inits <- function() {
    list(Z=z.init)
  }

  dd <- list(data=dd,
             inits=my.inits,
             params=c(get.params(), 'psi.beta.sp', 'p.0')) ### Updating here to track a few additional parameters
  
  run.R2jags.model(dd, ...)
}


source('analyses/occupancy-phylo/src/models/jags_lambda.R')


## Run first model in which lambda is free to vary
res.lambda <- analyse(dd.model, ni=10000, nb=5000, nt=5) 

round(res.lambda$BUGSoutput$summary[,c(1,2,3,7,8,9)], 2)

source('~/Dropbox/ccb_banding/analyses/occupancy-phylo/src/models/jags_null.R')
res.null <- analyse(dd.model, ni=10000, nb=5000, nt=1) 
round(res.null$BUGSoutput$summary[,c(1,2,3,7,8,9)], 5)


pom.coefs(dd.simulated, res.lambda, "psi.beta.sp", plot=T, color.gradient=c("gold", "darkgreen", "darkblue"), xlim=c(-10,7))

null.coefs<-pom.coefs(dd.simulated, res.null, "psi.beta.sp", plot=F)

points(y=1:dim(null.coefs)[1]-0.3, x= null.coefs[,1], pch=1, col="green", lwd=2)
arrows(y0=1:dim(null.coefs)[1]-0.3, x0= null.coefs[,2], x1= null.coefs[,3], length=0, col="green")

points(y=1:dd.simulated$inputs$nsp-0.15, x=dd.simulated$inputs$psi.beta.sp, pch=22, lwd=2, cex=1.3, col="red")





######################################################################
### Attempting computationally intensive parametric bootstrapping method. Shoot me now.
######################################################################

## ************************************************************
setwd('~/Dropbox/ccb_banding')
rm(list=ls())
source('analyses/occupancy-phylo/src/initialize.R')
## ************************************************************

## Set tree so as to re-run multiple times but control the tree

mcmc.null.dist<-function(x){
set.seed(x)
num.sp <- 32
tree <- balanced.tree(num.sp)
## pom.sim(tree=tree1)
dd.model <- make.data(nsp=num.sp,
                      tree=tree,
                      beta.trait=0,
                      lambda.trait=0,
                      lambda=0,
                      sigma.psi.beta=3,
                      sigma.psi.site=1,
                      sigma.psi.year=1,
                      sigma.p.0=1,
                      mu.p.0=0,
                      mu.psi.0=0.5,
                      sigma.psi.0=2,
                      imperfect.detection=TRUE)

phylosig(dd.model$inputs$tree,
         dd.model$inputs$psi.beta.sp, method="lambda", test=T)

analyse <- function(dd, ...) {

  z.init <- dd$Z
  my.inits <- function() {
    list(Z=z.init)
  }

  dd <- list(data=dd,
             inits=my.inits,
             params=get.params())
  
  run.R2jags.model(dd, ...)
}


## load model file
source('analyses/occupancy-phylo/src/models/jags_lambda.R')


## Run first model in which lambda is free to vary

res <- analyse(dd.model, ni=2000, nb=1000, nt=1) 
res
logit.lambda <- res$BUGSoutput$sims.matrix[,'logit.lambda']

hist(expit(logit.lambda), main=x)
abline(v=expit(mean(logit.lambda)), col='red', lwd=3)

## traceplot(res)

logit.lambda.obs<- mean(res$BUGSoutput$sims.matrix[,'logit.lambda'])

## These parameters are saved to input into the second model to help
## with model mixing.
source('~/Dropbox/ccb_banding/analyses/occupancy-phylo/src/models/jags_null.R')
res2 <- analyse(dd.model, ni=2000, nb=1000, nt=1) ###Run second model with indicator to switch lambda on and off
res2


res2$BUGSoutput$summary[,c('mean', 'Rhat')]
### traceplot(res2)

psi.0s<-res2$BUGSoutput$summary[grep('psi.0.sp',row.names(res2$BUGSoutput$summary)),'mean']



alpha<- 0.05
sim.max<-100
p.holding<-1
minimum<-10

min.regulator<-c(rep(1000, minimum), rep(1, sim.max-minimum))

i<-1
lambda.null.dist<-numeric()

while ((min.regulator[i]*p.holding) > 0.05 && i<sim.max) {

dd.model.hold <- make.data(nsp=length(psi.0s),
                      tree=tree,
                      beta.trait=0,
                      lambda.trait=0,
                      lambda=0,
                      sigma.psi.beta=res2$BUGSoutput$summary['sigma.psi.beta','mean'],
                      sigma.psi.site= res2$BUGSoutput$summary['sigma.psi.site','mean'],
                      sigma.psi.year= res2$BUGSoutput$summary['sigma.psi.year','mean'],
                      sigma.p.0= res2$BUGSoutput$summary['sigma.p.0','mean'],
                      mu.p.0= res2$BUGSoutput$summary['mu.p.0','mean'],
                      psi.0s=psi.0s,
                      use.random.int=FALSE,
                      imperfect.detection=TRUE)

## load model file
source('analyses/occupancy-phylo/src/models/jags_lambda.R')


## Run first model in which lambda is free to vary

res.temp <- analyse(dd.model.hold, ni=500, nb=100, nt=1) 

lambda.null.dist[i]<-res.temp$BUGSoutput$summary['logit.lambda','mean']
print(paste(date(),"------",i, "of model", x))

bigger<-length(which(lambda.null.dist>logit.lambda.obs))
total<-length(lambda.null.dist)

print(binom.test(bigger, total, p=alpha, conf.level=0.95))


p.holding<-binom.test(bigger, total, p=alpha, conf.level=0.99)$p.value

i<-i+1

}


return(list(dd.model,res, lambda.null.dist, bigger, total))

}



system.time(
	out<-lapply(1:100, mcmc.null.dist)
)

save(out, file="analyses/occupancy-phylo/saved/parametric bootstrap/lambda0_100.rdata")

extract.estimated.p<-function(data,x) {data[[x]][[4]]/data[[x]][[5]]}

ps<-sapply(1:100, extract.estimated.p, data=out)

hist(ps)

sigs<-ps [which(ps<0.05)] ### Not a complete success. Appearent false positive rate is 9% instead of 5%, but at least it's within the realm of reason. Possible that I'm not doing this correctly, and if the random effects were dealt with correctly then this could legitimatly work. Takes an obnoxious amount of time to run though.

length(sigs)

binom.test(9, 100, p=0.05)

ps [which(ps<0.2)]

######################################################################
### OLD
######################################################################



cols <- c('mean', '2.5%', '97.5%', 'Rhat', 'n.eff')
res$BUGSoutput$summary[c('beta.trait',
                         'lambda',
                         'sigma.psi.beta'), cols]







cols <- c('mean', '2.5%', '97.5%', 'Rhat', 'n.eff')
res$summary[c('sigma.psi.beta.non.phylo',
              'sigma.psi.beta'), cols]

hist(res$bugs$BUGSoutput$sims.matrix[,'lambda'],
     col='red', breaks=0:200/200)

hist(res$bugs$BUGSoutput$sims.matrix[,'sigma.psi.beta'],
     col='red', breaks=-50:300/300)

hist(res$bugs$BUGSoutput$sims.matrix[,'sigma.psi.beta']^2/
     (res$bugs$BUGSoutput$sims.matrix[,'sigma.psi.beta']^2 +
      res$bugs$BUGSoutput$sims.matrix[,'sigma.psi.beta.non.phylo']^2),
     col='red', breaks=100)



## simulate with no phylogenetic covariance in random effect of
## species but a phylogenetic signal in the species' trait values
##
## run a model without beta.trait and we should infer a large-ish
## lambda
##
## run a model with beta.trait and we should infer a smaller lambda









## runif(1)
## set.seed(1)

## run and save for various values of lambda
run.lambda <- function(nsp) {
  vary.param(save.dir=sprintf('nsp_%d_p_05', nsp),
             param='lambda',
             values=seq(from=0, to=1, length=11),
             nsp=nsp,
             nreps=10)
  NULL
}  
run.lambda(nsp=32)
run.lambda(nsp=64)
run.lambda(nsp=128)

## run and save for various values of sigma.psi.beta
run.sigma.psi.beta <- function(nsp, lambda) {
  if(lambda==0) save.dir <- sprintf('nsp_%s_p_05_lambda_0', nsp)
  if(lambda==0.5) save.dir <- sprintf('nsp_%s_p_05_lambda_05', nsp)
  if(lambda==1) save.dir <- sprintf('nsp_%s_p_05_lambda_1', nsp)
  vary.param(save.dir=save.dir,
             param='sigma.psi.beta',
             values=exp(seq(from=log(0.1), to=log(100), length=10)),
             nsp=nsp,
             nreps=10,
             lambda=lambda)
  NULL
}

## with lambda=0
run.sigma.psi.beta(32, lambda=0)
run.sigma.psi.beta(64, lambda=0)
run.sigma.psi.beta(128, lambda=0)

## with lambda=0.5
run.sigma.psi.beta(32, lambda=0.5)
run.sigma.psi.beta(64, lambda=0.5)
run.sigma.psi.beta(128, lambda=0.5)

## with lambda=1
run.sigma.psi.beta(32, lambda=1)
run.sigma.psi.beta(64, lambda=1)
run.sigma.psi.beta(128, lambda=1)







load('~/Dropbox/ccb_banding/analyses/occupancy-phylo/saved/simulated/sigma.psi.beta/nsp_64_p_05_lambda_1/10.000/10.RData')
cols <- c('mean', '2.5%', '97.5%', 'Rhat', 'n.eff')
show <- c('sigma.psi.beta',
          'lambda')
summary <- res$model.out$bugs$BUGSoutput$summary
summary[show,cols]

expit(summary['mu.p.0','mean'])


load('~/Dropbox/ccb_banding/analyses/occupancy-phylo/saved/simulated/lambda_64/0.200.RData')
cols <- c('mean', '2.5%', '97.5%', 'Rhat', 'n.eff')
show <- c('mu.psi.0',
          'mu.p.0',
          'psi.0.trait')
summary <- res$model.out$bugs$BUGSoutput$summary
summary[show,cols]
expit(summary['mu.p.0','mean'])


## TO DO:

## try with different slopes mu.psi.0, mu.p.0, mu.beta (mu.beta, being
## the most interesting)

## look at nsite, nsp, nyear, nrep once we know its working

## next consider non-zero beta.trait.  The interesting thing here is
## to contrast cases with and without phylogeny (lambda=0 vs
## lambda=1).  In doing this, it will be necessary to decide whether
## occurrence is simulated WITH additional phylogenetic structure,
## beyond that that would arise from the trait alone, or only with
## that which arises from the trait.
