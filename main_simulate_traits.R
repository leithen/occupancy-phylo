## ************************************************************
## run analysis on simulated data
## ************************************************************
setwd('~/Dropbox/ccb_banding')
rm(list=ls())
source('analyses/occupancy-phylo/src/initialize.R')
## ************************************************************

## Set tree so as to re-run multiple times but control the tree


master.simulator<-function(seed, lambda, lambda.trait, nsp){

set.seed(seed)
nsp<-nsp
num.sp <- nsp
tree <- balanced.tree(num.sp)
## pom.sim(tree=tree1)
dd.model <- make.data(nsp=num.sp,
                      tree=tree,
                      beta.trait=1,
                      psi.0.trait=0,
                      lambda.trait=lambda.trait,
                      lambda=lambda,
                      sigma.psi.beta=3,
                      sigma.psi.site=1,
                      sigma.psi.year=1,
                      sigma.p.0=1,
                      mu.p.0=-1,
                      mu.psi.0=0,
                      sigma.psi.0=2,
                      nyr=5,
                      nrep=5,
                      imperfect.detection=TRUE)

phylosig(dd.model$inputs$tree,
         dd.model$inputs$psi.beta.sp.rand, method="lambda", test=T)

phylosig(dd.model$inputs$tree,
         dd.model$inputs$beta.spp, method="lambda", test=T)


analyse <- function(dd, psi.beta.trait.init, ...) {
  z.init <- dd$Z
  hold<-as.numeric(psi.beta.trait.init)
  my.inits <- function() {
    list(Z=z.init,
    	psi.beta.trait=hold)
  }

  dd <- list(data=dd,
             inits=my.inits,
             params=get.params())
  
  run.R2jags.model(dd, ...)
}



### First run model without lambda, what's the effect of trait?
source('~/Dropbox/ccb_banding/analyses/occupancy-phylo/src/models/traits/jags_null.R')
res.null <- analyse(dd.model, ni=3000, nb=500, psi.beta.trait.init=0) 

sum.null<-res.null$BUGSoutput$summary; betaz.null<-sum.null[grep('psi.beta', row.names(sum.null)),]


## Then run with lambda is free to vary
source('~/Dropbox/ccb_banding/analyses/occupancy-phylo/src/models/traits/jags_lambda.R')
res.lam <- analyse(dd.model, ni=3000, nb=500, psi.beta.trait.init=0) 
res.lam
sum.lam <-res.lam$BUGSoutput$summary; betaz.lam<-sum.lam[grep('psi.beta', row.names(sum.lam)),]

betaz.null
betaz.lam


logit.lambda <- res.lam$BUGSoutput$sims.matrix[,'logit.lambda']
## hist(expit(logit.lambda))


## These parameters are saved to input into the second model to help
## with model mixing.
MU.lambda <- mean(logit.lambda)
TAU.lambda <- 1/sd(logit.lambda)^2

dd.model$MU.lambda <- MU.lambda 
dd.model$TAU.lambda <- TAU.lambda

###Run second model with indicator to switch lambda on and off
source('~/Dropbox/ccb_banding/analyses/occupancy-phylo/src/models/traits/jags_lambda_select.R')
res2 <- analyse(dd.model, ni=3000, nb=500, 
               psi.beta.trait.init=res.lam$BUGSoutput$summary['psi.beta.trait','mean']) 
               

## bf.pre <- res2$BUGSoutput$summary['IND.lambda','mean'] #Pull out indicator posterior

## bf <- bf.pre/(1-bf.pre) ### Calculate Bayes Factor (How many times more likely is the full model than the null?)
## bf

## pseudoLR <- 2*log(bf) ### Convert Bayes Factor to something that looks like a likelihood ratio
## pseudoLR

## pchisq(pseudoLR, df=1, lower.tail=F) ### Do a completely silly 'LR' test. This doesn't actually make any sense, but helps me visualize how 'significant' the parameter is.

## expit(MU.lambda)
## expit(res$BUGSoutput$summary['logit.lambda',c('mean','2.5%','97.5%')])


## hold<-round(res$BUGSoutput$summary[,c('mean', '2.5%', '97.5%', 'Rhat')] , 3) ; hold[grep('psi.beta', row.names(hold)),]

master.list<-list(res.null$BUGSoutput$summary, res.lam$BUGSoutput$summary, res2$BUGSoutput$summary)
return(master.list)

}


extract.summary<-function(x) {
	NOM<-x[[1]]['psi.beta.trait',c('mean', '2.5%', '97.5%', 'Rhat')]
	POM<-x[[2]]['psi.beta.trait',c('mean', '2.5%', '97.5%', 'Rhat')]

	MW<-x[[3]]['IND.lambda','mean']

	out<-data.frame(rbind(NOM, POM))
	out$lambda<-c(0, expit(x[[2]]['logit.lambda',c('mean')]))
	out$ModelWeight<-c(1-MW, MW)

	return(out)
}




system.time(
# out.list<-
master.simulator(seed=1, lambda=1, lambda.trait=0, nsp=32)
)

setwd('~/')

#file_name<-paste("~/Dropbox/ccb_banding/analyses/occupancy-phylo/saved/traits/nsp16.rdata")
file_name<-paste("~/Dropbox/ccb_banding/analyses/occupancy-phylo/saved/traits/nsp32.rdata")


seed.max<-100

108*100*4/60/60
9*100*1/60


nsp32_lambda1_lambdatrait_0<-lapply(1:seed.max, master.simulator, lambda=1, lambda.trait=0, nsp=32)
file_name<-paste("~/Dropbox/ccb_banding/analyses/occupancy-phylo/saved/traits/nsp32-1-0.rdata")
save(nsp32_lambda1_lambdatrait_0, file=file_name)


nsp32_lambda1_lambdatrait_1<-lapply(1:seed.max, master.simulator, lambda=1, lambda.trait=1, nsp=32)
file_name<-paste("~/Dropbox/ccb_banding/analyses/occupancy-phylo/saved/traits/nsp32-1-1.rdata")
save(nsp32_lambda1_lambdatrait_1, file=file_name)

nsp32_lambda0_lambdatrait_1<-lapply(1:seed.max, master.simulator, lambda=0, lambda.trait=1, nsp=32)
save(nsp32_lambda1_lambdatrait_0, 
     nsp32_lambda1_lambdatrait_1, 
     nsp32_lambda0_lambdatrait_1, file=file_name)

nsp32_lambda0_lambdatrait_0<-lapply(1:seed.max, master.simulator, lambda=0, lambda.trait=0, nsp=32)
save(nsp32_lambda1_lambdatrait_0, 
     nsp32_lambda1_lambdatrait_1, 
     nsp32_lambda0_lambdatrait_1,
     nsp32_lambda0_lambdatrait_0, file=file_name)

load(file_name, verbose=T)





### Look and analyze results
load('~/Dropbox/ccb_banding/analyses/occupancy-phylo/saved/traits/nsp16.rdata', verbose=T)

out.hold<-lapply(nsp16_lambda0_lambdatrait_0, extract.summary)

out.hold<-lapply(nsp16_lambda0_lambdatrait_1, extract.summary)

out.hold<-lapply(nsp16_lambda1_lambdatrait_0, extract.summary)

out.hold<-lapply(nsp16_lambda1_lambdatrait_1, extract.summary)

load('~/Dropbox/ccb_banding/analyses/occupancy-phylo/saved/traits/nsp32-1-1.rdata', verbose=T)
out.hold<-lapply(nsp32_lambda1_lambdatrait_1, extract.summary)

load('~/Dropbox/ccb_banding/analyses/occupancy-phylo/saved/traits/nsp32-1-0.rdata', verbose=T)
out.hold<-lapply(nsp32_lambda1_lambdatrait_0, extract.summary)

load('~/Dropbox/ccb_banding/analyses/occupancy-phylo/saved/traits/nsp32-0-1.rdata', verbose=T)
out.hold<-lapply(nsp32_lambda0_lambdatrait_1, extract.summary)

load('~/Dropbox/ccb_banding/analyses/occupancy-phylo/saved/traits/nsp32-0-0.rdata', verbose=T)
out.hold<-lapply(nsp32_lambda0_lambdatrait_0, extract.summary)

load('~/Dropbox/ccb_banding/analyses/occupancy-phylo/saved/traits/nsp64-1-0.rdata', verbose=T)
out.hold<-lapply(nsp64_lambda1_lambdatrait_0, extract.summary)



### What's the breadth of the confidence intervals
hold<-matrix(unlist(lapply(out.hold,function(x) {apply((x[,2:3]), 1, diff)})), ncol=2, byrow=T)
apply(hold, 2, mean)

bci.inf<-hold[,1]/hold[,2]
mean(bci.inf)

### Does the confidence interval contain the true value?
hold<-matrix(unlist(lapply(out.hold,function(x) {x[,2] < 1 & x[,3] > 1})), ncol=2, byrow=T)
hold<-data.frame(hold);names(hold)<-c("NOM","POM")
head(hold)
apply(hold, 2, sum)

### Which estimate is closer to the true estimate?
hold<-matrix(unlist(lapply(out.hold,function(x) {abs(x[,1]-1)})), ncol=2, byrow=T)
hold<-data.frame(hold);names(hold)<-c("NOM","POM")
apply(hold, 2, mean)

bci.inf<-hold[,1]/hold[,2]
mean(bci.inf)


### Which model is better
apply(t(apply(hold, 1, function (x) {x==min(x)})), 2, sum)



### POM absolutely necessary when there's phylogenetic signal in both trait, and residuals. Helps when there's phylogenetic signal in just the residuals


### Quick assay of type I error rate. What model weight cut off is right???
hold<-matrix(unlist(lapply(out.hold,function(x) {x[,6]})), ncol=2, byrow=T)
length(which(hold[,2]>0.95))/100
length(which(hold[,2]>0.75))/100



length(which(hold[,2]>0.5))/100
length(which(hold[,2]>0.6667))/100
length(which(hold[,2]>0.75))/100

length(which(hold[,2]>0.90))/100
length(which(hold[,2]>0.95))/100


