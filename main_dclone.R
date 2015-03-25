## ************************************************************
## run analysis on simulated data using data cloning
## ************************************************************
setwd('~/Dropbox/ccb_banding')
rm(list=ls())
source('analyses/occupancy-phylo/src/initialize.R')
## ************************************************************

## ************************************************************
## STEP 0: Simulate data.
set.seed(1)
nsp <- 8
lambda.input<-1
tree <- balanced.tree(nsp)

file.ext<-paste("lambda",lambda.input,"nsp",nsp,sep="_")

dd.simulated <- make.data(nsp=nsp,
                          tree=tree,
                          beta.trait=0,
                          lambda.trait=0,
                          lambda=lambda.input,
                          sigma.psi.beta=3,
                          sigma.psi.site=1,
                          sigma.psi.year=1,
                          sigma.p.0=1,
                          mu.p.0=1,
                          mu.psi.0=0.5,
                          sigma.psi.0=2,
                          imperfect.detection=TRUE)

phylosig(dd.simulated$inputs$tree, dd.simulated$inputs$beta.sp, 'lambda')
fitContinuous(dd.simulated$inputs$tree, dd.simulated$inputs$beta.sp, model='OU') #These parameter estimates will differ from the POM estimates because the tree is scaled to be a correlation matrix in POM, but is left with branch lengths in years (i.e. covariance matrix) here.
fitContinuous(dd.simulated$inputs$tree, dd.simulated$inputs$beta.sp, model='lambda')
############################################################
## *********************************************************
## *********************************************************
## RUN MODS AND CALC LL
## *********************************************************
## *********************************************************
############################################################
## Full phylogenetic signal model
############################################################

## Establish general dclone control
## paramaters to run for each model
nc<-c(2)
n.iter<-1000
n.adapt<-1000
n.update<-1000
num.cores<-3

## *** CRITICAL NOTE ::: KEEP IN LONG TERM KNOWLEDGE BASE *** ##
## I tried running a few models with short run times (n.iter, 
## n.adapt at 100) and was getting "Error in chol.default(W) : 
##  the leading minor of order 30 is not positive definite"
## I tracked the problem to psi.sp.0. The problem seemed to get
## worse the more species were simulated. Inference: It seems the 
## problem was with one or two species that were hard to estimate
## parameter values for. By increaseing n.iter, and n.adapt the
## problem went away. SO REMEMBER IF YOU GET AN chol.default(W)
## error ENSURE THAT THE n.iters etc ARE SUFFICIENT!

draws<-1e4 #For importance sampling

##Run dclone for MLEs
dclone.out.lam <- run.dc.fit.update(dd.model=dd.simulated,
                                    nc=nc,
                                    n.iter=n.iter,
                                    n.adapt=n.adapt,
                                    n.update=n.update,
                                    num.cores=num.cores,
                                    model.case='phylo.lam')

summary(dclone.out.lam)

##*******************************************************************
### LOAD OLD MODEL FOR DEBUGGING--This is to examine the need to 
### explicity model the Z state
##*******************************************************************
load('~/Dropbox/ccb_banding/analyses/occupancy-phylo/saved/dclone/POM_mods_lambda_1_nsp_32_V1.rdata', verbose=T)


str(mod.lam)
##Given MLEs calculate random effect distributions
mod.lam <- run.jags.after.dclone(dclone.mod= dclone.out.lam)

mod.lam <- run.jags.after.dclone(dclone.mod= mod.lam)
mod.lam.debug <- run.jags.after.dclone(dclone.mod= mod.lam, ALTERNATIVE_DELETE_ME=T)


debug.vec<-attr(mod.lam.debug, 'jags.mod')$bugs$BUGSoutput$summary[,'mean']
normal.vec<-attr(mod.lam, 'jags.mod')$bugs$BUGSoutput$summary[,'mean']

pz.debug<-debug.vec[grep("p.0", names(debug.vec))]
pz.normal<-normal.vec[grep("p.0", names(normal.vec))]

beta.debug<-debug.vec[grep("psi.beta.sp", names(debug.vec))]
beta.normal<-normal.vec[grep("psi.beta.sp", names(normal.vec))]


plot(pz.debug, pz.normal)
abline(a=0, b=1)


plot(beta.debug, beta.normal)
abline(a=0, b=1)

max(attr(mod.lam,"jags.mod")$bugs$BUGSoutput$summary[,'Rhat'])
max(attr(mod.lam.debug,"jags.mod")$bugs$BUGSoutput$summary[,'Rhat'])
traceplot(attr(mod.lam.debug,"jags.mod"))
traceplot(as.mcmc(attr(mod.lam,"jags.mod")$bugs))

round(attr(mod.lam, 'jags.mod')$bugs$BUGSoutput$summary[,c('mean', '2.5%', '97.5%', 'Rhat')], 3)
round(attr(mod.lam.debug, 'jags.mod')$bugs$BUGSoutput$summary[,c('mean', '2.5%', '97.5%', 'Rhat')], 3)



##*******************************************************************
###END DEBUGGING
### Conclusion: Z-state maybe there only to increase speed of 
### convergence resulting from setting initial values?
##*******************************************************************



##Given MLEs calculate random effect distributions
mod.lam <- run.jags.after.dclone(dclone.mod= dclone.out.lam)


## Check for convergence
max(attr(mod.lam,"jags.mod")$bugs$BUGSoutput$summary[,'Rhat'])

##Based on random effect distributions calculate LL
mod.lam<-jags.LL(mod.lam, draw.num=draws, num.cores=7, detection=T)

##Check statistics summary of the model
attr(mod.lam, "stat.sum")

### get mles of parameters
summary(mod.lam)

############################################################
## No phylogenetic signal
############################################################

dclone.out.null <- run.dc.fit.update(dd.model=dd.simulated,
                         nc=nc,
                         n.iter=n.iter,
                         n.adapt=n.adapt,
                         n.update=n.update,
                         num.cores=num.cores,
                         model.case='phylo.null')
                                                  
mod.null <- run.jags.after.dclone(dclone.mod= dclone.out.null)

## Check for convergence
max(attr(mod.null,"jags.mod")$bugs$BUGSoutput$summary[,'Rhat'])

##Calc LLs
mod.null<-jags.LL(mod.null, draw.num= draws, num.cores=7, detection=T)



############################################################
## Full phylogenetic signal model using OU model
############################################################

dclone.out.ou <- run.dc.fit.update(dd.model=dd.simulated,
                         nc=nc,
                         n.iter=n.iter,
                         n.adapt=n.adapt,
                         n.update=n.update,
                         num.cores=num.cores,
                         model.case='phylo.ou')

##Given MLEs calculate random effect distributions
mod.ou <- run.jags.after.dclone(dd=dd.simulated,
                             dclone.mod= dclone.out.ou)
## Check for convergence
max(attr(mod.ou,"jags.mod")$bugs$BUGSoutput$summary[,'Rhat'])

##Based on random effect distributions calculate LL
mod.ou<-jags.LL(mod.ou, draw.num=draws, num.cores=7, detection=T)


############################################################
## No detection but phylogenetic signal
############################################################

dclone.out.no.det <- run.dc.fit.update(dd.model=dd.simulated,
                         nc=nc,
                         n.iter=n.iter,
                         n.adapt=n.adapt,
                         n.update=n.update,
                         num.cores=num.cores,
                         model.case='phylo.lam.no.det')                         

mod.no.det <- run.jags.after.dclone(dd=dd.simulated,
                             dclone.mod= dclone.out.no.det)
max(attr(mod.no.det,"jags.mod")$bugs$BUGSoutput$summary[,'Rhat'])


mod.no.det<-jags.LL(mod.no.det, draw.num=draws, num.cores=7, detection=F)


###Save models
lr.save.name<-paste("analyses/occupancy-phylo/saved/dclone/POM_mods_",file.ext,".rdata", sep="")

save(mod.lam, mod.null, mod.ou , mod.no.det, file=lr.save.name)

##******************************************************
## Load And View Results
##******************************************************
load('~/Dropbox/ccb_banding/analyses/occupancy-phylo/saved/dclone/POM_mods_lambda_0.9_nsp_32.rdata', verbose=T)

### get mles of parameters
summary(mod.lam)
summary(mod.ou)
summary(mod.null)
summary(mod.no.det)

### Test for significant differences between models
anova.pom(mod.lam, mod.null)
anova.pom(mod.ou, mod.null)
anova.pom(mod.lam, mod.no.det)
anova.pom(mod.lam, mod.ou, mod.null, mod.no.det)


### Diagnostic plots--Check that likelihood means are assymptoting:
imp.samp.plot(mod.lam, mod.ou, mod.null, mod.no.det)





anova.pom(mod.lam, mod.ou)
imp.samp.plot(mod.lam, mod.null)

