setwd('~/Dropbox/ccb_banding/analyses/occupancy-phylo/src/from_perry')

## work on efficiency of the importance sampling step by Perry, started 2/3/15


## From "model summaries.R"
load('~/Dropbox/ccb_banding/analyses/occupancy-phylo/saved/dclone/POM_mods_lambda_0_nsp_32_V1.rdata')

## Parameter estimates for both models.
summary(mod.lam) # Lambda estimated
summary(mod.null) # Lambda=0

## Some Rhats are high, but are for spp with occupancy intercepts
## (psi.0.sp) near 1 or 0 when passed through logit link.

## Additional Model properties are stored as attributes
summary(attributes(mod.lam))

## Whether model has lambda paramter or not, used to decide which
## likelihood calculation to use in jags.LL function below
attr(mod.lam, "model.case")

## Post-dclone jags model output, monitoring random effects. Acts as
## input for jags.LL function below
attr(mod.lam, "jags.mod") ## huge

## Vector of simulated likelihood values. Output of jags.LL
## function. These models were run with 1 million draws.
attr(mod.lam, "lik.vals")

## Summary of likelihood statistics. Output of jags.LL function.
attr(mod.lam, "stat.sum")

## Running jags.LL function to caluclate LLs for specified number of
## random draws. Functions currently set to add attributes of
## likelihood values for all simulations, and a set of statistical
## summaries based on those likelihoods, to the original model
## object. Detection here should always be TRUE since imperfect
## detection is included in all these models.

### Read in functions.R file
source('functionsFaster.R')
## Run this line once to compile phyloOccLL.c
system('R CMD SHLIB phyloOccLL.c')
## Run this line each time to load the compilation result.  I think on Windows it will have a different extension.
dyn.load('phyloOccLL.so')

## debug(jags.LL.faster)
## debug(par.calc.LL.faster)

## There is now a setup step that creates some useful things.
jags.LL.setup <- setup.jags.LL.faster(dclone.mod=mod.lam,
                                      draw.num=100)
## Do this if you want profiling results
Rprof('impSampSpeed.Rout')

## try it
system.time(mod.lam <- jags.LL.faster(setup=jags.LL.setup,
                          dclone.mod=mod.lam,
                          draw.num = 100000,
                          detection=TRUE))
## If you did the Rprof call above, do the next two lines to look at results
Rprof(NULL)
summaryRprof('impSampSpeed.Rout')

## I DID NOT TEST THIS CASE
mod.null <- jags.LL(mod.null,
                    draw.num=100,
                    num.cores=4,
                    detection=T)



batch.wrapper(dclone.mod=mod.lam, batch.size=1000, n.batches=2)

library(parallel)

as.numeric(


n.cores<-7
n.batches<-5
batch.size<-100000

n.cores*n.batches*batch.size


system.time(
					
					all.valz<-mclapply(1: n.cores ,
                               batch.wrapper,
                               dclone.mod=mod.lam,
                               batch.size= batch.size, 
                               n.batches=n.batches,
                               mc.cores= n.cores)
)                               
                


mod.lam<-pom.LL(mod.lam, n.cores=7, n.batches=5)
mod.null<-pom.LL(mod.null, n.cores=7, n.batches=5)


attr(mod.lam, "stat.sum")
attr(mod.null, "stat.sum")

anova.pom(mod.lam, mod.null)


all.valz<-unlist(all.valz)

log(mean(exp(all.valz - mean(all.valz)))) + mean(all.valz)

log.mean<-function(x) {log(mean(exp(x - mean(x)))) + mean(x)}



### Random diagnostic plots

NN<-length(all.valz)
samples.subi<-matrix(NA, nrow=10000, ncol=NN)

for (N in 1:NN){
for (i in 1:10000) {
	samples.subi [i,N]<-log.mean(sample(all.valz, N, replace=T))
}}

plot(1:NN, apply(samples.subi, 2, mean))
abline(h=max(apply(samples.subi, 2, mean)), lty=2)


apply(samples.subi, 2, mean)
hist(all.valz)       

acc<-apply(samples.subi, 2, mean)

preds<-acc-mean(all.valz)

x<-1:N
mod<-lm(exp(preds)~(x)-1)
summary(mod)

log(3.49967*100000)


exp(preds)


diff(c(500, 550))

          
