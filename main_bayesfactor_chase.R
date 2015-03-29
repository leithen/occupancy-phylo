## ************************************************************
## prepare matrices
## ************************************************************
setwd('~/Dropbox/ccb_banding')
rm(list=ls())
source('analyses/multi-species-phylo/src/initialize.R')
## source('analyses/occupancy-phylo/src/initialize.R') ## The corresponding prep 
                                                ## file is screwing things up...
## ************************************************************

path <- 'data/for_analysis/chase'
dd.model <-
  prep.data(mat.path=file.path(path,
              'occupancy/ccb_2004_2012_mat.RData'),
            case.data='chase-late',
            case.site.traits='cat',
            case.sp.traits='all',
            threshold=25,
            reps=1:6,
            drop.migrants=FALSE)
## add VCOV mat
load('data/original/chase/ChaseVCOV.rdata')
spp.keep <- dimnames(dd.model$X)$species
VCOV.tmp <- ChaseVCOV[rownames(ChaseVCOV) %in% spp.keep,
                      colnames(ChaseVCOV) %in% spp.keep]
dd.model$VCOV <- VCOV.tmp[match(rownames(VCOV.tmp), spp.keep),
                          match(colnames(VCOV.tmp), spp.keep)]

dd.model$ID <- diag(1, dim(dd.model$VCOV)[1])

##Temporary fix for lack of env data
site.data<-read.csv("~/Dropbox/ccb_banding/data/original/chase/site_specific_spatial_covariates.csv")

dd.model$env<-site.data$Percent.forest.cover.within.100.meter.radius.of.center-0.5 ### Minus 0.5 'scales' deforestation

## Add info on coffee non-coffee for detection covar
dd.model$site.dummy<-rep(0, 18)
dd.model$site.dummy[dd.model$site.traits=="Coffee"]<-1

## Add info on net position for detection covar
terrain.info<-read.csv("~/Dropbox/ccb_banding/data/original/chase/Variable to explain detection bias.csv")
terrain.info<-terrain.info[terrain.info[,1] %in% dimnames(dd.model$X)$site,]
dd.model$site.cont.covar<-as.numeric(scale(terrain.info[match(terrain.info[,1], dimnames(dd.model$X)$site) , 2]))


save(dd.model, file=file.path(path, 'occupancy/dd.model.RData'))

str(dd.model)

## ************************************************************
## Configure the models files
## ************************************************************

source('analyses/occupancy-phylo/src/initialize.R')


dd.model$nsp<-dim(dd.model$X)[4]
dd.model$nsite<-dim(dd.model$X)[1]
dd.model$nrep<-array(6, dim=c(dim(dd.model$X)[1], dim(dd.model$X)[2] ,dim(dd.model$X)[4]))
dd.model$nyr<-dim(dd.model$X)[2]
dd.model$Z<-dd.model$z.init
dd.model$trait<-dd.model$sp.traits$average.mass #Trait holder since it won't run without it


## *******************************************************************
## Run the models --- Option 1. With detection [No detection at bottom]
## *******************************************************************


## load model file
## source('analyses/occupancy-phylo/src/models/jags_lambda.R')

## source('analyses/occupancy-phylo/src/models/det_covar/jags_lambda_det_cont_fact.R')

source('~/Dropbox/ccb_banding/analyses/occupancy-phylo/src/models/det_covar/jags_lambda_det_mu_only.R')

analyse <- function(dd, ...) {

  z.init <- dd$Z
  my.inits <- function() {
    list(Z=z.init)
  }

  dd <- list(data=dd,
             inits=my.inits,
             params=c(get.params(), 'psi.beta.sp')
             )
  
  run.R2jags.model(dd, ...)
}


## Run first model in which lambda is free to vary
res <- analyse(dd.model, ni=6000, nb=3000) 
res
logit.lambda <- res$BUGSoutput$sims.matrix[,'logit.lambda']
hist(expit(logit.lambda))
## traceplot(res)


res.full<-res
save(res.full, dd.model,file="~/Dropbox/ccb_banding/analyses/occupancy-phylo/saved/RJMCMC/Chase/POM_mods_ChaseData_nsp_123_det_mu_p.rdata")


## These parameters are saved to input into the second model to help
## with model mixing.
MU.lambda <- mean(logit.lambda)
TAU.lambda <- 1/sd(logit.lambda)^2

dd.model$MU.lambda <- MU.lambda 
dd.model$TAU.lambda <- TAU.lambda

source('~/Dropbox/ccb_banding/analyses/occupancy-phylo/src/models/jags_lambda_select.R')
source('~/Dropbox/ccb_banding/analyses/occupancy-phylo/src/models/det_covar/jags_lambda_select_det_cont_fact.R')


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

res.full<-res
res.select<-res2
res.psi.betas<-res
save(res.full, res.select, res.psi.betas, dd.model,file="~/Dropbox/ccb_banding/analyses/occupancy-phylo/saved/RJMCMC/Chase/POM_mods_ChaseData_nsp_123_det_mu_p.rdata")



##****************************************************************************************
## Load model and make figure (Still in option 1)
##****************************************************************************************
setwd('~/Dropbox/ccb_banding')
rm(list=ls())
source('analyses/multi-species-phylo/src/initialize.R')
source('analyses/occupancy-phylo/src/initialize.R')
## ************************************************************
load('~/Dropbox/ccb_banding/analyses/occupancy-phylo/saved/RJMCMC/Chase/POM_mods_ChaseData_nsp_123.rdata', verbose=T)

sumz<-res.psi.betas$BUGSoutput$summary
show<-c('mean', '2.5%', '97.5%')
sumz2<-round(sumz[,show], 2)
sumz2<-data.frame(sumz2)
sumz2$BCI<-paste("(", sumz2[,2], ", ",sumz2[,3], ")", sep="")
sumz3<-sumz2[,-c(2,3)]

mean(expit(sumz[grep("psi.0.sp", row.names(sumz)),])[,1])
sd(expit(sumz[grep("psi.0.sp", row.names(sumz)),])[,1])

### write.csv(sumz3, file='~/Desktop/ChaseParamEst.csv')

### Coef Fig.


nm<-read.csv("~/Dropbox/ccb_banding/data/original/NameMatch_Master_WorkingVersion.csv")

treez<-read.tree("~/Dropbox/ccb_banding/data/original/CRINDTreeALL.tre")

tree<-treez[[1]]


meanz<-pom.coefs(dd.model, res, phy=tree1,parameter.group='psi.beta.sp', plot=T, xlim=c(-20,20), tree.label.cex=0.3, color.gradient=c("red3", 'red' , "gold", "green3", 'darkgreen',"blue"))


meanz<-pom.coefs(dd.model, res, phy=treez[[5]],parameter.group='psi.beta.sp', plot=T, xlim=c(-20,20), tree.label.cex=0.3, cex.tip.labs=0.7,color.gradient=c("yellow",'gold', 'yellow3', 'green3','blue','darkblue'), xlab="Forest Affiliation (Slope coefficient; 100m buffer)")





##****************************************************************************************
## Option 2. Final model in which there is no detection (Naive model) 
##****************************************************************************************

source('~/Dropbox/ccb_banding/analyses/occupancy-phylo/src/models/det_covar/jags_lambda_no_det.R')


analyse <- function(dd, ...) {

  dd <- list(data=dd,
             params=c(get.params(), 'psi.beta.sp')
             )
  
  run.R2jags.model(dd, ...)
}


res <- analyse(dd.model, ni=6000, nb=3000, nt=10) 
res
logit.lambda <- res$BUGSoutput$sims.matrix[,'logit.lambda']
hist(expit(logit.lambda))

expit(res$BUGSoutput$summary['logit.lambda',c('mean','2.5%','97.5%')])

save(res, file="~/Dropbox/ccb_banding/analyses/occupancy-phylo/saved/RJMCMC/Chase/POM_mods_ChaseData_nsp_123_no_det.rdata")



params.naive<-res$BUGSoutput$summary[,c('mean','2.5%','97.5%')]

load('~/Dropbox/ccb_banding/analyses/occupancy-phylo/saved/RJMCMC/Chase/POM_mods_ChaseData_nsp_123_det_covar.rdata', verbose=T)

params.POM<-res.full$BUGSoutput$summary[,c('mean','2.5%','97.5%')]


params.naive[,'mean']
params.POM[,'mean']

params.POM.out<-params.POM[row.names(params.POM) %in% row.names(params.naive),]

round(data.frame(params.POM.out[,'mean'], params.naive[,'mean']), 3)

