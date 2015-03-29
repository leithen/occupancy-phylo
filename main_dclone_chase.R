## ************************************************************
## prepare matrices
## ************************************************************
setwd('~/Dropbox/ccb_banding')
rm(list=ls())
source('analyses/multi-species-phylo/src/initialize.R')
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

save(dd.model, file=file.path(path, 'occupancy/dd.model.RData'))

str(dd.model)

## ************************************************************
## Run the models
## ************************************************************

source('analyses/occupancy-phylo/src/initialize_temp.R')


dd.model$nsp<-dim(dd.model$X)[4]
dd.model$nsite<-dim(dd.model$X)[1]
dd.model$nrep<-array(6, dim=c(dim(dd.model$X)[1], dim(dd.model$X)[2] ,dim(dd.model$X)[4]))
dd.model$nyr<-dim(dd.model$X)[2]


lmer.df<-function(dd.model){
	df<-data.frame(
  			presence=as.vector(dd.model$X),
  			expand.grid(site=1:dd.model$nsite, year=1:dd.model$nyr, rep=1:dd.model$nrep[1], spp=1:dd.model$nsp),
  			env=rep(dd.model$env, dd.model$nyr*dd.model$nrep[1]*dd.model$nsp)
  			)
  			return(df)
}

dd.model$inputs$df<-lmer.df(dd.model)

dd.model$Z<-dd.model$z.init
dd.model$trait<-rep(0, dd.model$nsp)
nsp<-dd.model$nsp #update.priors lambda apparently needs nsp in global environment? Correct this later





nc<-c(2,5,20)
n.iter<-5000
n.adapt<-2000
n.update<-5000
num.cores<-3
draws<-1e3 #For importance sampling

##Run dclone for MLEs
dclone.out.lam <- run.dc.fit.update(dd.model=dd.model,
                                    nc=nc,
                                    n.iter=n.iter,
                                    n.adapt=n.adapt,
                                    n.update=n.update,
                                    num.cores=num.cores,
                                    model.case='phylo.lam')

#summary(dclone.out.lam)

file.ext<-paste("ChaseData_nsp",dd.model$nsp,sep="_")
lr.save.name.lam<-paste("analyses/occupancy-phylo/saved/dclone/POM_mods_",file.ext,"_lam_only.rdata", sep="")

save(dclone.out.lam, file= lr.save.name.lam)


##Given MLEs calculate random effect distributions
mod.lam <- run.jags.after.dclone(dclone.mod= dclone.out.lam)
## Check for convergence
max(attr(mod.lam,"jags.mod")$bugs$BUGSoutput$summary[,'Rhat'])

##Based on random effect distributions calculate LL
mod.lam<-jags.LL(mod.lam, draw.num=draws, num.cores=4, detection=T)

##Check statistics summary of the model
attr(mod.lam, "stat.sum")

############################################################
## No phylogenetic signal
############################################################

dclone.out.null <- run.dc.fit.update(dd.model=dd.model,
                         nc=nc,
                         n.iter=n.iter,
                         n.adapt=n.adapt,
                         n.update=n.update,
                         num.cores=num.cores,
                         model.case='phylo.null')

file.ext<-paste("ChaseData_nsp",dd.model$nsp,sep="_")
lr.save.name.null<-paste("analyses/occupancy-phylo/saved/dclone/POM_mods_",file.ext,"_null_only.rdata", sep="")

save(dclone.out.null, file= lr.save.name.null)
                                                  
mod.null <- run.jags.after.dclone(dclone.mod= dclone.out.null)

## Check for convergence
max(attr(mod.null,"jags.mod")$bugs$BUGSoutput$summary[,'Rhat'])

##Calc LLs
mod.null<-jags.LL(mod.null, draw.num= draws, num.cores=4, detection=T)

anova.pom(mod.lam, mod.null)

imp.samp.plot(mod.lam, mod.null)

attr(mod.lam, 'stat.sum')
attr(mod.null, 'stat.sum')

summary(mod.lam)
summary(mod.null)


file.ext<-paste("ChaseData_nsp",dd.model$nsp,sep="_")

lr.save.name<-paste("analyses/occupancy-phylo/saved/dclone/POM_mods_",file.ext,".rdata", sep="")
save(mod.lam, mod.null, file=lr.save.name)
