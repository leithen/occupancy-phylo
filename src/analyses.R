## a dual lambda (set to 0 and varying) runs with dc.fit
run.dc.fit.update <- function(dd.model,
                              nc,
                              n.iter,
                              n.adapt,
                              n.update,
                              num.cores,
                              model.case=c('phylo.lam','phylo.null','phylo.ou')) {

	raw.data<-dd.model #save raw data to append to model at the end

  ## load model file
  source(file.path('analyses/occupancy-phylo/src/models',
                   'clone_model_prior_update.R'))
  
  ## inits stuff
  require(lme4)
  df <- dd.model$inputs$df	
  lme4.out <- glmer(presence ~ env + (1|site) + (1|year) + (env|spp),
                    data=df, family='binomial')

  prm.est <- list(
    mu.psi.0=fixef(lme4.out)[[1]],
    sigma.psi.0 = attr(VarCorr(lme4.out)$spp, 'stddev')[1],
    sigma.psi.beta = attr(VarCorr(lme4.out)$spp, 'stddev')[1],
    sigma.site = attr(VarCorr(lme4.out)$site, 'stddev')[1],
    sigma.year = attr(VarCorr(lme4.out)$year, 'stddev')[1],
    mu.psi.beta=fixef(lme4.out)[[2]],
    psi.site= coef(lme4.out)$site[,1] - fixef(lme4.out)[1] ,
    psi.year= coef(lme4.out)$year[,1] - fixef(lme4.out)[1] ,
    psi.0.sp= coef(lme4.out)$spp[,1], # Because each spp has its own intercept no need to scale by the fixed effect
    psi.beta.sp= coef(lme4.out)$spp[,2])

  Z.init <- dd.model$Z
  Z <- dclone(dcdim(Z.init), n.clone=nc[1])
  psi.site.init <- dclone(dcdim(prm.est$psi.site), n.clone=nc[1])
  psi.year.init <- dclone(dcdim(prm.est$psi.year), n.clone=nc[1])
  psi.beta.sp.init <- dclone(dcdim(prm.est$psi.beta.sp), n.clone=nc[1])
  psi.0.sp.init <- prm.est$psi.0.sp ###This is a fixed effect so no cloning
  
  my.inits <- function() {
    list(Z=Z,
         psi.site=psi.site.init,
         psi.year=psi.year.init,
         psi.0.sp=psi.0.sp.init,
         psi.beta.sp=psi.beta.sp.init
         )
  }
  ## below is only used in dc.fit
  ## THIS HAS BEEN MODIFIED TO ACCEPT MEAN VALUES OF PREVIOUS RUNS
  ## If you put an init updater on lambda remember that it has to be
  ## for lambda.pre!!!
  ifun <- function(x, n.clones) {
  	
    Z <- dclone(dcdim(Z.init), n.clone=n.clones)
    psi.site.init <- dclone(dcdim(prm.est$psi.site), n.clone=n.clones)
    psi.year.init <- dclone(dcdim(prm.est$psi.year), n.clone=n.clones)
    psi.beta.sp.init <- dclone(dcdim(prm.est$psi.beta.sp), n.clone= n.clones)
    psi.0.sp.init <- prm.est$psi.0.sp ###This is a fixed effect so no cloning		
	
    list(Z=Z,
         psi.site=psi.site.init,
         psi.year=psi.year.init,		
         psi.0.sp=psi.0.sp.init,
         psi.beta.sp=psi.beta.sp.init,
         mu.psi.beta=coef(x)['mu.psi.beta'],
         mu.p.0=coef(x)['mu.p.0'],         
         sigma.psi.beta=coef(x)['sigma.psi.beta'],
         sigma.psi.site=coef(x)['sigma.psi.site'],
         sigma.psi.year=coef(x)['sigma.psi.year'],
         sigma.p.0=coef(x)['sigma.p.0']
         )
  }
  
  keep <- c('X', 'VCOV', 'ID',
            'nsite', 'nyr', 'nrep', 'nsp',
            'env', 'trait')
  dd.model <- dd.model[names(dd.model) %in% keep]

  dd.model$clones.index <- 1 ## dummy variable to scale by nc_i
  dd.model$X <- dcdim(dd.model$X)

  unchanged <- c('VCOV', 'ID', 'env', 'trait',
                 'nsite', 'nyr', 'nrep', 'nsp')
 
 	nsp<-dd.model$nsp   #Get nsp in easily accesible place for
 						#downstream functions
 
  ## ------------------------
	
  if (model.case=="phylo.lam") {
	  dd.model$pr <- update.priors.lam.free() ## Include initial prior values. 
	    									##Update priors function
  											## housed within main jags model file


  cl <- makePSOCKcluster(num.cores)
  ##   clusterExport(cl=cl, list="Z.init")

  phylo.lam <- dc.parfit(cl,
                        data=dd.model,
                        params=get.params(),
                        model=run.dclones.model.lam.free,
                        n.clones=nc,
                        multiply='clones.index',
                        n.iter=n.iter,
                        n.adapt=n.adapt,
                        n.update=n.update,
                        thin=10,
                        unchanged=unchanged,
                        inits=my.inits,
                        initsfun=ifun,
                        update="pr", 	### Remember to put this in quotes
                        updatefun=update.priors.lam.free,
                        partype='parchains')
	
	attr(phylo.lam, "model.case")<-model.case
	attr(phylo.lam, "data")<-raw.data

	stopCluster(cl)
	return(phylo.lam)
	}
	
	###For a single peak ornstein-ulhenbeck model. Work in progress.
  if (model.case=="phylo.ou") {
	  dd.model$pr <- update.priors.ou.free() ## Include initial prior values. 
	    									##Update priors function
  											## housed within main jags model file


  cl <- makePSOCKcluster(num.cores)
  ##   clusterExport(cl=cl, list="Z.init")

  phylo.ou <- dc.parfit(cl,
                        data=dd.model,
                        params=get.params.ou(),
                        model=run.dclones.model.ou.free,
                        n.clones=nc,
                        multiply='clones.index',
                        n.iter=n.iter,
                        n.adapt=n.adapt,
                        n.update=n.update,
                        thin=10,
                        unchanged=unchanged,
                        inits=my.inits,
                        initsfun=ifun,
#                        update="pr", 	### Remember to put this in quotes
#                        updatefun=update.priors.ou.free,
                        partype='parchains')
	
	attr(phylo.ou, "model.case")<-model.case
	attr(phylo.ou, "data")<-dd.model
	stopCluster(cl)
	
	return(phylo.ou)
	}	

if (model.case=="phylo.null") {
  dd.model$pr <- update.priors.lam.0() 	### Since lambda won't be tracked 
  ## need to eliminate it from the 
  ## prior inputs.

  cl <- makePSOCKcluster(num.cores)

  ##    clusterExport(cl=cl, list=list("Z.init"))

  phylo.null <- dc.parfit(cl,
                     data=dd.model,
                     params=get.params()[!get.params()=='lambda'],
                     model=run.dclones.model.lam.0,
                     n.clones=nc,
                     multiply='clones.index',
                     n.iter=n.iter,
                     n.adapt=n.adapt,
                     n.update=n.update,
                     thin=10,
                     unchanged=unchanged,
                     inits=my.inits,
                     initsfun=ifun,
                     update="pr", 	### Remember to put this in quotes
                     updatefun=update.priors.lam.0,
                     partype='parchains')
  
  attr(phylo.null, "model.case")<-model.case
  attr(phylo.null, "data")<-dd.model

  stopCluster(cl)
  return(phylo.null)
  }
  
  if (model.case=="phylo.lam.no.det") { ###MUCH OF THE FOLLOWING WILL BE REPEATED
  										## FROM ABOVE. THEREBY OVERWRITING OBJECTS

  
  my.inits <- function() {
    list(psi.site=psi.site.init,
         psi.year=psi.year.init,
         psi.0.sp=psi.0.sp.init,
         psi.beta.sp=psi.beta.sp.init
         )
  }
  ## below is only used in dc.fit
  ## THIS HAS BEEN MODIFIED TO ACCEPT MEAN VALUES OF PREVIOUS RUNS
  ## If you put an init updater on lambda remember that it has to be
  ## for lambda.pre!!!
  ifun <- function(x, n.clones) {
  	
    psi.site.init <- dclone(dcdim(prm.est$psi.site), n.clone=n.clones)
    psi.year.init <- dclone(dcdim(prm.est$psi.year), n.clone=n.clones)
    psi.beta.sp.init <- dclone(dcdim(prm.est$psi.beta.sp), n.clone= n.clones)
    psi.0.sp.init <- prm.est$psi.0.sp ###This is a fixed effect so no cloning		
	
    list(psi.site=psi.site.init,
         psi.year=psi.year.init,		
         psi.0.sp=psi.0.sp.init,
         psi.beta.sp=psi.beta.sp.init,
         mu.psi.beta=coef(x)['mu.psi.beta'],       
         sigma.psi.beta=coef(x)['sigma.psi.beta'],
         sigma.psi.site=coef(x)['sigma.psi.site'],
         sigma.psi.year=coef(x)['sigma.psi.year']
         )
  }
  
  keep <- c('X', 'VCOV', 'ID',
            'nsite', 'nyr', 'nrep', 'nsp',
            'env', 'trait')
  dd.model <- dd.model[names(dd.model) %in% keep]

  dd.model$clones.index <- 1 ## dummy variable to scale by nc_i
  dd.model$X <- dcdim(dd.model$X)
  dd.model$pr <- update.priors.no.det() ## Include initial prior values. 

  unchanged <- c('VCOV', 'ID', 'env', 'trait',
                 'nsite', 'nyr', 'nrep', 'nsp')
  ## ------------------------

  cl <- makePSOCKcluster(num.cores)
  ##   clusterExport(cl=cl, list="Z.init")

  no.det <- dc.parfit(cl,
                        data=dd.model,
                        params= get.params.no.det(),
                        model=run.dclones.model.no.det,
                        n.clones=nc,
                        multiply='clones.index',
                        n.iter=n.iter,
                        n.adapt=n.adapt,
                        n.update=n.update,
                        thin=10,
                        unchanged=unchanged,
                        inits=my.inits,
                        initsfun=ifun,
                        update="pr", 	### Remember to put this in quotes
                        updatefun=update.priors.no.det,
                        partype='parchains')
	
	
	attr(no.det, "model.case")<-model.case
	attr(no.det, "data")<-dd.model	
	stopCluster(cl)
	return(no.det)

	}  
  
  }




## **********************************************************************
## a single run with dc.fit
run.dc.fit <- function(dd.model,
                       nc,
                       n.iter,
                       n.adapt,
                       n.update,
                       num.cores) {

  ## load model file
  source(file.path('analyses/occupancy-phylo/src/models',
                   'clone_model.R'))
  
  ## inits stuff
  require(lme4)
  df <- dd.model$inputs$df	
  lme4.out <- glmer(presence~env+(1|site)+(1|year)+(env|spp), data=df, family='binomial')

  prm.est <- list(
    mu.psi.0=fixef(lme4.out)[[1]],
    sigma.psi.0 = attr(VarCorr(lme4.out)$spp, 'stddev')[1],
    sigma.psi.beta = attr(VarCorr(lme4.out)$spp, 'stddev')[1],
    sigma.site = attr(VarCorr(lme4.out)$site, 'stddev')[1],
    sigma.year = attr(VarCorr(lme4.out)$year, 'stddev')[1],
    mu.psi.beta=fixef(lme4.out)[[2]],
    psi.site= coef(lme4.out)$site[,1] - fixef(lme4.out)[1] ,
    psi.year= coef(lme4.out)$year[,1] - fixef(lme4.out)[1] ,
    psi.0.sp= coef(lme4.out)$spp[,1], #Because each spp has its own intercept no need 												for scale by the fixed effect
    psi.beta.sp= coef(lme4.out)$spp[,2]
    )


  Z.init <- dd.model$Z
  Z <- dclone(dcdim(Z.init), n.clone=nc[1])
  psi.site.init <- dclone(dcdim(prm.est$psi.site), n.clone=nc[1])
  my.inits <- function() {
    list(Z=Z
         ##	,psi.site=psi.site.init
         )
  }
  ## below is only used in dc.fit
  ifun <- function(x, n.clones) {
    Z <- dclone(dcdim(Z.init), n.clone=n.clones)
    psi.site.init <- dclone(dcdim(prm.est$psi.site), n.clone=n.clones)
    list(Z=Z
         ##	,psi.site=psi.site.init
         )
  }
  
  keep <- c('X', 'VCOV', 'ID',
            'nsite', 'nyr', 'nrep', 'nsp',
            'env', 'trait')
  dd.model <- dd.model[names(dd.model) %in% keep]

  ## ## --- for full version ---
  ## dd.model$clones.index <- 1 ## dummy variable to scale by nc_i
  ## dd.model$X <- dcdim(dd.model$X)
  ## dd.model$VCOV <- dcdim(dd.model$VCOV)
  ## dd.model$ID <- dcdim(dd.model$ID)
  ## dd.model$env <- dcdim(dd.model$env)
  ## dd.model$trait <- dcdim(dd.model$trait)

  ## unchanged <- c('nsite', 'nyr', 'nrep', 'nsp')
  ## ## ------------------------
  
  ## --- for fewer clones ---
  dd.model$clones.index <- 1 ## dummy variable to scale by nc_i
  dd.model$X <- dcdim(dd.model$X)

  unchanged <- c('VCOV', 'ID', 'env', 'trait',
                 'nsite', 'nyr', 'nrep', 'nsp')
  ## ------------------------

  ## browser()
  if(num.cores==1)
    res <- dc.fit(data=dd.model,
                  params=get.params(),
                  model=run.dclones.model,
                  n.clones=nc,
                  multiply='clones.index',
                  n.iter=n.iter,
                  n.adapt=n.adapt,
                  n.update=n.update,
                  thin=20,
                  unchanged=unchanged,
                  inits=my.inits,
                  initsfun=ifun)
  if(num.cores>1) {
    cl <- makePSOCKcluster(num.cores)
    res <- dc.parfit(cl,
                     data=dd.model,
                     params=get.params(),
                     model=run.dclones.model,
                     n.clones=nc,
                     multiply='clones.index',
                     n.iter=n.iter,
                     n.adapt=n.adapt,
                     n.update=n.update,
                     thin=20,
                     unchanged=unchanged,
                     inits=my.inits,
                     initsfun=ifun,
                     partype='balancing')
  }
  res
}


## a single run with jags.fit
run.jags.fit <- function(dd.model,
                         nc,
                         n.iter,
                         n.adapt,
                         n.update) {

  ## inits stuff
  Z.init <- dd.model$Z
  Z <- dclone(dcdim(Z.init), n.clone=nc)
  my.inits <- function() {
    list(Z=Z)
  }
  
  keep <- c('X', 'VCOV', 'ID',
            'nsite', 'nyr', 'nrep', 'nsp',
            'env', 'trait')
  dd.model <- dd.model[names(dd.model) %in% keep]

  dd.model$clones.index <- 1 ## dummy variable to scale by nc
  dd.model$X <- dcdim(dd.model$X)
  dd.model$VCOV <- dcdim(dd.model$VCOV)
  dd.model$ID <- dcdim(dd.model$ID)
  dd.model$env <- dcdim(dd.model$env)
  dd.model$trait <- dcdim(dd.model$trait)
  
  dd.model.dclone <- dclone(dd.model,
                            n.clones=nc,
                            unchanged=c('nsite', 'nyr', 'nrep', 'nsp'),
                            multiply='clones.index')
  
  jags.fit(data=dd.model.dclone,
           params=get.params(),
           inits=my.inits,
           model=run.dclones.model,
           n.iter=n.iter,
           n.adapt=n.adapt,
           n.update=n.update,
           thin=20)
  
}

## JAGS analysis for post declone draws
run.jags.after.dclone <- function(dclone.mod,
                                  ni=5000,
                                  nt=5,
                                  nb=1000,
                                  nc=3,
                                  ALTERNATIVE_DELETE_ME=F) { ## ADDED THIS ARGUMENT 
                                  							 ## TO ALLOW LOADING A MODEL
                                  							 ## IN WHICH Z STATES ARE
                                  							 ## NOT ESTIMATED

	model.case<-attr(dclone.mod, "model.case")
	dd<-attr(dclone.mod, "data")


  ##Extract parameter estimates from dclone.model
  ## Extract MLEs for lambda model

  if (model.case=="phylo.lam"){
summ.lam.free <- summary(dclone.mod)[[1]]
prms <- list(mu.p.0= summ.lam.free['mu.p.0','Mean'],
             sigma.p.0= summ.lam.free['sigma.p.0','Mean'],
             sigma.psi.site= summ.lam.free['sigma.psi.site','Mean'],
             sigma.psi.year= summ.lam.free['sigma.psi.year','Mean'],
             psi.0.sp= summ.lam.free[sapply(1:dd$nsp,
               function(s) sprintf('psi.0.sp[%s]', s)),'Mean'],
             mu.psi.beta= summ.lam.free['mu.psi.beta','Mean'],
             sigma.psi.beta= summ.lam.free['sigma.psi.beta','Mean'],
             lambda= summ.lam.free['lambda','Mean'])
             }

  if (model.case=="phylo.ou"){
summ.lam.free <- summary(dclone.mod)[[1]]
prms <- list(mu.p.0= summ.lam.free['mu.p.0','Mean'],
             sigma.p.0= summ.lam.free['sigma.p.0','Mean'],
             sigma.psi.site= summ.lam.free['sigma.psi.site','Mean'],
             sigma.psi.year= summ.lam.free['sigma.psi.year','Mean'],
             psi.0.sp= summ.lam.free[sapply(1:dd$nsp,
               function(s) sprintf('psi.0.sp[%s]', s)),'Mean'],
             mu.psi.beta= summ.lam.free['mu.psi.beta','Mean'],
             sigma.psi.beta= summ.lam.free['sigma.psi.beta','Mean'],
             alpha= summ.lam.free['alpha','Mean'])
             }


  if (model.case=="phylo.null") {
## Extract MLEs for lambda = 0 model
summ.lam.0 <- summary(dclone.mod)[[1]]
prms <- list(mu.p.0=summ.lam.0['mu.p.0','Mean'],
             sigma.p.0=summ.lam.0['sigma.p.0','Mean'],
             sigma.psi.site=summ.lam.0['sigma.psi.site','Mean'],
             sigma.psi.year=summ.lam.0['sigma.psi.year','Mean'],
             psi.0.sp=summ.lam.0[sapply(1:dd$nsp,
               function(s) sprintf('psi.0.sp[%s]', s)),'Mean'],
             mu.psi.beta=summ.lam.0['mu.psi.beta','Mean'],
             sigma.psi.beta=summ.lam.0['sigma.psi.beta','Mean'],
             lambda= 0 )
             }
             
  if (model.case=='phylo.lam.no.det') {             
## Extract MLEs for no detection model
summ.no.det <- summary(dclone.mod)[[1]]
prms<- list(mu.p.0=999,  ###CORECT THIS
             sigma.p.0=0.001,	###CORECT THIS
             sigma.psi.site= summ.no.det['sigma.psi.site','Mean'],
             sigma.psi.year= summ.no.det['sigma.psi.year','Mean'],
             psi.0.sp= summ.no.det[sapply(1:dd$nsp,
               function(s) sprintf('psi.0.sp[%s]', s)),'Mean'],
             mu.psi.beta= summ.no.det['mu.psi.beta','Mean'],
             sigma.psi.beta= summ.no.det['sigma.psi.beta','Mean'],
             lambda= summ.no.det['lambda','Mean'] )
             }




  z.init <- apply(dd$X, c(1,2,4),
                  function(x) (sum(x,na.rm=TRUE)>0)*1)
  z.init[apply(dd$X, c(1,2,4),
               function(x) !any(!is.na(x)))] <- NA
  
  
  if(model.case=='phylo.lam.no.det'){
  	  d <- list(X=dd$X,
            env=dd$env,
            trait=dd$trait,
            nsite=dim(dd$X)[1],
            nrep=dd$nrep,
            nsp=dim(dd$X)[4],
            nyr=dim(dd$X)[2],
            sigma.psi.site=prms$sigma.psi.site,
            sigma.psi.year=prms$sigma.psi.year,
            psi.0.sp=prms$psi.0.sp,
            mu.psi.beta=prms$mu.psi.beta,
            sigma.psi.beta=prms$sigma.psi.beta,
            lambda=prms$lambda,
            VCOV=dd$VCOV,
            ID=dd$ID)
  }
  
  
  
  if(model.case=='phylo.lam'|model.case=='phylo.null') {

  d <- list(X=dd$X,
            env=dd$env,
            trait=dd$trait,
            nsite=dim(dd$X)[1],
            nrep=dd$nrep,
            nsp=dim(dd$X)[4],
            nyr=dim(dd$X)[2],
            mu.p.0=prms$mu.p.0,
            sigma.p.0=prms$sigma.p.0,
            sigma.psi.site=prms$sigma.psi.site,
            sigma.psi.year=prms$sigma.psi.year,
            psi.0.sp=prms$psi.0.sp,
            mu.psi.beta=prms$mu.psi.beta,
            sigma.psi.beta=prms$sigma.psi.beta,
            lambda=prms$lambda,
            VCOV=dd$VCOV,
            ID=dd$ID)
			}
			
  if(model.case=='phylo.ou') {
  d <- list(X=dd$X,
            env=dd$env,
            trait=dd$trait,
            nsite=dim(dd$X)[1],
            nrep=dd$nrep,
            nsp=dim(dd$X)[4],
            nyr=dim(dd$X)[2],
            mu.p.0=prms$mu.p.0,
            sigma.p.0=prms$sigma.p.0,
            sigma.psi.site=prms$sigma.psi.site,
            sigma.psi.year=prms$sigma.psi.year,
            psi.0.sp=prms$psi.0.sp,
            mu.psi.beta=prms$mu.psi.beta,
            sigma.psi.beta=prms$sigma.psi.beta,
            alpha=prms$alpha,
            VCOV=dd$VCOV)
			}
						


  if (model.case=='phylo.lam.no.det'){
  	 source(file.path('analyses/occupancy-phylo/src/models',
                   'post_dclone_no_det.R'))
	  
	  dd <- list(data=d, params=get.params())

  	  res <- list(data=d, bugs=run.model.no.det(dd, ni, nt, nb, nc))
  	}
 
  if(model.case=='phylo.lam'|model.case=='phylo.null') {
	source(file.path('analyses/occupancy-phylo/src/models',
                 'post_dclone.R'))



	if (ALTERNATIVE_DELETE_ME==TRUE) {source(file.path('analyses/occupancy-phylo/src/models', ## DELETE
                 'post_dclone_VALIDATING.R'))}												  ## THESE LINES



  	  ## create inits
  		my.inits <- function() {
    			list(Z=z.init)
  			}

		## data
  		dd <- list(data=d, inits=my.inits, params=get.params())
  		res <- list(data=d, bugs=run.model(dd, ni, nt, nb, nc)) 
  		}
  		
   if(model.case=='phylo.ou') {
  	 source(file.path('analyses/occupancy-phylo/src/models',
                   'post_dclone_ou.R'))  	  
  	  
  	  ## create inits
  		my.inits <- function() {
    			list(Z=z.init)
  			}

		## data
  		dd <- list(data=d, inits=my.inits, params=get.params())
  		res <- list(data=d, bugs=run.model(dd, ni, nt, nb, nc)) 
  		} 		
  		
  
  	attr(dclone.mod, "jags.mod")<-res
   
   return(dclone.mod)
  
}
## ************************************************************

### LoF. 22-Dec-2014: Creating master LL calculation function that streamlines code.

jags.LL<-function(dclone.mod, draw.num, num.cores, detection=TRUE) {
  
require(mvtnorm)
expit <- function(x) 1/(1+exp(-x))

model.case<-attr(dclone.mod, "model.case")
jags.mod<-attr(dclone.mod, "jags.mod")

##Set prms
prms<-jags.mod$data

## Extract raw chains
chains <- jags.mod$bugs$BUGSoutput$sims.matrix

## drop deviance
chains <- chains[,-which(colnames(chains)=='deviance')]

## mean vector
mu.vec <- colMeans(chains)
## variance-covariance matrix
vcov.mat <- cov(chains)	

random.draws <- mvrnorm(n=draw.num, mu=mu.vec, Sigma=vcov.mat)

##

if (detection==TRUE){
	vals <-as.numeric(mclapply(1:draw.num, par.calc.LL, random.draws=random.draws, prms=prms, mu.vec=mu.vec, vcov.mat=vcov.mat, model.case=model.case, mc.cores=num.cores))
}

if (detection==FALSE){
	vals <-as.numeric(mclapply(1:draw.num, par.calc.LL.no.det, random.draws=random.draws, prms=prms, mu.vec=mu.vec, vcov.mat=vcov.mat, mc.cores=num.cores))
}


	attr(dclone.mod, "lik.vals")<-vals

	mean.vals <- log(mean(exp(vals - mean(vals)))) + mean(vals)
	se.vals <- sqrt(var(vals-mean(vals)))/sqrt(length(vals))
	AIC<-2*length(coef(dclone.mod)) - 2* mean.vals
	out<-c(mean.vals, se.vals, AIC, draw.num)
	names(out)<-c("LL", "SE Log Lik.", "AIC", "Sample Number")
	out<-round(out, 5)
	
	attr(dclone.mod, "stat.sum")<-out

	return(dclone.mod)

}





par.calc.LL <- function(chain, random.draws, prms=prms, mu.vec=mu.vec, vcov.mat=vcov.mat, model.case=model.case) {

  ## get important data into easily accessible form
  env <- prms$env
  nsite <- prms $nsite
  nyr <- prms $nyr
  nrep <- max(prms$nrep) ## note: assuming same num of reps
  nsp <- prms $nsp
  X <- prms$X
  
  if(model.case=='phylo.lam'|model.case=='phylo.null'){
  lambda.mat <- prms$lambda* prms$VCOV +
    (1-prms$lambda)* prms $ID
	}

  ## get mle estimates into easily accessible form
  mu.p.0 <- prms$mu.p.0
  sigma.p.0 <- prms$sigma.p.0
  mu.psi.beta <- prms$mu.psi.beta
  sigma.psi.beta <- prms$sigma.psi.beta
  sigma.psi.site <- prms$sigma.psi.site
  sigma.psi.year <- prms$sigma.psi.year
  psi.0.sp <- prms$psi.0.sp
  
  if(model.case=='phylo.lam'|model.case=='phylo.null'){  
    lambda <- prms$lambda
  }
  
  if(model.case=='phylo.ou'){  
    VCOV<-prms$VCOV
    alpha <- prms$alpha
  }

  p.0.terms <- grep('p.0',         names(random.draws[1,]))
  psi.beta.terms <- grep('psi.beta.sp', names(random.draws[1,]))
  psi.site.terms <- grep('psi.site',    names(random.draws[1,]))
  psi.year.terms <- grep('psi.year',    names(random.draws[1,]))



  p.0.vec      <- random.draws[chain, p.0.terms]
  psi.0.vec <- psi.0.sp
  psi.beta.vec <- random.draws[chain, psi.beta.terms]
  psi.site.vec <- random.draws[chain, psi.site.terms]
  psi.year.vec <- random.draws[chain, psi.year.terms]

if (model.case=="phylo.lam"|model.case=="phylo.null"){
  RR.psi.beta <- dmvnorm(psi.beta.vec,
                         mean=rep(mu.psi.beta, nsp),
                         sigma=(sigma.psi.beta^2)*lambda.mat, 
                         log=TRUE)
      }

if (model.case=="phylo.ou"){
  RR.psi.beta <- dmvnorm(psi.beta.vec,
                         mean=rep(mu.psi.beta, nsp),
                         sigma=((sigma.psi.beta^2)/(2*alpha))*exp(-2*alpha*(1-VCOV))*(1-exp(-2*alpha*VCOV)),
                         log=TRUE)
      }

                                                  
  LL.psi.beta.given.theta <- sum(RR.psi.beta) 	#The summation is unnessesary as the dmvnorm
                                        # already is looking across the entire
                                        # psi.beta.vec

  RR.site <- dnorm(psi.site.vec,
                   mean=0,
                   sd=sigma.psi.site,
                   log=TRUE)
  LL.psi.site.given.theta <- sum(RR.site)

  RR.year <- dnorm(psi.year.vec,
                   mean=0,
                   sd=sigma.psi.year,
                   log=TRUE)
  LL.psi.year.given.theta <- sum(RR.year) 

  RR.p.0 <- dnorm(p.0.vec,
                  mean=mu.p.0,
                  sd=sigma.p.0,
                  log=TRUE)
  LL.p.given.theta <- sum(RR.p.0)

  ## Calculates probability of obtaining a certain occupancy status at a
  ## site, at a year, (and if this were actually on observations rather
  ## than occupancy, during a rep), given the expected psi for the
  ## species at the site.

  ## vector of species specific effects on occurrence
  psi.sp.int.vec <- rowSums(expand.grid(psi.site.vec,
                                        psi.year.vec,
                                        rep(0,nrep),
                                        psi.0.vec))
  ## vector of species-specific environmental slopes
  psi.sp.slope.vec <- apply(expand.grid(env,
                                        rep(1,nyr),
                                        rep(1,nrep),
                                        psi.beta.vec), 1, prod)
  ## occupancy matrix
  psi.mat <- expit(array(psi.sp.int.vec +
                         psi.sp.slope.vec, dim=dim(X)))

  ## Calculate the probability of obtianing the the detection history
  ## given the expected occupancy at a site, and the estimated
  ## detectability of each spp
  ##
  ## As long as species detectability does not depend on site or year
  ## this should work fine. Spp as final layer.
  
  ## detectability probs
  p.sp.vec <- rowSums(expand.grid(rep(0,nsite),
                                  rep(0,nyr),
                                  rep(0,nrep),
                                  p.0.vec))
  p.mat <- expit(array(p.sp.vec, dim=dim(X)))
  
  ## calculate probability of X given psi.mat and p.mat
  QQ <- log(psi.mat*dbinom(X, size=1, prob=p.mat, log=FALSE) +
            (1-psi.mat)*(X==0))
  LL.X.given.psi.and.p <- sum(QQ)

  ## Calculate probability of random effect parameters given the
  ## multi-variate normal we used to draw them from
  used <- c(p.0.terms,
            psi.beta.terms,
            psi.site.terms,
            psi.year.terms) 

  TT <- dmvnorm(random.draws[chain, used],
                mean=mu.vec[used],
                sigma=vcov.mat[used, used],
                log=TRUE)
  LL.theta.given.MU.and.SIGMA <- TT
  
  ## calculate the total log likelihood.
  LL.psi.beta.given.theta +
    LL.psi.site.given.theta +
      LL.psi.year.given.theta +
        LL.p.given.theta + 
          LL.X.given.psi.and.p -
            LL.theta.given.MU.and.SIGMA
  
}



par.calc.LL.no.det <- function(chain, random.draws, prms=prms, mu.vec=mu.vec, vcov.mat=vcov.mat) {
  ## get important data into easily accessible form
  env <- prms$env
  nsite <- prms $nsite
  nyr <- prms $nyr
  nrep <- max(prms$nrep) ## note: assuming same num of reps
  nsp <- prms $nsp
  X <- prms$X
  lambda.mat <- prms$lambda* prms$VCOV +
    (1-prms$lambda)* prms $ID

  ## get mle estimates into easily accessible form

  mu.psi.beta <- prms$mu.psi.beta
  sigma.psi.beta <- prms$sigma.psi.beta
  sigma.psi.site <- prms$sigma.psi.site
  sigma.psi.year <- prms$sigma.psi.year
  psi.0.sp <- prms$psi.0.sp
  lambda <- prms$lambda
  


  psi.beta.terms <- grep('psi.beta.sp', names(random.draws[1,]))
  psi.site.terms <- grep('psi.site',    names(random.draws[1,]))
  psi.year.terms <- grep('psi.year',    names(random.draws[1,]))

  psi.0.vec <- psi.0.sp
  psi.beta.vec <- random.draws[chain, psi.beta.terms]
  psi.site.vec <- random.draws[chain, psi.site.terms]
  psi.year.vec <- random.draws[chain, psi.year.terms]

  RR.psi.beta <- dmvnorm(psi.beta.vec,
                         mean=rep(mu.psi.beta, nsp),
                         sigma=(sigma.psi.beta^2)*lambda.mat, log=TRUE)

  LL.psi.beta.given.theta <- sum(RR.psi.beta) 	#The summation is unnessesary as the dmvnorm
                                        # already is looking across the entire
                                        # psi.beta.vec

  RR.site <- dnorm(psi.site.vec,
                   mean=0,
                   sd=sigma.psi.site,
                   log=TRUE)
  LL.psi.site.given.theta <- sum(RR.site)

  RR.year <- dnorm(psi.year.vec,
                   mean=0,
                   sd=sigma.psi.year,
                   log=TRUE)
  LL.psi.year.given.theta <- sum(RR.year) 


  ## Calculates probability of obtaining a certain occupancy status at a
  ## site, at a year, (and if this were actually on observations rather
  ## than occupancy, during a rep), given the expected psi for the
  ## species at the site.

  ## vector of species specific effects on occurrence
  psi.sp.int.vec <- rowSums(expand.grid(psi.site.vec,
                                        psi.year.vec,
                                        rep(0,nrep),
                                        psi.0.vec))
  ## vector of species-specific environmental slopes
  psi.sp.slope.vec <- apply(expand.grid(env,
                                        rep(1,nyr),
                                        rep(1,nrep),
                                        psi.beta.vec), 1, prod)
  ## occupancy matrix
  psi.mat <- expit(array(psi.sp.int.vec +
                         psi.sp.slope.vec, dim=dim(X)))


  
  ## I think this is correct for the case when we ignore detectability
  QQ <- dbinom(X, size=1, prob=psi.mat, log=T)
  LL.X.given.psi<- sum(QQ)

  ## Calculate probability of random effect parameters given the
  ## multi-variate normal we used to draw them from
  used <- c(
            psi.beta.terms,
            psi.site.terms,
            psi.year.terms) 


  TT <- dmvnorm(random.draws[chain, used],
                mean=mu.vec[used],
                sigma=vcov.mat[used, used],
                log=TRUE)
  LL.theta.given.MU.and.SIGMA <- TT
  
  ## calculate the total log likelihood.
  LL.psi.beta.given.theta +
    LL.psi.site.given.theta +
      LL.psi.year.given.theta +
          LL.X.given.psi -
            LL.theta.given.MU.and.SIGMA

  
}





anova.pom<-function(object, ...){
	  ancall <- sys.call()
      ancall$verbose <- ancall$test <- NULL

	  
	  aux <- list(object, ...)
	  deparse(substitute(aux[[1]]))
		rt<-length(aux)
				
		statz<-simplify2array(lapply(aux, attr, "stat.sum"))
		dfModel <- unlist(lapply(lapply(aux, coef), length))
		logLik<-statz[1,]
		AIC<-statz[3,]
		
		aod <- data.frame(df = dfModel, 
            AIC = AIC, logLik = logLik)
            
        ddf <- diff(dfModel)
            if (sum(abs(ddf)) > 0) {
                effects <- rep("", rt)
                for (i in 2:rt) {
                  if (ddf[i - 1] != 0) {
                    effects[i] <- paste(i - 1, i, sep = " vs ")
                  }
                }

        pval <- rep(NA, rt - 1)
        ldf <- as.logical(ddf)
        lratio <- 2 * abs(diff(logLik))
        lratio[!ldf] <- NA
        pval[ldf] <- 1 - pchisq(lratio[ldf], abs(ddf[ldf]))
        aod <- data.frame(aod, Test = effects, L.Ratio = c(NA, 
            lratio), `p-value` = c(NA, pval), check.names = FALSE)		
		}
        
        row.names(aod) <- unlist(lapply(as.list(ancall[-1L]), 
            deparse))
		
	aod
}


imp.samp.plot<-function(object ,...) {
	  require(scales)
	  aux <- list(object, ...)
	  ancall <- sys.call()
      ancall$verbose <- ancall$test <- NULL

	
	x.most<-max(unlist(lapply(1:length(aux), function(x) {length(attr(aux[[x]], 'lik.vals'))})))
	y.range<-range(unlist(lapply(1:length(aux), function(x) {(attr(aux[[x]], 'lik.vals'))})))
	colz<-rainbow(length(aux))


	plot(0,0, ylab="Log Likelihoods", xlab='Draw Number', 
		ylim=y.range, xlim=c(1,x.most), type='n')

		ot<-list()
		
	for(i in 1:length(aux)) {
		vals<-attr(aux[[i]], 'lik.vals')
		ot[[i]]<-numeric()
		notch<-length(vals)/1000
		
		for (j in 1:1000){
			ot[[i]][j]<-log(mean(exp(vals[1:(j*notch)]-mean(vals[1:(j*notch)]))))+ mean(vals[1:(j*notch)])
		}
	x<-seq(notch, length(vals) , by=notch)

	points(vals, col=alpha(colz[i], 0.2), pch=1)	
	}

for (i in 1:length(aux)){
	points(x, ot[[i]], type='l', col="black", lwd=6)
	points(x, ot[[i]], type='l', col=colz[i], lwd=3)
}

namez<-unlist(lapply(as.list(ancall[-1L]), 
            deparse))
legend('bottomleft', namez, lwd=3, col=colz)
}


