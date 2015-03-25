## a dual lambda (set to 0 and varying) runs with dc.fit
run.dc.fit.update <- function(dd.model,
                              nc,
                              n.iter,
                              n.adapt,
                              n.update,
                              num.cores) {
  ##	dd.model <- dd.simulated
  ##	n.iter=100
  ##	n.adapt=100
  ##	n.update=100
  ##	num.cores=3
  ##	nc=c(2,3)	


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
    sigma.psi.beta = attr(VarCorr(lme4.out)$spp, 'stddev')[2], #sigma.psi.beta is actually the variance term (because of how the data are simulated in the mvnorm). So don't use this in downstream estimates.
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
  dd.model$pr <- update.priors.lam.free() ## Include initial prior values. 
  ##Update priors function
  ## housed within main jags model file

  unchanged <- c('VCOV', 'ID', 'env', 'trait',
                 'nsite', 'nyr', 'nrep', 'nsp')
  ## ------------------------

  cl <- makePSOCKcluster(num.cores)
  ##   clusterExport(cl=cl, list="Z.init")

  lam.free <- dc.parfit(cl,
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

  dd.model$pr <- update.priors.lam.0() 	### Since lambda won't be tracked 
  ## need to eliminate it from the 
  ## prior inputs.

  cl <- makePSOCKcluster(num.cores)

  ##    clusterExport(cl=cl, list=list("Z.init"))

  lam.0 <- dc.parfit(cl,
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
  
  list(lam.free, lam.0)
}


## **********************************************************************
## **********************************************************************
run.dc.fit.no.det <- function(dd.model,
                              nc,
                              n.iter,
                              n.adapt,
                              n.update,
                              num.cores) {
  ##	dd.model <- dd.simulated
  ##	n.iter=100
  ##	n.adapt=100
  ##	n.update=100
  ##	num.cores=3
  ##	nc=c(2,3)	


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
    sigma.psi.beta = attr(VarCorr(lme4.out)$spp, 'stddev')[2], #sigma.psi.beta is actually the variance term (because of how the data are simulated in the mvnorm). So don't use this in downstream estimates.
    sigma.site = attr(VarCorr(lme4.out)$site, 'stddev')[1],
    sigma.year = attr(VarCorr(lme4.out)$year, 'stddev')[1],
    mu.psi.beta=fixef(lme4.out)[[2]],
    psi.site= coef(lme4.out)$site[,1] - fixef(lme4.out)[1] ,
    psi.year= coef(lme4.out)$year[,1] - fixef(lme4.out)[1] ,
    psi.0.sp= coef(lme4.out)$spp[,1], # Because each spp has its own intercept no need to scale by the fixed effect
    psi.beta.sp= coef(lme4.out)$spp[,2])

  psi.site.init <- dclone(dcdim(prm.est$psi.site), n.clone=nc[1])
  psi.year.init <- dclone(dcdim(prm.est$psi.year), n.clone=nc[1])
  psi.beta.sp.init <- dclone(dcdim(prm.est$psi.beta.sp), n.clone=nc[1])
  psi.0.sp.init <- prm.est$psi.0.sp ###This is a fixed effect so no cloning
  
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
  ##Update priors function
  ## housed within main jags model file

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
                        n.iter=100,    #n.iter,
                        n.adapt=100,   #n.adapt,
                        n.update=100,  #n.update,
                        thin=10,
                        unchanged=unchanged,
#                        inits=my.inits,
#                        initsfun=ifun,
#                        update="pr", 	### Remember to put this in quotes
#                        updatefun=update.priors.no.det,
                        partype='parchains')

	return(no.det)
}


## **********************************************************************

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
    sigma.psi.beta = attr(VarCorr(lme4.out)$spp, 'stddev')[2],
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
run.jags.after.dclone <- function(dd,
                                  mu.p.0,
                                  sigma.p.0,
                                  sigma.psi.site,
                                  sigma.psi.year,
                                  psi.0.sp,
                                  mu.psi.beta,
                                  sigma.psi.beta,
                                  lambda,
                                  VCOV,
                                  ID,
                                  ni,
                                  nt,
                                  nb,
                                  nc) {

  z.init <- apply(dd$X, c(1,2,4),
                  function(x) (sum(x,na.rm=TRUE)>0)*1)
  z.init[apply(dd$X, c(1,2,4),
               function(x) !any(!is.na(x)))] <- NA
  
  d <- list(X=dd$X,
            env=dd$env,
            trait=dd.simulated$trait,
            nsite=dim(dd$X)[1],
            nrep=dd$nrep,
            nsp=dim(dd$X)[4],
            nyr=dim(dd$X)[2],
            mu.p.0=mu.p.0,
            sigma.p.0=sigma.p.0,
            sigma.psi.site=sigma.psi.site,
            sigma.psi.year=sigma.psi.year,
            psi.0.sp=psi.0.sp,
            mu.psi.beta=mu.psi.beta,
            sigma.psi.beta=sigma.psi.beta,
            lambda=lambda,
            VCOV=VCOV,
            ID=ID)

  ## create inits
  my.inits <- function() {
    list(Z=z.init)
  }

  ## data
  dd <- list(data=d, inits=my.inits, params=get.params())

  res <- list(data=d, bugs=run.model(dd, ni, nt, nb, nc))
  res
  
}
## ************************************************************

 
## JAGS analysis for post declone draws for no detectability case
run.jags.after.dclone.no.det <- function(dd,
                                  sigma.psi.site,
                                  sigma.psi.year,
                                  psi.0.sp,
                                  mu.psi.beta,
                                  sigma.psi.beta,
                                  lambda,
                                  VCOV,
                                  ID,
                                  ni,
                                  nt,
                                  nb,
                                  nc) {

 source(file.path('analyses/occupancy-phylo/src/models',
                   'post_dclone_no_det.R'))


  d <- list(X=dd$X,
            env=dd$env,
            trait=dd.simulated$trait,
            nsite=dim(dd$X)[1],
            nrep=dd$nrep,
            nsp=dim(dd$X)[4],
            nyr=dim(dd$X)[2],
            sigma.psi.site=sigma.psi.site,
            sigma.psi.year=sigma.psi.year,
            psi.0.sp=psi.0.sp,
            mu.psi.beta=mu.psi.beta,
            sigma.psi.beta=sigma.psi.beta,
            lambda=lambda,
            VCOV=VCOV,
            ID=ID)


  ## data
  dd <- list(data=d, params=get.params())

  res <- list(data=d, bugs=run.model.no.det(dd, ni, nt, nb, nc))
  res
  
}
## ************************************************************


### LoF. 19-Dec-2014: Calculates log likelihood for each random 
## draw deriving from the jags random effect estimates after dclone
## determines MLEs of main model parameters. Input random.draws, 
## MLE parameter estimates (prms), and the relevant mean vector 
## and varaiance covariance matrix of all random effect parameters 
## from jags. This function can likely be streamlined to increase
## speed and efficiency

par.calc.LL <- function(chain, random.draws, prms, mu.vec, vcov.mat) {

  require(mvtnorm)
  expit <- function(x) 1/(1+exp(-x))

  ## get important data into easily accessible form
  env <- dd.simulated$env
  nsite <- dd.simulated$nsite
  nyr <- dd.simulated$nyr
  nrep <- max(dd.simulated$nrep) ## note: assuming same num of reps
  nsp <- dd.simulated$nsp
  X <- dd.simulated$X
  lambda.mat <- prms$lambda*dd.simulated$VCOV +
    (1-prms$lambda)*dd.simulated$ID

  ## get mle estimates into easily accessible form
  mu.p.0 <- prms$mu.p.0
  sigma.p.0 <- prms$sigma.p.0
  mu.psi.beta <- prms$mu.psi.beta
  sigma.psi.beta <- prms$sigma.psi.beta
  sigma.psi.site <- prms$sigma.psi.site
  sigma.psi.year <- prms$sigma.psi.year
  psi.0.sp <- prms$psi.0.sp
  lambda <- prms$lambda
  
  p.0.terms      <- grep('p.0',         names(random.draws[1,]))
  ## psi.0.terms    <- grep('psi.0',       names(random.draws[1,]))
  psi.beta.terms <- grep('psi.beta.sp', names(random.draws[1,]))
  psi.site.terms <- grep('psi.site',    names(random.draws[1,]))
  psi.year.terms <- grep('psi.year',    names(random.draws[1,]))

  p.0.vec      <- random.draws[chain, p.0.terms]
  ## psi.0.vec    <- random.draws[chain, psi.0.terms]
  psi.0.vec <- psi.0.sp
  psi.beta.vec <- random.draws[chain, psi.beta.terms]
  psi.site.vec <- random.draws[chain, psi.site.terms]
  psi.year.vec <- random.draws[chain, psi.year.terms]

  ## ## because we do not have a random effect of species, the below is
  ## ## no longer needed
  ## RR.psi.0 <- dnorm(psi.0.vec,
  ##                   mean=mu.psi.0,
  ##                   sd=sigma.psi.0,
  ##                   log=TRUE)
  ## LL.psi.0.given.theta <- sum(RR.psi.0)

  RR.psi.beta <- dmvnorm(psi.beta.vec,
                         mean=rep(mu.psi.beta, nsp),
                         sigma=(sigma.psi.beta)*lambda.mat, #Will need to revisit this
                                        #to make sure it's consistent
                                        #throughout entire script
                                        #(i.e. not sigma.psi.beta^2)
                                        #as coded this should not match
                                        #lme4 output
                         log=TRUE)
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
  used <- c(## psi.0.terms,
            p.0.terms,
            psi.beta.terms,
            psi.site.terms,
            psi.year.terms) 

  TT <- dmvnorm(random.draws[chain, used],
                mean=mu.vec[used],
                sigma=vcov.mat[used, used],
                log=TRUE)
  LL.theta.given.MU.and.SIGMA <- TT
  
  ## calculate the total log likelihood.
  ## LL.psi.0.given.theta +
  LL.psi.beta.given.theta +
    LL.psi.site.given.theta +
      LL.psi.year.given.theta +
        LL.p.given.theta + 
          LL.X.given.psi.and.p -
            LL.theta.given.MU.and.SIGMA
  
}




### LoF. 22-Dec-2014: Adding version of calc.LL function for a simple
## model that does not estimate detectability

par.calc.LL.no.det <- function(chain, random.draws, prms, mu.vec, vcov.mat) {

  require(mvtnorm)
  expit <- function(x) 1/(1+exp(-x))

  ## get important data into easily accessible form
  env <- dd.simulated$env
  nsite <- dd.simulated$nsite
  nyr <- dd.simulated$nyr
  nrep <- max(dd.simulated$nrep) ## note: assuming same num of reps
  nsp <- dd.simulated$nsp
  X <- dd.simulated$X
  lambda.mat <- prms$lambda*dd.simulated$VCOV +
    (1-prms$lambda)*dd.simulated$ID

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
                         sigma=(sigma.psi.beta)*lambda.mat, #Will need to revisit this
                                        #to make sure it's consistent
                                        #throughout entire script
                                        #(i.e. not sigma.psi.beta^2)
                                        #as coded this should not match
                                        #lme4 output
                         log=TRUE)

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

###################################
## Currently the importance sampling is giving me -Inf.
## Maybe an error related to incomplete sampling???
## Keep an eye on this....
###################################
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


