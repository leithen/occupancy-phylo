## ************************************************************
## simulate data
## ************************************************************

source('~/Dropbox/ccb_banding/analyses/Phylogenetic Mixed Effect Model/myfunctions.R')
library(ape)
library(MASS)

##########################################
###  THE  SIMULATION  FUNCTION
##########################################

pom.sim <- function(nsite=15, ## number of sites
                    nsp=64, ## number of species
                    nyr=5, ## number of years
                    nrep=5, ## number of replicates per year
                    mu.psi.0=0, ## mean occupancy
                    sigma.psi.0=0.5, ## variation in occupancy
                    mu.psi.beta=0, ## mean response to environment
                    sigma.beta.non.phylo=1, ## non-phylogenetic
                                            ## variation in response
                                            ## to environment.  Not
                                            ## used with lambda!
                    sigma.psi.site=0.1, ## width of site random effect
                                        ## in occupancy 
                    sigma.psi.year=0.1, ## width of year random effect
                                        ## in occupancy 
                    mu.p.0=0, ## mean detectability
                    sigma.p.0=0.5, ## variation in detectability
                    psi.0.trait=0, ## effect of trait on mean
                                   ## occupancy (i.e., slope of
                                   ## occupancy with regard to
                                   ## species' traits)
                    beta.trait=0, ## effect of trait on species'
                                  ## responses to the environment
                    tree=NULL, ## input a tree
                    sigma.psi.beta=1, ## variation in species'
                                     ## responses to the environment.
                                     ## This scales the width of the
                                     ## covariance matrix that is used
                                     ## to draw species' responses to
                                     ## the environment. sigma.psi.beta
                                     ## is squred within the mvnorm
                                     ## function to convert it to a
                                     ## variance
                    lambda=1, ## weighting on phylogenetic covariance
                              ## matrix (e.g., with lambda=1, the
                              ## covariance matrix is given full
                              ## weight and with lambda=0, the model
                              ## reduces to a simple random effect of
                              ## species)
                    lambda.trait=1, ## the degree to wich the trait is
                                    ## phylogenetically conserved
                    psi.0s=NULL,                
                    use.random.int=TRUE,
                    use.lambda=TRUE,
                    imperfect.detection=TRUE, ## turn on or off
                                              ## imperfect detection
                    site.structure='linear', ## environmental gradient
                                             ## type 
                    psi.image=FALSE, ## show plots while creating data
                    occ.image=FALSE,
                    det.image=FALSE,
                    plot.tree=FALSE) {

  inputs <- list(nsite=nsite,
                 nsp=nsp,
                 nyr=nyr,
                 nrep=nrep,
                 mu.psi.0=mu.psi.0,
                 sigma.psi.0=sigma.psi.0,
                 mu.psi.beta=mu.psi.beta,
                 sigma.beta.non.phylo=sigma.beta.non.phylo,
                 sigma.psi.site=sigma.psi.site,
                 sigma.psi.year=sigma.psi.year,
                 mu.p.0=mu.p.0,
                 sigma.p.0=sigma.p.0,
                 psi.0.trait=psi.0.trait,
                 beta.trait=beta.trait,
                 tree=tree,
                 sigma.psi.beta=sigma.psi.beta,
                 lambda=lambda,
                 lambda.trait=lambda.trait,
                 psi.0s=psi.0s,                
                 use.random.int=use.random.int,
                 use.lambda=use.lambda,
                 imperfect.detection=imperfect.detection,
                 site.structure=site.structure,
                 psi.image=psi.image,
                 occ.image=occ.image,
                 det.image=det.image,
                 plot.tree=plot.tree)
  
  if(site.structure=='linear')
    site.env <- seq(-1,1, length.out=nsite)
  site.env <- as.numeric(scale(site.env))
  
  env.pre <- rep(site.env, each=nsp)
  env <- rep(env.pre, nyr)
  
  ## Generate phylogenetic relatedness for all speceis
  if(is.null(tree)==TRUE) {
    tree <- rcoal(nsp, tip.label=1:nsp)
  }

  tree.vcv <- vcv(tree, corr=T)
  i <- order(as.numeric(rownames(tree.vcv)))
  j <- order(as.numeric(colnames(tree.vcv)))
  tree.vcv <- tree.vcv[i,j]

  ## Generate trait values for all species
  VCOV.trait.mat <- tree.vcv*lambda.trait + (1-lambda.trait)*diag(nsp)
  spp.trait <- mvrnorm(mu=rep(0, nsp),
                       Sigma=VCOV.trait.mat)
  spp.trait <- (spp.trait-mean(spp.trait))/sd(spp.trait)
  
  ## Species specific random intercept in occupancy
  ##
  ## Generate species intercepts based on draws from a normal
  ## distribution with mean mu.psi.0 and standard deviation
  ## sigma.psi.0
  if (use.random.int==T) {
  alpha.spp <- rnorm(nsp, mean=mu.psi.0, sd=sigma.psi.0) +
    psi.0.trait*spp.trait 
	}
	
	else {
		if(length(psi.0s)!=nsp) {warning("Length of psi.0s does not equal number of species!")}
		
		alpha.spp<-psi.0s
		}

  alpha.pre <- rep(alpha.spp, nsite)
  alpha <- rep(alpha.pre, nyr)

  ## Generate species response slopes to environmental gradient, based
  ## on random draws from a normal distribution with mean mu.psi.beta
  ## and standard deviation sigma.beta.non.phylo
  if(use.lambda) {
    VCOV.mat <- tree.vcv*lambda + (1-lambda)*diag(nsp)
    ran.beta.spp <- 0
  } else {
    VCOV.mat <- tree.vcv
    ran.beta.spp <- rnorm(n=nsp, mean=0, sd=sigma.beta.non.phylo)
  }

  ####### EDITS March 24th 2015 ################################
  ## There are multiple ways to think about phylogenetic effects.
  ## Originally we coded them as follows:
  ##
  ## mu.beta.sp <- rep(mu.psi.beta, nsp) + beta.trait*spp.trait
  ## psi.beta.sp <- mvrnorm(mu=mu.beta.sp,
  ##                       Sigma=(sigma.psi.beta^2)*VCOV.mat) 
  ##
  ## ## ran.beta.spp below is 0 in the case of use.lambda
  ## beta.spp <- ran.beta.spp + psi.beta.sp 

  ## To conform with how they're usually thought of (phylogenetic
  ## correlation in the residuals after accounting for other 
  ## variation) I have reparameterized things as follows:

  zeros.sp <- rep(0, nsp)
  psi.beta.sp.rand<- mvrnorm(mu=zeros.sp, 
                             Sigma=(sigma.psi.beta^2)*VCOV.mat)
  
  beta.spp <- mu.psi.beta + beta.trait*spp.trait + psi.beta.sp.rand
  
  mu.beta.sp <- mu.psi.beta + 
  					beta.trait*spp.trait # This is the expectation
  										 # prior to random effects
  										 # of spp and phylogeney
  													
  
  ## Moving mu.psi.beta out of the random effect should have not 
  ## real effect on things, but is heuristically conventient if 
  ## we're thinking of phylogeney as a purely random effect to 
  ## explain correlation structure. Models in 
  ## ccb_banding/analyses/occupancy-phylo/src/models/traits have 
  ## been updated to match this parameterization, though non-trait 
  ## models in ccb_banding/analyses/occupancy-phylo/src/models have 
  ## not been. (Though it should not matter, as these parameterization 
  ## should be mathamatically equivalent when there is no trait)
  ####### End Edits March 24 2015 ################################ 
  
  
  beta.pre <- rep(beta.spp, nsite)
  beta <- rep(beta.pre, nyr)

  ## Generate random effect of site overall occupancy, drawing from
  ## random distribution with mean zero and standard deviation
  ## sigma.psi.site
  site.effect <- rnorm(nsite, mean=0, sd=sigma.psi.site)
  ran.site.pre <- rep(site.effect, each=nsp)
  ran.site <- rep(ran.site.pre, nyr)

  ## Generate random effect of sampling year, drawn from a normal
  ## distribution with mean zero and standard deviation sigma.psi.year
  yr.effect <- rnorm(nyr, mean=0, sd=sigma.psi.year)
  ran.yr <- rep(yr.effect, each=nsp*nsite)

  ## Calculate occupancy probability for each site in each year for
  ## every species
  psi <- invlogit(alpha + beta*env + ran.site + ran.yr)
  
  ## Convert to a spp by site by yr array for downstream
  psi.A <- array(psi, dim=c(nsp,nsite,nyr))

  if(psi.image==TRUE) {
    quartz()
    image(t(apply(psi.A, c(1,2) , mean)),
          main='Mean Psi across yrs',
          xlab='Sites', ylab='Spp', col=rev(gray.colors(12))) 
  }

  ## Simulate true occupancy states for each spp at each site for each
  ## year. Drawn from binomial dist with size 1 and mean psi.
  ZZ <- rbinom(nsp*nsite*nyr, 1, psi)
  ## Convert to array for downstream
  ZZ.A <- array(ZZ, dim=c(nsp,nsite,nyr)) 

  if(occ.image==TRUE) {
    quartz()
    image(t(apply(ZZ.A, c(1,2) , mean)),
          main='Mean Occupancy States accros yrs',
          xlab='Sites', ylab='Spp',
          col=rev(gray.colors(12))) 
  }

  ## Simulate detection proccess, based on draws from a normal
  ## distribution with mean mu.p.0 and standard deviation sigma.p.0. Run
  ## through an inverse logit link to make sure it goes from 0 to 1.
  det.spp <- invlogit(rnorm(n=nsp, mean=mu.p.0, sd=sigma.p.0))
  det.pre <- rep(det.spp, nsite)
  ## set detectability to 1 for the case where we do not want
  ## imperfect detetion
  if(imperfect.detection==FALSE)
    det.pre <- rep(1, length(det.pre))
  
  ## Generate Observation Matrix based on the number of replicates in
  ## each site in each year. Based on true occupancy states
  ## discounting detection probability
  XX <- array(NA, dim=c(nsp, nsite, nyr, nrep))
  for(yr in 1:nyr) {
    for(rep in 1:nrep) {
      ## Simulate observation proceess
      hold <- rbinom(nsp*nsite, 1, prob=as.vector(ZZ.A[,,yr])*det.pre)
      XX[,,yr,rep] <- hold
    }
  }

  dimnames(XX) <- list(species=paste('sp.', 1:nsp, sep=''),
                       site=paste('site.', 1:nsite, sep=''),
                       year=2000+1:nyr,
                       rep=1:nrep)

  ## permute above so that it conforms with JAGS models
  XX <- aperm(XX, c(2,3,4,1))
  
  if(det.image==TRUE) {
    quartz()
    image(t(apply(XX, c(4,1) , mean)),
          main='Mean Detection States across replicates and yrs',
          xlab='Sites', ylab='Spp', col=rev(gray.colors(12))) 
  }
  
  if(plot.tree==T) {
    quartz()
    
    colsTips <- colorize(beta.spp,
                         colors=c('red','gold','darkgreen'),
                         min= min(beta.spp), max=max(beta.spp)) 

    phycomp <- psi.beta.sp.rand ###NOTE: March 24th 2015. This may
    							### In some cases need to be psi.beta.sp
    							
    suppressWarnings(anres <- ace(phycomp, phy=tree, method='REML', CI=FALSE))
    
    
    allnodes <- c(phycomp, anres$ace)

    colsEdges <- colorize(allnodes,
                          colors=c('red','gold','darkgreen'),
                          min=min(beta.spp), max=max(beta.spp))

    plot.phylo(tree, type='p', show.tip.label=F,
               edge.color= colsEdges[match(tree$edge[,1], names(allnodes))]) 
    tiplabels(pch=16, col=colsTips[tree$tip.label])

    print('Tree tips colored according to overall species response to environmental gradient (i.e., beta.spp), branches colored according to the phylogenetic component of (i.e., sigma.psi.beta * tree.vcv)')
  }

## ********************************************************************
  ## Creating Dataframe to plug into lmer

	df<-data.frame(
  			presence=as.vector(XX),
  			expand.grid(site=1:nsite, year=1:nyr, rep=1:nrep, spp=1:nsp),
  			env=rep(site.env, nyr*nrep*nsp)
  			)
  			
## ********************************************************************  
  
  
  if(use.lambda) {  ###NOTE March 24th 2015. These returned inputs may need to be more comphehensive
    return(c(inputs,
             list(alpha.spp=alpha.spp,
                  beta.spp=beta.spp,
                  det.spp= det.spp,
                  env=site.env,
                  mu.beta.sp=mu.beta.sp,
                  psi=psi.A,
                  psi.beta.sp.rand=psi.beta.sp.rand,
                  site.effect=site.effect,
                  trait=spp.trait,
                  tree=tree,
                  tree.vcv=tree.vcv,
                  yr.effect=yr.effect,
                  XX=XX,
                  ZZ=ZZ.A,
                  df=df)))
  } else {
    return(c(inputs,
             list(alpha.spp=alpha.spp,
                  beta.spp=beta.spp,
                  det.spp= det.spp,
                  env=site.env,
                  mu.beta.sp=mu.beta.sp,
                  psi=psi.A,
                  psi.beta.sp=psi.beta.sp,
                  ran.beta.spp=ran.beta.spp,
                  site.effect=site.effect,
                  trait=spp.trait,
                  tree=tree,
                  tree.vcv=tree.vcv,
                  yr.effect=yr.effect,
                  XX=XX,
                  ZZ=ZZ.A,
                  df=df)))
  }
}

## package up the data and create additional structures so that it is
## ready for use in the JAGS model
make.data <- function(...) {

  test <- pom.sim(...)
  
  ## create site x species matrix indicating site presence
  site.presence <- (apply(test$X, c(1,4), sum)>0)*1

  ## create site x species matrix indicating fraction of times a
  ## species was present
  frac.presence.num <- apply(test$X, c(1,4), sum, na.rm=TRUE)
  frac.presence.denom <- apply(test$X, c(1,4), function(x) sum(!is.na(x)))
  frac.presence <- frac.presence.num/frac.presence.denom
  pres.yr.1 <- (apply(test$X[,1,,], c(1,3), sum, na.rm=TRUE)>0)*1
  frac.presence[pres.yr.1==1] <- 1

  nrep <- apply(test$XX, c(1,2,4),
                function(x) sum(x>=0,na.rm=TRUE))

  Z <- apply(test$X, c(1,2,4),
             function(x) (sum(x,na.rm=TRUE)>0)*1)
  Z[apply(test$X, c(1,2,4), function(x) !any(!is.na(x)))] <- NA

  ## package up data so that we can run model on it
  list(inputs=test,
       date.mat=array(0, dim=c(test$nsite, test$nyr, max(nrep))),  ## not used
       env=test$env,
       frac.presence=frac.presence,
       ID=diag(1, nrow=test$nsp, ncol=test$nsp),
       nrep=nrep,
       nsite=test$nsite, 
       nsp=test$nsp,
       nyr=test$nyr,
       site.presence=site.presence,
       trait=test$trait,
       VCOV=test$tree.vcv,
       VV.mu=matrix(0, nrow=2, ncol=test$nsp),
       X=test$XX,
       Z=Z)
}
