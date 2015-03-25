##########################################
###  PHYLOGENETIC OCCUPANCY MODEL SIMULATION
###
### August 14 2014
###
### August 25 2014: Updated notation, and changed image colors to grey
###  scale
###
### August 26 2014: Included optional scaling argument, to determine
###  whether the mvnorm distribution is first scaled before plugging
###  into the psi function
###
### September 4 2014: Updated XX so that the dimensions are the same
###  as that required for the JAGS models.
###
### September 9 2014: Updated to remove scaling business, and to
### include optional lambda parameter (activated by setting type="lambda"
### and then providing a lambda value) rather than
### phylogenetic effect versus random effect of species competing with
### one another. More work needed to update beta summary to return the
### relevant portions, but that can be done another day
###
### September 12 2014: updated notation to match with JAGS model

### September 19 2014: Swapped use.lambda for type=c("lambda",
### "OU","relative") included specification for OU parameters
##########################################

source('~/Dropbox/ccb_banding/analyses/Phylogenetic Mixed Effect Model/myfunctions.R')
library(ape)
library(MASS)

##########################################
###  THE  SIMULATION  FUNCTION
##########################################

pom.sim <- function(nsite=15, ## number of sites
                    nsp=50, ## number of species
                    nyr=5, ## number of years
                    nrep=5, ## number of replicates per year
                    mu.psi.0=0, ## mean occupancy
                    sigma.psi.0=0.5, ## variation in occupancy
                    mu.beta=0, ## mean response to environment
                    sigma.beta.spp=1, ## non-phylogenetic variation in
                                      ## response to environment.  Not
                                      ## used with lambda! Rename
                                      ## later on
                    sigma.psi.site=0.1, ## width of site random effect
                                        ## in occupancy 
                    sigma.psi.year=0.1, ## width of year random effect
                                        ## in occupancy 
                    mu.p.0=0, ## mean detectability
                    sigma.p.0=0.5, ## variation in detectability
                    sigma.spp.trait=1, ## variation in species' trait
                                       ## values 
                    psi.0.trait=0, ## effect of trait on mean
                                   ## occupancy (i.e., slope of
                                   ## occupancy with regard to
                                   ## species' traits)
                    beta.trait=0, ## effect of trait on species'
                                  ## responses to the environment
                    tree=NULL, ## input a tree
                    sigma.beta.sp=1, ## variation in species'
                                     ## responses to the environment.
                                     ## This scales the width of the
                                     ## covariance matrix that is used
                                     ## to draw species' responses to
                                     ## the environment
                    lambda=1, ## weighting on phylogenetic covariance
                              ##  matrix (e.g., with lambda=1, the
                              ##  covariance matrix is given full
                              ##  weight and with lambda=0, the model
                              ##  reduces to a simple random effect of
                              ##  species)
                    type=c("lambda", "OU", "relative"),
                    alpha.ou=1, ## Parameter dictating decay of
                                ## phylogenetic signal under the OU
                                ## model. Used only when type="OU"
                    site.structure='linear', ## environmental gradient
                                             ## type 
                    psi.image=FALSE, ## show plots while creating data
                    occ.image=FALSE,
                    det.image=FALSE,
                    plot.tree=FALSE) {

  if(site.structure=='linear')
    site.env <- seq(-1,1, length.out=nsite)
  
  env.pre <- rep(site.env, each=nsp)
  env <- rep(env.pre, nyr)

  ## Generate trait values for all species
  spp.trait <- rnorm(nsp, sd=sigma.spp.trait)
  
  ## Generate phylogenetic relatedness for all speceis
  if(is.null(tree)==TRUE) {
    tree <- rcoal(nsp, tip.label=1:nsp)
  }

  tree.vcv <- vcv(tree, corr=T)
  i <- order(as.numeric(rownames(tree.vcv)))
  j <- order(as.numeric(colnames(tree.vcv)))
  tree.vcv <- tree.vcv[i,j]
  
  ## Species specific random intercept in occupancy
  ##
  ## Generate species intercepts based on draws from a normal
  ## distribution with mean mu.psi.0 and standard deviation
  ## sigma.psi.0
  ## alpha.spp <- rnorm(nsp, mean=mu.psi.0, sd=sigma.psi.0) +
  ##   psi.0.trait*spp.trait 
  ## alpha.pre <- rep(alpha.spp, nsite)
  ## alpha <- rep(alpha.pre, nyr)
  alpha <- rep(mu.psi.0, nyr*nsp*nsite)

  ## Generate species response slopes to environmental gradient, based
  ## on random draws from a normal distribution with mean mu.beta
  ## and standard deviation sigma.beta.spp
  if(type=="lambda") {
    VCOV.mat <- sigma.beta.sp*(tree.vcv*lambda + (1-lambda)*diag(nsp))
    ran.beta.spp <- 0
  } 

  if(type=="OU") {
    OU <- function(sig, alpha, vcv)
      ((sig^2)/(2*alpha)) * exp(-2*alpha*(max(vcv)-vcv)) *
        (1-exp(-2*alpha*vcv)) 
    VCOV.mat <- OU(sig=sigma.beta.sp, alpha=alpha.ou, vcv=tree.vcv)
    ran.beta.spp <- 0
  } 
    
  if(type=="relative") {
    VCOV.mat <- tree.vcv
    ran.beta.spp <- rnorm(nsp, mean=0, sd=sigma.beta.spp)
  }
  
 
  
  mu.beta.sp <- rep(mu.beta, nsp) + beta.trait*spp.trait
  psi.beta.sp <- mvrnorm(mu=mu.beta.sp,
                         Sigma=VCOV.mat) 

  ## ran.beta.spp below is 0 in the case of type="lambda" or "OU"
  beta.spp <- ran.beta.spp + psi.beta.sp 
  
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
          main="Mean Psi across yrs",
          xlab="Sites", ylab="Spp", col=rev(gray.colors(12))) 
  }

  ## Simulate true occupancy states for each spp at each site for each
  ## year. Drawn from binomial dist with size 1 and mean psi.
  ZZ <- rbinom(nsp*nsite*nyr, 1, psi)
  ## Convert to array for downstream
  ZZ.A <- array(ZZ, dim=c(nsp,nsite,nyr)) 

  if(occ.image==TRUE) {
    quartz()
    image(t(apply(ZZ.A, c(1,2) , mean)),
          main="Mean Occupancy States accros yrs",
          xlab="Sites", ylab="Spp",
          col=rev(gray.colors(12))) 
  }

  ## Simulate detection proccess, based on draws from a normal
  ## distribution with mean mu.p.0 and standard deviation sigma.p.0. Run
  ## through an inverse logit link to make sure it goes from 0 to 1.
  det.spp <- invlogit( rnorm(nsp, mean=mu.p.0, sd=sigma.p.0) )
  det.pre <- rep(det.spp, nsite)
  det <- rep(det.pre, nyr)

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
          main="Mean Detection States across replicates and yrs",
          xlab="Sites", ylab="Spp", col=rev(gray.colors(12))) 
  }
  
  if(plot.tree==T) {
    quartz()
    
    colsTips <- colorize(beta.spp,
                         colors=c("red","gold","darkgreen"),
                         min= min(beta.spp), max=max(beta.spp)) 

    phycomp <- psi.beta.sp
    suppressWarnings(anres <- ace(phycomp, phy=tree, method="REML", CI=FALSE))
    
    
    allnodes <- c(phycomp, anres$ace)

    colsEdges <- colorize(allnodes,
                          colors=c("red","gold","darkgreen"),
                          min=min(beta.spp), max=max(beta.spp))

    plot.phylo(tree, type="p", show.tip.label=F,
               edge.color= colsEdges[match(tree$edge[,1], names(allnodes))]) 
    tiplabels(pch=16, col=colsTips[tree$tip.label])

    print("Tree tips colored according to overall species response to environmental gradient (i.e., beta.spp), branches colored according to the phylogenetic component of (i.e., sigma.beta.sp * tree.vcv)")
  }
  
  if(type=="lambda") {
    return(list(ZZ=ZZ.A,
                XX=XX,
                nsite=nsite,
                nyr=nyr,
                nsp=nsp,
                psi=psi.A,
                env=site.env,
                trait=spp.trait,
                tree=tree,
                tree.vcv=tree.vcv,
                det.spp= det.spp,
                beta.summary=data.frame(beta.spp,
                  mu.beta,
                  mu.beta.sp,
                  beta.trait,
                  spp.trait
                  )))
  } 
  
    if(type=="OU") {
    return(list(ZZ=ZZ.A,
                XX=XX,
                nsite=nsite,
                nyr=nyr,
                nsp=nsp,
                psi=psi.A,
                env=site.env,
                trait=spp.trait,
                tree=tree,
                tree.vcv=tree.vcv,
                det.spp= det.spp,
                beta.summary=data.frame(beta.spp,
                  mu.beta,
                  mu.beta.sp,
                  beta.trait,
                  spp.trait
                  )))
  } 
  
  if (type=="relative") {
    return(list(ZZ=ZZ.A,
                XX=XX,
                nsite=nsite,
                nyr=nyr,
                nsp=nsp,
                psi=psi.A,
                env=site.env,
                trait=spp.trait,
                tree=tree,
                tree.vcv=tree.vcv,
                det.spp= det.spp,
                beta.summary=data.frame(beta.spp,
                  mu.beta,
                  beta.trait,spp.trait,
                  ran.beta.spp,
                  psi.beta.sp)))
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

  env <- matrix(test$env, ncol=1)
  colnames(env) <- 'fc.100'

  nrep <- apply(test$XX, c(1,2,4),
                function(x) sum(x>=0,na.rm=TRUE))

  z.init <- apply(test$X, c(1,2,4),
                  function(x) (sum(x,na.rm=TRUE)>0)*1)
  z.init[apply(test$X, c(1,2,4), function(x) !any(!is.na(x)))] <- NA

  ## package up data so that we can run model on it
  list(X=test$XX,
       VCOV=test$tree.vcv,
       ID=diag(1, nrow=test$nsp, ncol=test$nsp),
       VV.mu=matrix(0, nrow=2, ncol=test$nsp),
       z.init=z.init,
       nsite=test$nsite, 
       nyr=test$nyr,
       nrep=nrep,
       nsp=test$nsp,
       date.mat=array(0, dim=c(test$nsite, test$nyr, max(nrep))),  ## not used
       site.presence=site.presence,
       frac.presence=frac.presence,
       landcover=env,
       sim.trait=test$trait)
}
