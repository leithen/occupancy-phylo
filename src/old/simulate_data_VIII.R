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
### include optional lambda parameter (activated by setting use.lambda
### to TRUE, and then providing a lambda value) rather than
### phylogenetic effect versus random effect of species competing with
### one another. More work needed to update beta summary to return the
### relevant portions, but that can be done another day

###September 12 2014: Updated notation to match with JAGS model
##########################################

source('~/Dropbox/ccb_banding/analyses/Phylogenetic Mixed Effect Model/myfunctions.R')
library(ape)
library(MASS)

##########################################
###  THE  SIMULATION  FUNCTION
##########################################

pom.sim <- function(nsite=10,
                    nsp=20,
                    nyr=2,
                    nrep=6,
                    mu.psi.0=0,
                    sigma.psi.0=1,
                    mu.beta=0,
                    sd.beta.spp=10, #Not used with lambda! Rename later on
                    sigma.psi.site=0,
                    sigma.psi.year=0,
                    mu.p.0=10,
                    sigma.p.0=0,
                    sd.spp.trait=1,
                    psi.0.trait=0,
                    beta.trait=0,
                    tree=NULL,
                    sigma.beta.sp=1,
                    psi.image=TRUE,
                    occ.image=TRUE,
                    det.image=TRUE,
                    plot.tree=TRUE,
                    use.lambda=FALSE,
                    lambda=1) {

  ## site.env <- rnorm(nsite) #Randomly select site environmental gradient values
  site.env <- seq(-3,3, length.out=nsite) #Or keep it structured for ease
  env.pre <- rep(site.env, each=nsp)
  env <- rep(env.pre, nyr)

  ## Generate trait values for all species
  spp.trait <- rnorm(nsp, sd=sd.spp.trait)
  
  ## Generate phylogenetic relatedness for all speceis
  if(is.null(tree)==TRUE) {
    tree <- rcoal(nsp, tip.label=1:nsp)
  }

  tree.vcv <- vcv(tree, corr=T)
  i <- order(as.numeric(rownames(tree.vcv)))
  j <- order(as.numeric(colnames(tree.vcv)))
  tree.vcv <- tree.vcv[i,j]
  
  #### Species specific random intercept in occupancy
  ## Generate species intercepts based on draws from a normal
  ## distribution with mean mu.psi.0 and standard deviation
  ## sigma.psi.0
  alpha.spp <- rnorm(nsp, mean=mu.psi.0, sd=sigma.psi.0) +
    psi.0.trait*spp.trait 
  alpha.pre <- rep(alpha.spp, nsite)
  alpha <- rep(alpha.pre, nyr)

  ## Generate species response slopes to environmental gradient, based
  ## on random draws from a normal distribution with mean mu.beta
  ## and standard deviation sd.beta.spp
  if(use.lambda) {
    VCOV.mat <- tree.vcv*lambda + (1-lambda)*diag(nsp)
    ran.beta.spp <- 0
  } else {
    VCOV.mat <- tree.vcv
    ran.beta.spp <- rnorm(nsp, mean=0, sd=sd.beta.spp)
  }
  
  mu.beta.sp <- rep(mu.beta, nsp) + beta.trait*spp.trait
  psi.beta.sp <- mvrnorm(mu=mu.beta.sp,
                            Sigma=sigma.beta.sp*VCOV.mat) 

###ran.beta.spp is 0 in the case of use.lambda
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

    print("Tree tips colored according to overal species response to environmental gradient (i.e. beta.spp), branches colored according to the phylogenetic component of (i.e. sigma.beta.sp * tree.vcv)")
  }
  
  if(use.lambda) {
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
  } else {
    return(list(ZZ=ZZ.A,
                XX=XX,
                nsite=nsite,
                nyr=nyr,
                nsp=nsp,
                psi=psi.A,
                env=site.env,
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
