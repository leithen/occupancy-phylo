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
##########################################

rm(list=ls())
setwd('~/Dropbox/ccb_banding/analyses')
source('Phylogenetic Mixed Effect Model/myfunctions.R', chdir = TRUE)
library(ape)
library(MASS)

##########################################
###  THE  SIMULATION  FUNCTION
##########################################

pom.sim <- function(nsite=10,
                    nsp=20,
                    nyr=2,
                    nrep=6,
                    mu.alpha.spp=0,
                    sd.alpha.spp=1,
                    mu.beta.spp=0,
                    sd.beta.spp=10,
                    sd.ran.site=0,
                    sd.ran.yr=0,
                    mu.det=10,
                    sd.det=0,
                    sd.spp.trait=1,
                    alpha.trait.effect=0,
                    beta.trait.effect=0,
                    tree=NULL,
                    beta.phylo.effect=1,
                    psi.image=TRUE,
                    occ.image=TRUE,
                    det.image=TRUE,
                    plot.tree=TRUE,
                    use.lambda=FALSE,
                    lambda=1) {

  ## site.env <- rnorm(nsite) #Randomly select site environmental gradient values
  site.env <- seq(-3,3, length.out=nsite) #Or keep it structured for ease
  landcover.pre <- rep(site.env, each=nsp)
  landcover <- rep(landcover.pre, nyr)

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
  
  ## Generate species intercepts based on draws from a normal
  ## distribution with mean mu.alpha.spp and standard deviation
  ## sd.alpha.spp
  alpha.spp <- rnorm(nsp, mean=mu.alpha.spp, sd=sd.alpha.spp) +
    alpha.trait.effect*spp.trait 
  alpha.pre <- rep(alpha.spp, nsite)
  alpha <- rep(alpha.pre, nyr)

  ## Generate species response slopes to environmental gradient, based
  ## on random draws from a normal distribution with mean mu.beta.spp
  ## and standard deviation sd.beta.spp
  if(use.lambda) {
    VCOV.mat <- tree.vcv*lambda + (1-lambda)*diag(nsp)
    ran.beta.spp <- 0
  } else {
    VCOV.mat <- tree.vcv
    ran.beta.spp <- rnorm(nsp, mean=0, sd=sd.beta.spp)
  }
  
  mu.phylo <- rep(mu.beta.spp, nsp) + beta.trait.effect*spp.trait
  ran.beta.phylo <- mvrnorm(mu=mu.phylo,
                            Sigma=beta.phylo.effect*VCOV.mat) 
  beta.spp <- ran.beta.spp + ran.beta.phylo 
  
  beta.pre <- rep(beta.spp, nsite)
  beta <- rep(beta.pre, nyr)

  ## Generate random effect of site overall occupancy, drawing from
  ## random distribution with mean zero and standard deviation
  ## sd.ran.site
  site.effect <- rnorm(nsite, mean=0, sd=sd.ran.site)
  ran.site.pre <- rep(site.effect, each=nsp)
  ran.site <- rep(ran.site.pre, nyr)

  ## Generate random effect of sampling year, drawn from a normal
  ## distribution with mean zero and standard deviation sd.ran.yr
  yr.effect <- rnorm(nyr, mean=0, sd=sd.ran.yr)
  ran.yr <- rep(yr.effect, each=nsp*nsite)

  ## Calculate occupancy probability for each site in each year for
  ## every species
  psi <- invlogit(alpha + beta*landcover + ran.site + ran.yr)

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
  ## distribution with mean mu.det and standard deviation sd.det. Run
  ## through an inverse logit link to make sure it goes from 0 to 1.
  det.spp <- invlogit( rnorm(nsp, mean=mu.det, sd=sd.det) )
  det.pre <- rep(det.spp, nsite)
  det <- rep(det.pre, nyr)

  ## Generate Observation Matrix based on the number of replicates in
  ## each site in each year. Based on true occupancy states
  ## discounting detection probability
  XX <- array(NA, dim=c(nsp, nsite, nyr, nrep))
  for(yr in 1:nyr) {
    for(rep in 1:nrep) {
      ## Simulate observation proceess
      hold <- rbinom(nsp*nsite, 1, as.vector(ZZ.A[,,yr])*det) 
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

    phycomp <- ran.beta.phylo
    suppressWarnings(anres <- ace(phycomp, phy=tree, method="REML", CI=FALSE))
    
    
    allnodes <- c(phycomp, anres$ace)

    colsEdges <- colorize(allnodes,
                          colors=c("red","gold","darkgreen"),
                          min=min(beta.spp), max=max(beta.spp))

    plot.phylo(tree, type="p", show.tip.label=F,
               edge.color= colsEdges[match(tree$edge[,1], names(allnodes))]) 
    tiplabels(pch=16, col=colsTips[tree$tip.label])

    print("Tree tips colored according to overal species response to environmental gradient (i.e. beta.spp), branches colored according to the phylogenetic component of (i.e. beta.phylo.effect * tree.vcv)")
  }
  
  if(use.lambda) {
    return(list(ZZ=ZZ.A,
                XX=XX,
                nsite=nsite,
                nyr=nyr,
                nsp=nsp,
                psi=psi.A,
                env=site.env,
                landcover=site.env,
                trait=spp.trait,
                tree=tree,
                tree.vcv=tree.vcv,
                det.spp= det.spp,
                beta.summary=data.frame(beta.spp,
                  mu.beta.spp,
                  mu.phylo,
                  beta.trait.effect,
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
                landcover=site.env,
                trait=spp.trait,
                tree=tree,
                tree.vcv=tree.vcv,
                det.spp= det.spp,
                beta.summary=data.frame(beta.spp,
                  mu.beta.spp,
                  beta.trait.effect,spp.trait,
                  ran.beta.spp,
                  ran.beta.phylo)))
  }
}


##########################################
###  SIMULATE  THE  DATA 
##########################################
set.seed(1)

## Set tree so as to rerun multiple times but control the tree
nsp <- 50
tree1 <- rcoal(nsp, tip.label=1:nsp) 

test <- pom.sim(
### Sampling set up
  nsite=18, 
  nsp=nsp,
  nyr=6,
  nrep=6,
  sd.spp.trait=1, ## Set variation in trait values for species. 1 corresponds to a random standardized continuous trait
  tree=tree1, ## Input tree to use, or set as NULL to use a random tree  
  
### Effects relating to species intercepts (i.e. mean occupancy values)         
  mu.alpha.spp=0,
  sd.alpha.spp=5,
  alpha.trait.effect=0, 

### Effects relating to species slopes (i.e. response to the environmental gradient)         
  mu.beta.spp=0,
  sd.beta.spp=1, #Not used if use.lambda=T, only pertains to species specific random effect
  beta.trait.effect=0, #Set trait influence to zero for trouble shooting stage.
  beta.phylo.effect=50, #This is maybe more correctly thought of as sigma.beta.phylo now, given that it's a multiplier for the correlation matrix.
  
### Random effects of site and year        
  sd.ran.site=0.2, 
  sd.ran.yr=0.1,
  
###Effects related to detectability        
  mu.det=2, ## Detectability, on logit scale!
  sd.det=1,
  
###Use lambda specification if TRUE. Otherwise a random effect of species competing with a phylogenetic effect of species (as in Ives and Helmus)	
  use.lambda=TRUE,
  lambda=1,
  
###Set image printing         
  occ.image=F,                 
  psi.image=F,
  det.image=T,
  plot.tree=T) 


##########################################
###  LOOK  AT  THE  DATA 
##########################################

test$psi #Psi array
test$ZZ #True occupancy array
test$XX #Detection array
test$env #Environmental gradient corresponding to sites
test$trait #Trait values corresponding to species
test$tree #The phylogenetic tree used
test$tree.vcv #The correlation matrix from the phylogenetic tree used
test$det.spp #Species specific detectibilities
test$beta.summary #Summary of model components per species that go into generating beta.spp, or overall species response to the environmental gradient


## create site x species matrix indicating site presence
site.presence <- (apply(test$X, c(1,4), sum)>0)*1

## create site x species matrix indicating fraction of times a
## species was present
frac.presence.num <- apply(test$X, c(1,4), sum, na.rm=TRUE)
frac.presence.denom <- apply(test$X, c(1,4), function(x) sum(!is.na(x)))
frac.presence <- frac.presence.num/frac.presence.denom
pres.yr.1 <- (apply(test$X[,1,,], c(1,3), sum, na.rm=TRUE)>0)*1
frac.presence[pres.yr.1==1] <- 1

landcover <- matrix(test$landcover, ncol=1)
colnames(landcover) <- 'fc.100'

nrep <- apply(test$XX, c(1,2,4),
              function(x) sum(x>=0,na.rm=TRUE))

z.init <- apply(test$X, c(1,2,4),
                function(x) (sum(x,na.rm=TRUE)>0)*1)
z.init[apply(test$X, c(1,2,4), function(x) !any(!is.na(x)))] <- NA

## package up data so that we can run model on it
dd.model <- list(X=test$XX,
                 VCOV=test$tree.vcv,
                 ID=diag(1, nrow=test$nsp, ncol=test$nsp),
                 VV.mu=matrix(0, nrow=2, ncol=test$nsp),
                 z.init=z.init,
                 nsite=test$nsite, 
                 nyr=test$nyr,
                 nrep=nrep,
                 nsp=test$nsp,
                 date.mat=array(0, dim=c(test$nsite,
                                     test$nyr,
                                     max(nrep))),
                 site.presence=site.presence,
                 frac.presence=frac.presence,
                 landcover=landcover,
                 sim.trait=test$trait)

save(dd.model,
     file='multi-species-chase/data/simulated/dd.model.RData')
