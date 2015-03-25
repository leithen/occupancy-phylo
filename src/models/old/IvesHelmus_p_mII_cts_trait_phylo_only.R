ms.ms <- function(d,
                  ni=1100,
                  nt=50,
                  nb=100,
                  nc=3) {
  
  sink('model.jags')
  cat('model {

    ## multi-species priors
    mu.p.0  ~ dnorm(0,0.001)
    sigma.p.0 ~ dunif(0,100)
    tau.p.0 <- 1/(sigma.p.0*sigma.p.0)

    ## species-specific detectability
    for(sp in 1:nsp) {
      p.0[sp] ~ dnorm(mu.p.0,tau.p.0)
      logit(p[sp]) <- p.0[sp]
    }

    ## random intercept of species
    mu.psi.0 ~ dnorm(0,0.001)
    sigma.psi.0 ~ dunif(0,100)
    tau.psi.0 <- 1/(sigma.psi.0*sigma.psi.0)

    ## effect of species-level covariate on mean occupancy
    psi.0.trait ~ dnorm(0,0.001)
    for(sp in 1:nsp) {
      psi.0.mean[sp] ~ dnorm(mu.psi.0,tau.psi.0)
      psi.0[sp] <- psi.0.mean[sp] + psi.0.trait*trait.value[sp]
    }
    
    ## random effect of site
    mu.psi.site ~ dnorm(0,0.001)
    sigma.psi.site ~ dunif(0,100)
    tau.psi.site <- 1/(sigma.psi.site*sigma.psi.site)
    for (site in 1:nsite) {
      psi.site[site] ~ dnorm(0, tau.psi.site)
    }

    ## random effect of year
    mu.psi.year ~ dnorm(0,0.001)
    sigma.psi.year ~ dunif(0,100)
    tau.psi.year <- 1/(sigma.psi.year*sigma.psi.year)
    for (yr in 1:nyr) {
      psi.year[yr] ~ dnorm(0, tau.psi.year)
    }

    ## intercept for model on slope
    mu.beta.landcover ~ dnorm(0,0.001)

    ## create species mean vectors which are required for
    ## multi-variate normal and then create VCOV matrices, etc
    beta.lambda <- 1
    beta.VV.mu.zero <- 0
    for (sp in 1:nsp) {
      beta.VV.mu[sp] <- beta.VV.mu.zero
    }

    sigma.beta.spp ~ dunif(0,100)
    tau.beta.spp <- 1/(sigma.beta.spp*sigma.beta.spp)
    beta.lambdaMat[1:nsp,1:nsp] <-
      beta.lambda*VCOV[,] + (1-beta.lambda)*ID[,]
    tau.beta.spp.VV[1:nsp,1:nsp] <-
      tau.beta.spp*inverse(beta.lambdaMat[,])
    beta.VVphylo[1:nsp] ~ dmnorm(beta.VV.mu[], tau.beta.spp.VV[,])

    ## effect of species-level covariate on slope
    beta.trait ~ dnorm(0,0.001)
    
    for(sp in 1:nsp) {

      ## combine species specific portion of slope with phylogenetic
      ## portion of slope to get overall slope
      psi.landcover[sp] <-
        mu.beta.landcover +
          beta.VVphylo[sp] +
            beta.trait*trait.value[sp]
      
      for(site in 1:nsite) {
        for(yr in 1:nyr) {        
          
          logit(psi[site,yr,sp]) <-
            psi.0[sp] +
              psi.site[site] +
                psi.year[yr] + 
                  psi.landcover[sp]*landcover[site]

          Z[site,yr,sp] ~ dbern(psi[site,yr,sp])
          mu.p[site,yr,sp] <- Z[site,yr,sp]*p[sp]
          for(rep in 1:nrep[site,yr,sp]) {
            X[site,yr,rep,sp] ~ dbern(mu.p[site,yr,sp])
          }
        }
      }
    }
    
    zsig2.p.0                  <- 1/tau.p.0
    zsig2.psi.0                <- 1/tau.psi.0
    zsig2.psi.site             <- 1/tau.psi.site
    zsig2.psi.year             <- 1/tau.psi.year
    zsig2.beta.landcover.phylo <- 1/tau.beta.spp
    
  }', fill = TRUE)
  sink()

  attach(d$data)
  res <- jags.parallel(data=list('X', 'nsp', 'nsite', 'nyr', 'nrep', 
                         'landcover', 'VCOV', 'ID', 'trait.value'),
                       inits=d$inits,
                       parameters.to.save=d$params,
                       model.file='model.jags',
                       n.chains=3,
                       n.thin=50,
                       n.iter=11000,
                       n.burnin=1000,
                       working.directory=NULL)
  detach(d$data)
  res
}

## specify the parameters to be monitored
get.params <- function()
  c('psi.0',
    'mu.psi.0',
    'mu.beta.landcover',
    'psi.lambda',
    'mu.p.0',
    'zsig2.p.0',
    'zsig2.psi.0',
    'zsig2.psi.site',
    'zsig2.psi.year',
    'zsig2.beta.landcover.phylo',
    'psi.landcover',
    'beta.VVphylo',
    'beta.0.trait',
    'beta.trait')
