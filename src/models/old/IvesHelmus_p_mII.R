ms.ms <- function(d,
                  ni=11000,
                  nt=50,
                  nb=1000,
                  nc=3) {
  
  sink('model.jags')
  cat('model {

    ## multi-species priors
    mu.p.0  ~ dnorm(0,0.001)
    sigma.p.0 ~ dunif(0,100)
    zsig2.p.0 <- sigma.p.0*sigma.p.0
    tau.p.0 <- 1/(sigma.p.0*sigma.p.0)

    ## species-specific detectability
    for(sp in 1:nsp) {
      p.0[sp] ~ dnorm(mu.p.0,tau.p.0)
      logit(p[sp]) <- p.0[sp]
    }

    ## Random intercept of species
    mu.psi.0 ~ dnorm(0,0.001)
    sigma.psi.0 ~ dunif(0,100)
    zsig2.psi.0 <- sigma.psi.0*sigma.psi.0
    tau.psi.0 <- 1/(sigma.psi.0*sigma.psi.0)

    for(sp in 1:nsp) {
      psi.0[sp] ~ dnorm(mu.psi.0,tau.psi.0)
    }
    
    ## Random landcover slope of species
    mu.psi.landcover ~ dnorm(0,0.001)
    sigma.psi.landcover ~ dunif(0,100)
    zsig2.psi.landcover <- sigma.psi.landcover* sigma.psi.landcover
    tau.psi.landcover <- 1/(sigma.psi.landcover*sigma.psi.landcover)

    for(sp in 1:nsp) {
      psi.landcover[sp] ~ dnorm(0, tau.psi.landcover)
    }
    
    ## Random effect of site
    tau.psi.site ~ dgamma(1,1)
    zsig2.psi.site<-1/tau.psi.site
    
    for (site in 1:nsite) {
      psi.site[site] ~ dnorm(0, tau.psi.site)
    }

    ## Random effect of year
    tau.psi.year ~ dgamma(1,1)
    zsig2.psi.year <- 1/tau.psi.year
    
    for (yr in 1:nyear) {
      psi.year[yr] ~ dnorm(0, tau.psi.year)
    }

    ## create species mean vectors which are required for
    ## multi-variate normal and then create VCOV matrices, etc
    
    psi.VV.mu.zero <- 0
    
    for (sp in 1:nsp) {
      psi.VV.mu[sp]<-psi.VV.mu.zero
    }
    
    tau.psi.spp ~ dgamma(1,1)
    zsig2.psi.landcover.phylo<-1/tau.psi.spp
    
    psi.lambda <-1

    psi.lambdaMat[1:nsp,1:nsp] <-
      psi.lambda*VCOV[,] + (1-psi.lambda)*ID[,]
    
    tau.psi.spp.VV[1:nsp,1:nsp] <-
      tau.psi.spp*inverse(psi.lambdaMat[,])

    psi.VVphylo[1:nsp] ~ dmnorm(psi.VV.mu[], tau.psi.spp.VV[,])

    ## Create tracking variable for the proportion of variance in
    ## response to landcover that comes from phylogeney 
    
    Prop.Phylo <- zsig2.psi.landcover.phylo / (zsig2.psi.landcover +
                                               zsig2.psi.landcover.phylo)

    ## Convert species specific portion of slope with phylogenetic
    ## portion of slope to get overall slope 
    
    beta.landcover[1:nsp] <-
      mu.psi.landcover +
        psi.landcover[1:nsp] +
          psi.VVphylo[1:nsp]
    
    for(sp in 1:nsp) {
      for(site in 1:nsite) {
        for(yr in 1:nyear) {        
          
          logit(psi[site,yr,sp]) <- psi.0[sp] + 
            beta.landcover[sp]*landcover[site] +
              psi.site[site] +
                psi.year[yr]

          Z[site,yr,sp] ~ dbern(psi[site,yr,sp])

          mu.p[site,yr,sp] <- Z[site,yr,sp]*p[sp]
          
          for(rep in 1:nrep[site,yr,sp]) {
            X[site,yr,rep,sp] ~ dbern(mu.p[site,yr,sp])
          }

          
        }
      }
    }
    
    
  }', fill = TRUE)
  sink()

  attach(d$data)
  jags.parallel(data=list('X', 'nsp', 'nsite', 'nyear', 'nrep',
                  'date.mat', 'landcover', 'frac.presence',  'VCOV', 'ID'),
                inits=d$inits,
                parameters.to.save=d$params,
                model.file='model.jags',
                n.chains=3,
                n.thin=50,
                n.iter=11000,
                n.burnin=1000,
                working.directory=NULL)
}

## specify the parameters to be monitored
get.params <- function()
  c('psi.0',
    'mu.psi.0',
    'mu.psi.landcover',
                                        #    'psi',
    'psi.lambda',
    'mu.p.0',
    'zsig2.p.0',
    'zsig2.psi.0',
    'zsig2.psi.landcover',
    'zsig2.psi.landcover.phylo',
    'zsig2.psi.site',
    'zsig2.psi.year',
    'Prop.Phylo',
    'beta.landcover',
    'psi.landcover',
    'psi.VVphylo')



