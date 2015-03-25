ms.ms <- function(d,
                  ni=5000,
                  nt=100,
                  nb=1000,
                  nc=3) {
  
  sink('model.jags')
  cat('model {

    ## multi-species priors
    mu.p.0  ~ dnorm(0,0.001)
    sigma.p.0 ~ dunif(0,100)
    tau.p.0 <- 1/(sigma.p.0*sigma.p.0)

    mu.phi.aa ~ dgamma(1,1)
    mu.gam.aa ~ dgamma(1,1)
    sigma.phi.aa ~ dunif(0,100)
    sigma.gam.aa ~ dunif(0,100)
    tau.phi.aa <- 1/(sigma.phi.aa*sigma.phi.aa)
    tau.gam.aa <- 1/(sigma.gam.aa*sigma.gam.aa)

    mu.phi.cc ~ dnorm(0,0.001)
    mu.gam.cc ~ dnorm(0,0.001)
    sigma.phi.cc ~ dunif(0,100)
    sigma.gam.cc ~ dunif(0,100)
    tau.phi.cc <- 1/(sigma.phi.cc*sigma.phi.cc)
    tau.gam.cc <- 1/(sigma.gam.cc*sigma.gam.cc)

    for(sp in 1:nsp) {
      phi.bb[sp] ~ dunif(0,1)
      gam.bb[sp] ~ dunif(0,1)
    }

    ## species-specific detectability
    for(sp in 1:nsp) {
      p.0[sp] ~ dnorm(mu.p.0,tau.p.0)
      logit(p[sp]) <- p.0[sp]
    }
    
    for(sp in 1:nsp) {

      phi.aa[sp] ~ dnorm(mu.phi.aa,tau.phi.aa)
      gam.aa[sp] ~ dnorm(mu.gam.aa,tau.gam.aa)
      phi.cc[sp] ~ dnorm(mu.phi.cc,tau.phi.cc)
      gam.cc[sp] ~ dnorm(mu.gam.cc,tau.gam.cc)

      for(site in 1:nsite) {
        
        ## occupancy in year 1

        mu.psi.1[site,1,sp] ~ dunif(0,1)
        psi[site,1,sp] <- mu.psi.1[site,1,sp] *
          site.presence[site,sp]
        Z[site,1,sp] ~ dbern(psi[site,1,sp])

        ## detectability in year 1
        mu.p[site,1,sp] <- Z[site,1,sp]*p[sp]
        for(rep in 1:nrep[site,1,sp]) {
          X[site,1,rep,sp] ~ dbern(mu.p[site,1,sp])
        }

        logit(phi[site,sp]) <-
          (-phi.aa[sp])*(landcover[site]-phi.bb[sp])^2 + phi.cc[sp]

        logit(gam[site,sp]) <-
          (-gam.aa[sp])*(landcover[site]-gam.bb[sp])^2 + gam.cc[sp]

        ## occupancy and detectability in subsequent years
        for(yr in 1:(nyear-1)) {

          psi[site,yr+1,sp] <-
            Z[site,yr,sp] * phi[site,sp] + 
              (1-Z[site,yr,sp]) * gam[site,sp]

          Z[site,yr+1,sp] ~ dbern(psi[site,yr+1,sp])

          mu.p[site,yr+1,sp] <- Z[site,yr+1,sp]*p[sp]
          for(rep in 1:nrep[site,yr+1,sp]) {
            X[site,yr+1,rep,sp] ~ dbern(mu.p[site,yr+1,sp])
          }
        }
      }
    }

    for(site in 1:nsite) {
      for(sp in 1:nsp) {
        psi.mean[site,sp] <- mean(psi[site,,sp])
      }
    }
  }', fill = TRUE)
  sink()

  attach(d$data)
  jags.parallel(data=list('X', 'nsp', 'nsite', 'nyear', 'nrep',
                  'date.mat', 'landcover', 'site.presence'),
                inits=d$inits,
                parameters.to.save=d$params,
                model.file='model.jags',
                n.chains=3,
                n.thin=100,
                n.iter=10100,
                n.burnin=100,
                working.directory=NULL)
}

## specify the parameters to be monitored
get.params <- function()
  c('mu.p.0',
    'tau.p.0',
    'phi.bb',
    'gam.bb',
    'mu.phi.aa',
    'tau.phi.aa',
    'mu.gam.aa',
    'tau.gam.aa',
    'mu.phi.cc',
    'tau.phi.cc',
    'mu.gam.cc',
    'tau.gam.cc',
    'phi',
    'gam',
    'psi.mean')
