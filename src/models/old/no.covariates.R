ms.ms <- function(d,
                  ni=1100,
                  nt=100,
                  nb=100,
                  nc=3) {
  
  sink('model.jags')
  cat('model {

    ## multi-species priors
    mu.p.0  ~ dnorm(0,0.001)
    sigma.p.0 ~ dunif(0,100)
    tau.p.0 <- 1/(sigma.p.0*sigma.p.0)

    mu.phi.0 ~ dnorm(0,0.001)
    mu.gam.0 ~ dnorm(0,0.001)
    sigma.phi.0 ~ dunif(0,100)
    sigma.gam.0 ~ dunif(0,100)
    tau.phi.0 <- 1/(sigma.phi.0*sigma.phi.0)
    tau.gam.0 <- 1/(sigma.gam.0*sigma.gam.0)

    phi.landcover ~ dnorm(0,0.001)
    gam.landcover ~ dnorm(0,0.001)

    ## species-specific detectability
    for(sp in 1:nsp) {
      p.0[sp] ~ dnorm(mu.p.0,tau.p.0)
      logit(p[sp]) <- p.0[sp]
    }
    
    for(sp in 1:nsp) {

      phi.0[sp] ~ dnorm(mu.phi.0,tau.phi.0)
      gam.0[sp] ~ dnorm(mu.gam.0,tau.gam.0)

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
          phi.0[sp] + phi.landcover*landcover[site]
        
        logit(gam[site,sp]) <-
          gam.0[sp] + gam.landcover*landcover[site]

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
      for(yr in 1:nyear) {
        N[site,yr] <- sum(Z[site,yr,1:nsp])
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
    'mu.phi.0',
    'mu.gam.0',
    'tau.phi.0',
    'tau.gam.0',
    'phi',
    'gam',
    'phi.landcover',
    'gam.landcover',
    'psi.mean',
    'N')
