ms.ms <- function(d,
                  ni=600,
                  nt=10,
                  nb=100,
                  nc=3) {
  
  sink('model.jags')
  cat('model {

    ## multi-species priors
    mu.p.0  ~ dnorm(0,0.001)
    sigma.p.0 ~ dunif(0,100)
    tau.p.0 <- 1/(sigma.p.0*sigma.p.0)

	## Random effect of species
    for(i in 1:2) {
      mu.phi.0[i] ~ dnorm(0,0.001)
      mu.gam.0[i] ~ dnorm(0,0.001)
      sigma.phi.0[i] ~ dunif(0,100)
      sigma.gam.0[i] ~ dunif(0,100)
      tau.phi.0[i] <- 1/(sigma.phi.0[i]*sigma.phi.0[i])
      tau.gam.0[i] <- 1/(sigma.gam.0[i]*sigma.gam.0[i])
    }

    for(sp in 1:nsp) {
      for(i in 1:2) {
        phi.0[i,sp] ~ dnorm(mu.phi.0[i],tau.phi.0[i])
        gam.0[i,sp] ~ dnorm(mu.gam.0[i],tau.gam.0[i])
      }
	}
	
	##Random effect of site
	tau.phi.site ~ dgamma(1,1)
	tau.gam.site ~ dgamma(1,1)
	
	for (site in 1:nsite) {
		phi.site[site] ~ dnorm(0, tau.phi.site)
		gam.site[site] ~ dnorm(0, tau.gam.site)
	}
			
    ## species-specific detectability
    for(sp in 1:nsp) {
      p.0[sp] ~ dnorm(mu.p.0,tau.p.0)
      logit(p[sp]) <- p.0[sp]
    }
    

## create species mean vectors which are required for
    ## multi-variate normal and then create VCOV matrices, etc
    for(i in 1:2) {
      for(sp in 1:nsp){
      	phi.VV.mu[i,sp]<-phi.0[i,sp]
      	gam.VV.mu[i,sp]<-gam.0[i,sp]
      }
      
      tau.phi.spp[i] ~ dgamma(1,1)
      tau.gam.spp[i] ~ dgamma(1,1)

      phi.lambda[i] <-1
      gam.lambda[i] <-1

      phi.lambdaMat[i,1:nsp,1:nsp] <-
        phi.lambda[i]*VCOV[,] + (1-phi.lambda[i])*ID[,]
      gam.lambdaMat[i,1:nsp,1:nsp] <-
        gam.lambda[i]*VCOV[,] + (1-gam.lambda[i])*ID[,]

      tau.phi.spp.0[i,1:nsp,1:nsp] <-
        tau.phi.spp[i]*inverse(phi.lambdaMat[i,,])
      tau.gam.spp.0[i,1:nsp,1:nsp] <-
        tau.gam.spp[i]*inverse(gam.lambdaMat[i,,])

      phi.VVphylo[i,1:nsp] ~ dmnorm(phi.VV.mu[i,], tau.phi.spp.0[i,,])
      gam.VVphylo[i,1:nsp] ~ dmnorm(phi.VV.mu[i,], tau.gam.spp.0[i,,])
    }
        
    for(sp in 1:nsp) {

      for(site in 1:nsite) {
        
        ## occupancy in year 1
        psi[site,1,sp] <- frac.presence[site,sp]
        Z[site,1,sp] ~ dbern(psi[site,1,sp])

        ## detectability in year 1
        mu.p[site,1,sp] <- Z[site,1,sp]*p[sp]
        for(rep in 1:nrep[site,1,sp]) {
          X[site,1,rep,sp] ~ dbern(mu.p[site,1,sp])
        }

        logit(phi[site,sp]) <- phi.VVphylo[landcover[site],sp]
        								## + 	phi.site[site]
        
        logit(gam[site,sp]) <- gam.VVphylo[landcover[site],sp]
        								##	+	gam.site[site]
       
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
                  'date.mat', 'landcover', 'frac.presence',  'VCOV', 'ID'),
                inits=d$inits,
                parameters.to.save=d$params,
                model.file='model.jags',
                n.chains=3,
                n.thin=10,
                n.iter=600,
                n.burnin=100,
                working.directory=NULL)
}

## specify the parameters to be monitored
get.params <- function()
  c('phi.0',
    'gam.0',
    'mu.phi.0',
    'mu.gam.0',
    'sigma.phi.0',
    'sigma.gam.0',
    'tau.phi.site',
    'tau.gam.site',
    'phi',
    'gam',
    'tau.phi.spp',
    'tau.gam.spp',
    'phi.lambda',
    'gam.lambda',
    'sigma.phi.spp',
    'sigma.gam.spp')
    
    
    