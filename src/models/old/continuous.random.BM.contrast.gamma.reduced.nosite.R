ms.ms <- function(d,
                  ni=3100,
                  nt=50,
                  nb=100,
                  nc=3) {
  
  sink('model.jags')
  cat('model {

    ## multi-species priors
    mu.p.0  ~ dnorm(0,0.001)
    sigma.p.0 ~ dunif(0,100)
    tau.p.0 <- 1/(sigma.p.0*sigma.p.0)


	## Random intercept of species
      mu.phi.0 ~ dnorm(0,0.001)
      mu.gam.0 ~ dnorm(0,0.001)
      sigma.phi.0 ~ dunif(0,100)
      sigma.gam.0 ~ dunif(0,100)
      tau.phi.0 <- 1/(sigma.phi.0*sigma.phi.0)
      tau.gam.0 <- 1/(sigma.gam.0*sigma.gam.0)
    

    for(sp in 1:nsp) {
        phi.0[sp] ~ dnorm(mu.phi.0,tau.phi.0)
        gam.0[sp] ~ dnorm(mu.gam.0,tau.gam.0)
      }
	
	
	## Fixed landcover slope across species
     phi.landcover ~ dnorm(0, 0.0001)
     gam.landcover ~ dnorm(0, 0.0001)
		
	
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
      
      phi.VV.mu.0 <- 0
      gam.VV.mu.0 <- 0
      
      for (sp in 1:nsp) {
      phi.VV.mu[sp]<-phi.VV.mu.0
      gam.VV.mu[sp]<-gam.VV.mu.0
      }
      
      tau.phi.spp ~ dgamma(1,1)
      tau.gam.spp ~ dgamma(1,1)
      phi.lambda <-1
      gam.lambda <-1

      phi.lambdaMat[1:nsp,1:nsp] <-
       			 phi.lambda*VCOV[,] + (1-phi.lambda)*ID[,]
      
      gam.lambdaMat[1:nsp,1:nsp] <-
       			 gam.lambda*VCOV[,] + (1-gam.lambda)*ID[,]

      tau.phi.spp.0[1:nsp,1:nsp] <-
      			  tau.phi.spp*inverse(phi.lambdaMat[,])
      
      tau.gam.spp.0[1:nsp,1:nsp] <-
        			tau.gam.spp*inverse(gam.lambdaMat[,])

      phi.VVphylo[1:nsp] ~ dmnorm(phi.VV.mu[], tau.phi.spp.0[,])
      gam.VVphylo[1:nsp] ~ dmnorm(gam.VV.mu[], tau.gam.spp.0[,])

    
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

        logit(phi[site,sp]) <- phi.0[sp]
        									+ phi.landcover*landcover[site]
        										+ phi.VVphylo[sp]*landcover[site]
        		##									+ phi.site[site]
        logit(gam[site,sp]) <- gam.0[sp]
        									+ gam.landcover*landcover[site]
        										+ gam.VVphylo[sp]*landcover[site]
        		##									+ gam.site[site]
        
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
                n.thin=50,
                n.iter=3100,
                n.burnin=100,
                working.directory=NULL)
}

## specify the parameters to be monitored
get.params <- function()
  c('phi.0',
    'gam.0',
    'mu.phi.0',
    'mu.gam.0',
    'phi.landcover',
    'gam.landcover',
    'sigma.phi.0',
    'sigma.gam.0',
    'sigma.phi.landcover',
    'sigma.gam.landcover',
    'tau.phi.site',
    'tau.gam.site',
    'phi',
    'gam',
    'psi.mean',
    'N',
    'tau.phi.spp',
    'tau.gam.spp',
    'phi.lambda',
    'gam.lambda',
    'sigma.phi.spp',
    'sigma.gam.spp')
    
    
    
    