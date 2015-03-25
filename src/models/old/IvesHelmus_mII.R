ms.ms <- function(d,
                  ni=1000,
                  nt=10,
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
      }
	
	
	## Random landcover slope of species
      mu.phi.landcover ~ dnorm(0,0.001)
      sigma.phi.landcover ~ dunif(0,100)
      tau.phi.landcover <- 1/(sigma.phi.landcover*sigma.phi.landcover)
    

    for(sp in 1:nsp) {
        phi.landcover[sp] ~ dnorm(mu.phi.landcover,tau.phi.landcover)
      }
		
	
	##Random effect of site
	tau.phi.site ~ dgamma(1,1)
	
	for (site in 1:nsite) {
		phi.site[site] ~ dnorm(0, tau.phi.site)
	}

	

    ## create species mean vectors which are required for
    ## multi-variate normal and then create VCOV matrices, etc
      
      phi.VV.mu.0 <- 0
      
      for (sp in 1:nsp) {
     	 phi.VV.mu[sp]<-phi.VV.mu.0
      }
      
      tau.phi.spp ~ dgamma(1,1)
      phi.lambda <-1


      phi.lambdaMat[1:nsp,1:nsp] <-
       			 phi.lambda*VCOV[,] + (1-phi.lambda)*ID[,]
      
      tau.phi.spp.0[1:nsp,1:nsp] <-
      			  tau.phi.spp*inverse(phi.lambdaMat[,])

      phi.VVphylo[1:nsp] ~ dmnorm(phi.VV.mu[], tau.phi.spp.0[,])

    
    for(sp in 1:nsp) {
    	
      for(site in 1:nsite) {
		
			for(yr in 1:nyear) {        

				for(rep in 1:6){        
					
          logit(psi[site,yr,rep,sp]) <-phi.0[sp] + 
        									phi.landcover[sp]*landcover[site] +
        										phi.VVphylo[sp]*landcover[site] +
        											phi.site[site]
          
          X[site,yr,rep,sp] ~ dbern(psi[site,yr,rep,sp])
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
                n.thin=10,
                n.iter=1000,
                n.burnin=100,
                working.directory=NULL)
}

## specify the parameters to be monitored
get.params <- function()
  c('phi.0',
    'gam.0',
    'mu.phi.0',
    'mu.gam.0',
    'mu.phi.landcover',
    'mu.gam.landcover',
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
    
    
    
    