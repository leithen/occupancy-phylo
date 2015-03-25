ms.ms <- function(d,
                  ni=1200,
                  nt=10,
                  nb=200,
                  nc=3) {
  
  sink('model.jags')
  cat('model {

    ## multi-species priors
    mu.p.0  ~ dnorm(0,0.001)
    zsig2.p.0~ dunif(0,100)
    tau.p.0 <- 1/zsig2.p.0

  ## species-specific detectability
    for(sp in 1:nsp) {
      p.0[sp] ~ dnorm(mu.p.0,tau.p.0)
      logit(p[sp]) <- p.0[sp]
    }



	## Fixed intercept of species

    for(sp in 1:nsp) {
        psi.0[sp] ~ dnorm(0, 0.0001)
      }
	
	zsig2.psi.0<-sd(psi.0)*sd(psi.0)
	
		
	##Random effect of site

	zsig2.psi.site  ~ dunif(0,100)	
	tau.psi.site <- 1/zsig2.psi.site
	
	for (site in 1:nsite) {
		psi.site[site] ~ dnorm(0, tau.psi.site)
	}
		
	##Random effect of year
	zsig2.psi.year ~	dunif(0,100)
	tau.psi.year <- 1/zsig2.psi.year
	
	for (yr in 1:nyear) {
		psi.year[yr] ~ dnorm(0, tau.psi.year)
	}

    ## create species mean vectors which are required for
    ## multi-variate normal and then create VCOV matrices, etc
      
      psi.VV.mu.zero <- 0
      
      for (sp in 1:nsp) {
     	 psi.VV.mu[sp]<-psi.VV.mu.zero
      }
     
      zsig2.psi.phylo ~ dunif(0,100)     
      tau.psi.spp <- 1/zsig2.psi.phylo
      
      psi.lambda <-1

      psi.lambdaMat[1:nsp,1:nsp] <-
       			 psi.lambda*VCOV[,] + (1-psi.lambda)*ID[,]
      
      tau.psi.spp.VV[1:nsp,1:nsp] <-
      			  tau.psi.spp*inverse(psi.lambdaMat[,])

      psi.VVphylo[1:nsp] ~ dmnorm(psi.VV.mu[], tau.psi.spp.VV[,])

	##Create tracking variables for the proportion of variance in species occupancy that comes from phylogeney and other sources
	 
	 Prop.Phylo.Spp<-zsig2.psi.phylo/(zsig2.psi.0+zsig2.psi.phylo)	
	 Prop.Phylo<-zsig2.psi.phylo/(zsig2.psi.0+zsig2.psi.phylo+ zsig2.psi.site+zsig2.psi.year)
	 Prop.Spp<-zsig2.psi.0/(zsig2.psi.0+zsig2.psi.phylo+ zsig2.psi.site+zsig2.psi.year)
	 Prop.Site<-zsig2.psi.site/(zsig2.psi.0+zsig2.psi.phylo+ zsig2.psi.site+zsig2.psi.year)
	 Prop.Year<-zsig2.psi.year/(zsig2.psi.0+zsig2.psi.phylo+ zsig2.psi.site+zsig2.psi.year)
    
    for(sp in 1:nsp) {
    	
      for(site in 1:nsite) {
		
			for(yr in 1:nyear) {        
       
          logit(psi[site,yr,sp]) <-psi.0[sp] + 
        										psi.VVphylo[sp] +
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
                n.thin=10,
                n.iter=1200,
                n.burnin=200,
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
    'zsig2.psi.phylo',
    'zsig2.psi.site',
    'zsig2.psi.year',
    'Prop.Phylo',
    'Prop.Site',
    'Prop.Year',
    'Prop.Spp',
    'Prop.Phylo.Spp')
    
    
    
    