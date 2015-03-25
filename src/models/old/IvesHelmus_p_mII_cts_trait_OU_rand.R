ms.ms <- function(d,
                  ni=11000,
                  nt=100,
                  nb=1000,
                  nc=3) {
  
  sink('model.jags')
  cat('model {

    ## detectability
    mu.p.0  ~ dnorm(0,0.001)
    sigma.p.0 ~ dunif(0,100)
    tau.p.0 <- 1/(sigma.p.0*sigma.p.0)

    ## species-specific detectability
    for(sp in 1:nsp) {
      p.0[sp] ~ dnorm(mu.p.0,tau.p.0)
      logit(p[sp]) <- p.0[sp]
    }
    
    ## random effect of site on occupancy
    mu.psi.site ~ dnorm(0,0.001)
    sigma.psi.site ~ dunif(0,100)
    tau.psi.site <- 1/(sigma.psi.site*sigma.psi.site)
    for (site in 1:nsite) {
      psi.site[site] ~ dnorm(0, tau.psi.site)
    }

    ## random effect of year on occupancy
    mu.psi.year ~ dnorm(0,0.001)
    sigma.psi.year ~ dunif(0,100)
    tau.psi.year <- 1/(sigma.psi.year*sigma.psi.year)
    for (yr in 1:nyr) {
      psi.year[yr] ~ dnorm(0, tau.psi.year)
    }

    ## construct species-specific occupancy intercepts which include
    ## both a random species effect and fixed effect based on the
    ## species trait value
    ##
    ## ## species-specific random effect in intercept
    ## mu.psi.0 ~ dnorm(0,0.001)
    ## sigma.psi.0 ~ dunif(0,100)
    ## tau.psi.0 <- 1/(sigma.psi.0*sigma.psi.0)

    psi.0.mean ~ dnorm(0,0.001)
    ## fixed effect of species-level covariate on mean occupancy
    psi.0.trait ~ dnorm(0,0.001)
    
    for(sp in 1:nsp) {
      ## psi.0.mean[sp] ~ dnorm(mu.psi.0,tau.psi.0)
      ## psi.0[sp] <- psi.0.mean[sp] + psi.0.trait*trait[sp]
      psi.0[sp] <- psi.0.mean + psi.0.trait*trait[sp]
    }

    ## intercept for model on slope
    mu.beta ~ dnorm(0,0.001)
    
    ## create species-specific slopes with regard to the environment
    ## (phylogeny is incorporated here)
    ##
    ## species-specific slope expectation and precision terms
    ## (beta.trait is the effect of species-level covariate on slope)
    beta.trait ~ dnorm(0,0.001)
    for(sp in 1:nsp) {
      mu.beta.sp[sp] <- mu.beta + beta.trait*trait[sp]
    }
    sigma.beta.sp ~ dunif(0,100)
   
    ## incorporate phylogenetic covariance structure into
    ## species-specific slope expectations.  alpha.ou allows us to scale
    ## the phylogenetic correlation matrix
    alpha.ou ~ dunif(0,100)
    
    for (i in 1:nsp){
    	for (j in 1:nsp){
    		 beta.mat[i,j] <- ((sigma.beta.sp^2)/(2*alpha.ou))*exp(-2*alpha.ou*(1-VCOV[i,j]))*(1-exp(-2*alpha.ou*VCOV[i,j]))
    	}
    }
   
  
    ## convert to precision (or co-precision) matrix - While alpha.ou
    ## above affects the relative weighting on the diagonal in
    ## precision matrix, tau.beta.sp below affects the overall
    ## magnitudes of the entries and thus the spread of the
    ## species-specific responses.
    tau.beta.sp.mat[1:nsp,1:nsp] <- inverse(beta.mat[,])
    ## draw species-specific slopes
    psi.beta.sp[1:nsp] ~ dmnorm(mu.beta.sp[], tau.beta.sp.mat[,])
    
    for(sp in 1:nsp) {
      for(site in 1:nsite) {
        for(yr in 1:nyr) {        
          
          logit(psi[site,yr,sp]) <-
            psi.0[sp] + 
              psi.beta.sp[sp]*env[site] +
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
  res <- jags.parallel(data=list('X', 'nsp', 'nsite', 'nyr', 'nrep', 
                         'env', 'VCOV', 'ID', 'trait'),
                       inits=d$inits,
                       parameters.to.save=d$params,
                       model.file='model.jags',
                       n.chains=3,
                       n.thin=100,
                       n.iter=11000,
                       n.burnin=1000,
                       working.directory=NULL)
  detach(d$data)
  res
}

## specify the parameters to be monitored
get.params <- function()
  c('alpha.ou',
    'mu.psi.0',
    'mu.beta',
    'mu.p.0',
    'sigma.p.0',
    ## 'sigma.psi.0',
    'sigma.psi.site',
    'sigma.psi.year',
    'sigma.beta.sp',
    'psi.0',
    'psi.0.mean',
    'psi.0.trait',
    'psi.beta.sp',
    'beta.trait')
