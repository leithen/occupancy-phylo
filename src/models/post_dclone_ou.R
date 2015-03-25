run.model <- function(d,
                      ni=1100,
                      nt=100,
                      nb=100,
                      nc=3) {
  
  sink('model.jags')
  cat('model {

    tau.p.0 <- 1/(sigma.p.0*sigma.p.0)
    tau.psi.site <- 1/(sigma.psi.site*sigma.psi.site)
    tau.psi.year <- 1/(sigma.psi.year*sigma.psi.year)

    psi.0.trait <- 0
    beta.trait <- 0
    
    for(sp in 1:nsp) {
      ## construct species-specific occupancy intercepts which include
      ## both a random species effect and fixed effect based on the
      ## species trait value
      psi.0[sp] <- psi.0.sp[sp] + psi.0.trait*trait[sp]
      ## species-specific slope expectation terms
      mu.beta.sp[sp] <- mu.psi.beta + beta.trait*trait[sp]
    }

      ## incorporate phylogenetic covariance structure into
      ## species-specific slope expectations.
  for (i in 1:nsp) {
  	for (j in 1:nsp) {
          beta.mat[i,j] <- ((sigma.psi.beta^2)/(2*alpha))*exp(-2*alpha*(1-VCOV[i,j]))*(1-exp(-2*alpha*VCOV[i,j])) 
		}
	}
 
  ## convert to precision (or co-precision) matrix
  tau.beta.sp.mat[1:nsp,1:nsp] <- inverse(beta.mat[,])    

    ## random effect of site on occupancy
    for(site in 1:nsite) {
      psi.site[site] ~ dnorm(0, tau.psi.site)
    }

    ## random effect of year on occupancy
    for(yr in 1:nyr) {
      psi.year[yr] ~ dnorm(0, tau.psi.year)
    }
    
    ## species-specific slopes with regard to the environment
    ## (phylogeny is incorporated here)
    psi.beta.sp[1:nsp] ~ dmnorm(mu.beta.sp[],
                                tau.beta.sp.mat[,])


    for(sp in 1:nsp) {
      
      ## species-specific detectability
      p.0[sp] ~ dnorm(mu.p.0, tau.p.0)
      logit(p[sp]) <- p.0[sp]
      
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

  jags(data=d$data,
       d$inits,
       d$params,
       'model.jags',
       n.chains=nc,
       n.thin=nt,
       n.iter=ni,
       n.burnin=nb,
       working.directory=NULL)
}

## specify the parameters to be monitored
get.params <- function()
  c('p.0',
    ## 'psi.0',
    'psi.site',
    'psi.year',
    'psi.beta.sp')
