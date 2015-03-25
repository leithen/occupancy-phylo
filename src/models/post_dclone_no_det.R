run.model.no.det <- function(d,
                      ni=1100,
                      nt=100,
                      nb=100,
                      nc=3) {
  
  sink('model.jags')
  cat('model {

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
    beta.mat[1:nsp,1:nsp] <- lambda*VCOV[,] + (1-lambda)*ID[,]
    ## convert to precision (or co-precision) matrix - While lambda
    ## above affects the relative weighting on the diagonal in precision
    ## matrix, sigma.psi.beta below affects the overall magnitudes of the
    ## entries and thus the spread of the species-specific responses.
    tau.beta.sp.mat[1:nsp,1:nsp] <- inverse((sigma.psi.beta^2) * beta.mat[,])
    
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
      
      for(site in 1:nsite) {
        for(yr in 1:nyr) {        
          
          logit(psi[site,yr,sp]) <-
            psi.0[sp] + 
              psi.beta.sp[sp]*env[site] +
                psi.site[site] +
                  psi.year[yr]
          
          for(rep in 1:nrep[site,yr,sp]) {
            X[site,yr,rep,sp] ~ dbern(psi[site,yr,sp])
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
  c(
    'psi.site',
    'psi.year',
    'psi.beta.sp')
