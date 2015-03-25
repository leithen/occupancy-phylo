run.dclones.model <- function() {

  ## priors and direct derivatives  *** not to be cloned ***

  ## detectability
  mu.p.0 ~ dnorm(0,0.001)
  sigma.p.0 ~ dunif(0,100)
  tau.p.0 <- 1/(sigma.p.0*sigma.p.0)

  ## random effect of site on occupancy
  sigma.psi.site ~ dunif(0,100)
  tau.psi.site <- 1/(sigma.psi.site*sigma.psi.site)

  ## random effect of year on occupancy
  sigma.psi.year ~ dunif(0,100)
  tau.psi.year <- 1/(sigma.psi.year*sigma.psi.year)

  ## fixed effect of species-level covariate on mean occupancy
  ## psi.0.trait ~ dnorm(0,0.001)
  psi.0.trait <- 0
  ## species-specific slope expectation
  mu.psi.beta ~ dnorm(0,0.001)
  ## effect of species-level covariate on slope
  ## beta.trait ~ dnorm(0,0.001)
  beta.trait <- 0
  
  for(sp in 1:nsp) {
    ## construct species-specific occupancy intercepts which include
    ## both a random species effect and fixed effect based on the
    ## species trait value
    psi.0.sp[sp] ~ dnorm(0,0.001)
    psi.0[sp] <- psi.0.sp[sp] + psi.0.trait*trait[sp]
    ## species-specific slope expectation terms
    mu.beta.sp[sp] <- mu.psi.beta + beta.trait*trait[sp]
  }
  
  ## lambda allows us to scale the phylogenetic correlation matrix
  lambda ~ dunif(0,1)

  ## incorporate phylogenetic covariance structure into
  ## species-specific slope expectations.
  beta.mat[1:nsp,1:nsp] <- lambda*VCOV[,] + (1-lambda)*ID[,]
  ## convert to precision (or co-precision) matrix - While lambda
  ## above affects the relative weighting on the diagonal in precision
  ## matrix, sigma.psi.beta below affects the overall magnitudes of the
  ## entries and thus the spread of the species-specific responses.
  sigma.psi.beta ~ dunif(0.01,100)
  tau.beta.sp.mat[1:nsp,1:nsp] <- inverse(sigma.psi.beta * beta.mat[,])

  ## construct likelihood function *** cloned ***
  for(clone in 1:clones.index) {
    
    ## random effect of site on occupancy
    for(site in 1:nsite) {
      psi.site[site,clone] ~ dnorm(0, tau.psi.site)
    }

    ## random effect of year on occupancy
    for(yr in 1:nyr) {
      psi.year[yr,clone] ~ dnorm(0, tau.psi.year)
    }
    
    ## species-specific slopes with regard to the environment
    ## (phylogeny is incorporated here)
    psi.beta.sp[1:nsp,clone] ~ dmnorm(mu.beta.sp[],
                                      tau.beta.sp.mat[,])


    for(sp in 1:nsp) {
    
      ## species-specific detectability
      p.0[sp,clone] ~ dnorm(mu.p.0, tau.p.0)
      logit(p[sp,clone]) <- p.0[sp,clone]
      
      for(site in 1:nsite) {
        for(yr in 1:nyr) {        
          
          logit(psi[site,yr,sp,clone]) <-
            psi.0[sp] + 
              psi.beta.sp[sp,clone]*env[site] +
                psi.site[site,clone] +
                  psi.year[yr,clone]
          
          Z[site,yr,sp,clone] ~ dbern(psi[site,yr,sp,clone])
          mu.p[site,yr,sp,clone] <- Z[site,yr,sp,clone]*p[sp,clone]
          for(rep in 1:nrep[site,yr,sp]) {
            X[site,yr,rep,sp,clone] ~ dbern(mu.p[site,yr,sp,clone])
          }
        }
      }
    }
  }
}

## specify the parameters to be monitored
get.params <- function()
  c('mu.p.0',
    'sigma.p.0',
    'sigma.psi.site',
    'sigma.psi.year',
    'psi.0.sp',
    'mu.psi.beta',
    'sigma.psi.beta',
    'lambda')
