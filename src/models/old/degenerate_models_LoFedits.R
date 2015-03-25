## specify the parameters to be monitored
get.params <- function()
  c('lambda',
    'sigma.beta.sp')

run.dclones.model <- function() {
  
  ## species-specific intercepts and the effect of the environment are
  ## modeled as fixed effects

  ## THIS GIVES COMPARABLE RESULTS
  for(sp in 1:nsp) {
    psi.0[sp] ~ dnorm(0,0.001)
    mu.beta.sp[sp] ~ dnorm(0,0.001)
  }

  ## ## THIS DOES NOT
  ## mu.beta ~ dnorm(0,0.001)
  ## for(sp in 1:nsp) {
  ##   psi.0[sp] ~ dnorm(0,0.001)
  ##   mu.beta.sp[sp] <- mu.beta
  ## }

  
  ## incorporate phylogenetic covariance structure into
  ## species-specific slope expectations.  lambda allows us to scale
  ## the phylogenetic correlation matrix
  lambda ~ dunif(0,1)
  beta.mat[1:nsp,1:nsp] <- lambda*VCOV[,] + (1-lambda)*ID[,]
  ## convert to precision (or co-precision) matrix
  sigma.beta.sp ~ dunif(0.01,100)
  tau.beta.sp.mat[1:nsp,1:nsp] <-  inverse(sigma.beta.sp *beta.mat[,])
  ## draw species-specific slopes
  psi.beta.sp[1:nsp] ~ dmnorm(mu.beta.sp[], tau.beta.sp.mat[,])
  
  for(sp in 1:nsp) {
    for(site in 1:nsite) {
      for(yr in 1:nyr) {        
        
        logit(psi[site,yr,sp]) <-
          psi.0[sp] + 
            psi.beta.sp[sp]*env[site]
        Z[site,yr,sp] ~ dbern(psi[site,yr,sp])
      }
    }
  }
}

run.R2jags.model <- function(d) {
  
  sink('model.jags')
  cat('model {

    ## species-specific intercepts and the effect of the environment are
    ## modeled as fixed effects

    ## THIS GIVES COMPARABLE RESULTS
    for(sp in 1:nsp) {
      psi.0[sp] ~ dnorm(0,0.001)
      mu.beta.sp[sp] ~ dnorm(0,0.001)
    }

    ## ## THIS DOES NOT
    ## mu.beta ~ dnorm(0,0.001)
    ## for(sp in 1:nsp) {
    ##   psi.0[sp] ~ dnorm(0,0.001)
    ##   mu.beta.sp[sp] <- mu.beta
    ## }
    
    ## incorporate phylogenetic covariance structure into
    ## species-specific slope expectations.  lambda allows us to scale
    ## the phylogenetic correlation matrix
    lambda ~ dunif(0,1)
    beta.mat[1:nsp,1:nsp] <- lambda*VCOV[,] + (1-lambda)*ID[,]
    ## convert to precision (or co-precision) matrix
    sigma.beta.sp ~ dunif(0.01,100)
    tau.beta.sp.mat[1:nsp,1:nsp] <-  inverse(sigma.beta.sp * beta.mat[,])
    ## draw species-specific slopes
    psi.beta.sp[1:nsp] ~ dmnorm(mu.beta.sp[], tau.beta.sp.mat[,])
    
    for(sp in 1:nsp) {
      for(site in 1:nsite) {
        for(yr in 1:nyr) {        
          
          logit(psi[site,yr,sp]) <-
            psi.0[sp] + 
              psi.beta.sp[sp]*env[site]
          Z[site,yr,sp] ~ dbern(psi[site,yr,sp])
        }
      }
    }
  }', fill = TRUE)
  sink()
  
  attach(d$data)
  res <- jags.parallel(data=list('Z', 'nsp', 'nsite', 'nyr', 
                         'env', 'VCOV', 'ID'),
                       parameters.to.save=d$params,
                       model.file='model.jags',
                       n.thin=10,
                       n.iter=1100,
                       n.burnin=100,
                       n.chains=3,
                       working.directory=NULL)
  detach(d$data)
  res
}
