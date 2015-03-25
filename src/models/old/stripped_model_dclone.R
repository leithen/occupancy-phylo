## specify the parameters to be monitored
get.params <- function()
  c('lambda',
    'sigma.beta.sp',
    'mu.beta')


## specify how to update priors over multiple iterations of the model
update.priors <- function(x) {
  if(missing(x)) {
    np <- length(get.params())
    return(cbind(c(0,0,0.01), c(1, 0.01, 10)))  #here are the first, and second distribution parameters for lambda, mu.beta, and simga.beta.sp respectivly. Set here to be for the dunif in the case of lambda and sigma.beta.sp
  } else {
    ncl <- nclones(x)
    if(is.null(ncl)) ncl <- 1

    first.param <- coef(x)
    second.param <- dcsd(x)
    
    ## ## Dealing with lambda if it's normal through a link	
    ## logit.lambda <- mcmcapply(x[,'lambda'], logit)
    ## estimates['lambda'] <- mean(logit.lambda)
    ## se['lambda'] <- sd(logit.lambda) * sqrt(ncl)

    ## Dealing with lambda if it's a dunif
    first.param ['lambda'] <- mcmcapply(x[,'lambda'], min)
    second.param['lambda'] <- mcmcapply(x[,'lambda'], max)

    ## ## Dealing with sigma.beta.sp if it's normal through a link
    ##	log.sigma <- mcmcapply(x[,'sigma.beta.sp'], log)
    ##	estimates['sigma.beta.sp'] <- mean(log.sigma)
    ##	se['sigma.beta.sp'] <- sd(log.sigma) * sqrt(ncl)

    ## Dealing with sigma.beta.sp if it's a dunif
    first.param ['sigma.beta.sp'] <- mcmcapply(x[,'sigma.beta.sp'], min)
    second.param['sigma.beta.sp'] <- mcmcapply(x[,'sigma.beta.sp'], max)
    
    return(cbind(first.param, second.param))
  }
}


## Using logit lambda to update priors
## Specify how to update priors over multiple iterations of the model

update.priors.logit.lambda <- function(x) {
  if(missing(x)) {
    np <- length(get.params())
    return(cbind(c(0,0,0.01), c(1, 0.01, 10)))  #here are the first, and second distribution parameters for lambda, mu.beta, and simga.beta.sp respectivly. Set here to be for the dunif in the case of lambda and sigma.beta.sp
  } else {
    ncl <- nclones(x)
    if(is.null(ncl)) ncl <- 1

    first.param <- coef(x)
    second.param <- 1/dcsd(x)
    
    ## Dealing with lambda if it's normal through a link	
    logit.lambda <- mcmcapply(x[,'lambda'], logit)
    first.param['lambda'] <- mean(logit.lambda)
    second.param['lambda'] <- 1/(sd(logit.lambda) * sqrt(ncl))

    ## ## Dealing with lambda if it's a dunif
    ## first.param ['lambda'] <- mcmcapply(x[,'lambda'], min)
    ## second.param['lambda'] <- mcmcapply(x[,'lambda'], max)

    ## ## Dealing with sigma.beta.sp if it's normal through a link
    ##	log.sigma <- mcmcapply(x[,'sigma.beta.sp'], log)
    ##	estimates['sigma.beta.sp'] <- mean(log.sigma)
    ##	se['sigma.beta.sp'] <- sd(log.sigma) * sqrt(ncl)

    ## Dealing with sigma.beta.sp if it's a dunif
    first.param ['sigma.beta.sp'] <- mcmcapply(x[,'sigma.beta.sp'], min)
    second.param['sigma.beta.sp'] <- mcmcapply(x[,'sigma.beta.sp'], max)
    
    return(cbind(first.param, second.param))
  }
}


run.dclones.model <- function() {
  
  ## species-specific intercepts and the effect of the environment are
  ## modeled as fixed effects
  mu.beta ~ dnorm(0, 0.001 )
  for(sp in 1:nsp) {
    psi.0[sp] ~ dnorm(0,0.001)
    mu.beta.sp[sp] <- mu.beta
  }
  
  ## incorporate phylogenetic covariance structure into
  ## species-specific slope expectations.  lambda allows us to scale
  ## the phylogenetic correlation matrix
  ##specifying logit lambda prioir for ease of updating priors

  ## lambda.prior ~ dnorm(pr[1,1], (pr[1,2]) )
  ## logit(lambda) <- lambda.prior

  lambda ~ dunif(0,1)
  
  beta.mat[1:nsp,1:nsp] <- lambda*VCOV[,] + (1-lambda)*ID[,]

  ## convert to precision (or co-precision) matrix
  ## Here using log sigma prior for ease of updating
  ## log.sigma.beta.sp ~ dnorm(pr[3,1], (pr[3,2]) )
  ## sigma.beta.sp <- exp(log.sigma.beta.sp)

  sigma.beta.sp ~ dunif(0.01,10)
  
  tau.beta.sp.mat[1:nsp,1:nsp] <- inverse(sigma.beta.sp * beta.mat[,])
  ## draw species-specific slopes
  psi.beta.sp[1:nsp] ~ dmnorm(mu.beta.sp[], tau.beta.sp.mat[,])

  for(clone in 1:k) {  
    for(sp in 1:nsp) {
      for(site in 1:nsite) {
        for(yr in 1:nyr) {        
          
          logit(psi[site,yr,sp,clone])  <- 
            psi.0[sp] + 
              psi.beta.sp[sp]*env[site]
          Z[site,yr,sp,clone] ~ dbern(psi[site,yr,sp,clone])
        }
      }
    }
  }
}

## Trying to make dclone work better. Apparently cloning ALL data is
## somehow important. 
run.dclones.model.II <- function() {
  
  ## species-specific intercepts and the effect of the environment are
  ## modeled as fixed effects
  mu.beta ~ dnorm(0, 0.001)
  for(sp in 1:nsp) {
    psi.0[sp] ~ dnorm(0,0.001)

    for(clone in 1:k) {
      mu.beta.sp[sp,clone] <- mu.beta
    }
  }
  
  ## incorporate phylogenetic covariance structure into
  ## species-specific slope expectations.  lambda allows us to scale
  ## the phylogenetic correlation matrix specifying logit lambda
  ## prioir for ease of updating priors

  ##  lambda.prior ~ dnorm(pr[1,1], (pr[1,2]) )
  ##  logit(lambda) <- lambda.prior
  lambda ~ dunif(0,1)

  ## convert to precision (or co-precision) matrix
  ## Here using log sigma prior for ease of updating
  ##  	log.sigma.beta.sp ~ dnorm(pr[3,1], (pr[3,2]) )
  ##	sigma.beta.sp <- exp(log.sigma.beta.sp)

  sigma.beta.sp ~ dunif(0.01,10)

  for(clone in 1:k) {  
    beta.mat[1:nsp,1:nsp,clone] <-
      lambda*VCOV[,,clone] + (1-lambda)*ID[,,clone]	 
    tau.beta.sp.mat[1:nsp,1:nsp,clone] <-
      inverse(sigma.beta.sp * beta.mat[,,clone])
    ## draw species-specific slopes
    
    psi.beta.sp[1:nsp,clone] ~ dmnorm(mu.beta.sp[,clone],
                                       tau.beta.sp.mat[,,clone])
    ## In truth I think this line might be the only one that needs to
    ## becloned. or rather psi.beta.sp needs to becloned, but not
    ## mu.beta.sp or tau.beta.sp.mat. It somehow makes sense to me
    ## that everything coming from a distribution that is not a prior
    ## would need to becloned.

    for(sp in 1:nsp) {
      for(site in 1:nsite) {
        for(yr in 1:nyr) {        
          
          logit(psi[site,yr,sp,clone])  <- 
            psi.0[sp] + 
              psi.beta.sp[sp,clone]*env[site]
          Z[site,yr,sp,clone] ~ dbern(psi[site,yr,sp,clone])
        }
      }
    }
  }
}


## Final test to see whether env also needs to becloned. I can't
## imagine this matters at all, but since the other ones did for no
## easily explicable reason it might be the case that cloning env is
## crucial too 

run.dclones.model.III <- function() {
  
  ## species-specific intercepts and the effect of the environment are
  ## modeled as fixed effects
  mu.beta ~ dnorm(0, 0.001 )
  for(sp in 1:nsp) {
    psi.0[sp] ~ dnorm(0,0.001)

    for(clone in 1:k) {
      mu.beta.sp[sp,clone] <- mu.beta
    }
  }

  
  ## incorporate phylogenetic covariance structure into
  ## species-specific slope expectations.  lambda allows us to scale
  ## the phylogenetic correlation matrix
  ##specifying logit lambda prioir for ease of updating priors

  ##  lambda.prior ~ dnorm(pr[1,1], (pr[1,2]) )
  ##  logit(lambda) <- lambda.prior

  lambda ~ dunif(0,1)

  ## convert to precision (or co-precision) matrix
  ## Here using log sigma prior for ease of updating
  ##  	log.sigma.beta.sp ~ dnorm(pr[3,1], (pr[3,2]) )
  ##	sigma.beta.sp <- exp(log.sigma.beta.sp)

  sigma.beta.sp ~ dunif(0.01,10)

  for(clone in 1:k) {  
    beta.mat[1:nsp,1:nsp,clone] <- lambda*VCOV[,,clone] + (1-lambda)*ID[,,clone]	
    tau.beta.sp.mat[1:nsp,1:nsp,clone] <- inverse(sigma.beta.sp * beta.mat[,,clone])
    ## draw species-specific slopes
    psi.beta.sp[1:nsp,clone] ~ dmnorm(mu.beta.sp[,clone], tau.beta.sp.mat[,,clone])

    for(sp in 1:nsp) {
      for(site in 1:nsite) {
        for(yr in 1:nyr) {        
          
          logit(psi[site,yr,sp,clone])  <- 
            psi.0[sp] + 
              psi.beta.sp[sp,clone]*env[site,clone]
          Z[site,yr,sp,clone] ~ dbern(psi[site,yr,sp,clone])
        }
      }
    }
  }
}



############################################################
#### Update priors to speed convergence in successive iterations
############################################################



###Trying to make dclone work better. Apparently cloning ALL data is somehow important.
run.dclones.model.II.update <- function() {
  
  ## species-specific intercepts and the effect of the environment are
  ## modeled as fixed effects


  ## LONG CONVERGENCE TIMES NEEDED!
  mu.beta ~ dnorm(pr[2,1], pr[2,2] )
  for(sp in 1:nsp) {
    psi.0[sp] ~ dnorm(0,0.001)

    for(clone in 1:k) {
      mu.beta.sp[sp,clone] <- mu.beta
    }
  }

  
  ## incorporate phylogenetic covariance structure into
  ## species-specific slope expectations.  lambda allows us to scale
  ## the phylogenetic correlation matrix
  ##specifying logit lambda prioir for ease of updating priors

  ##  lambda.prior ~ dnorm(pr[1,1], (pr[1,2]) )
  ##  logit(lambda) <- lambda.prior

  lambda ~ dunif(pr[1,1], (pr[1,2]) )

  


  ## convert to precision (or co-precision) matrix
  ## Here using log sigma prior for ease of updating
  ##  	log.sigma.beta.sp ~ dnorm(pr[3,1], (pr[3,2]) )
  ##	sigma.beta.sp <- exp(log.sigma.beta.sp)

  sigma.beta.sp ~ dunif(pr[3,1], (pr[3,2]) )
  


  for(clone in 1:k) {  
    beta.mat[1:nsp,1:nsp,clone] <- lambda*VCOV[,,clone] + (1-lambda)*ID[,,clone]	
    tau.beta.sp.mat[1:nsp,1:nsp,clone] <- inverse(sigma.beta.sp * beta.mat[,,clone])
    ## draw species-specific slopes
    
    psi.beta.sp[1:nsp,clone] ~ dmnorm(mu.beta.sp[,clone], tau.beta.sp.mat[,,clone]) ##In truth I think this line might be the only one that needs to becloned. or rather psi.beta.sp needs to becloned, but not mu.beta.sp or tau.beta.sp.mat. It somehow makes sense to me that everything coming from a distribution that is not a prior would need to becloned.

    for(sp in 1:nsp) {
      for(site in 1:nsite) {
        for(yr in 1:nyr) {        
          
          logit(psi[site,yr,sp,clone])  <- 
            psi.0[sp] + 
              psi.beta.sp[sp,clone]*env[site]
          Z[site,yr,sp,clone] ~ dbern(psi[site,yr,sp,clone])
        }
      }
    }
  }
}

############################################################


############################################################
#### Update priors to speed convergence in successive iterations. Using logit(lambda)
############################################################



###Trying to make dclone work better. Apparently cloning ALL data is somehow important.
run.dclones.model.II.update.logit.lambda <- function() {
  
  ## species-specific intercepts and the effect of the environment are
  ## modeled as fixed effects


  ## LONG CONVERGENCE TIMES NEEDED!
  mu.beta ~ dnorm(pr[2,1], pr[2,2] )
  for(sp in 1:nsp) {
    psi.0[sp] ~ dnorm(0,0.001)

    for(clone in 1:k) {
      mu.beta.sp[sp,clone] <- mu.beta
    }
  }

  
  ## incorporate phylogenetic covariance structure into
  ## species-specific slope expectations.  lambda allows us to scale
  ## the phylogenetic correlation matrix
  ##specifying logit lambda prioir for ease of updating priors

  lambda.prior ~ dnorm(pr[1,1], (pr[1,2]) )
  logit(lambda) <- lambda.prior

  ##	lambda ~ dunif(pr[1,1], (pr[1,2]) )

  


  ## convert to precision (or co-precision) matrix
  ## Here using log sigma prior for ease of updating
  ##  	log.sigma.beta.sp ~ dnorm(pr[3,1], (pr[3,2]) )
  ##	sigma.beta.sp <- exp(log.sigma.beta.sp)

  sigma.beta.sp ~ dunif(pr[3,1], (pr[3,2]) )
  


  for(clone in 1:k) {  
    beta.mat[1:nsp,1:nsp,clone] <- lambda*VCOV[,,clone] + (1-lambda)*ID[,,clone]	
    tau.beta.sp.mat[1:nsp,1:nsp,clone] <- inverse(sigma.beta.sp * beta.mat[,,clone])
    ## draw species-specific slopes
    
    psi.beta.sp[1:nsp,clone] ~ dmnorm(mu.beta.sp[,clone], tau.beta.sp.mat[,,clone]) ##In truth I think this line might be the only one that needs to becloned. or rather psi.beta.sp needs to becloned, but not mu.beta.sp or tau.beta.sp.mat. It somehow makes sense to me that everything coming from a distribution that is not a prior would need to becloned.

    for(sp in 1:nsp) {
      for(site in 1:nsite) {
        for(yr in 1:nyr) {        
          
          logit(psi[site,yr,sp,clone])  <- 
            psi.0[sp] + 
              psi.beta.sp[sp,clone]*env[site]
          Z[site,yr,sp,clone] ~ dbern(psi[site,yr,sp,clone])
        }
      }
    }
  }
}

############################################################



## **************************************************
run.R2jags.model <- function(d) {
  
  sink('model.jags')
  cat('model {

    ## species-specific intercepts and the effect of the environment are
    ## modeled as fixed effects
    mu.beta ~ dnorm(0,0.001)
    for(sp in 1:nsp) {
      psi.0[sp] ~ dnorm(0,0.001)
      mu.beta.sp[sp] <- mu.beta
    }
    
    ## incorporate phylogenetic covariance structure into
    ## species-specific slope expectations.  lambda allows us to scale
    ## the phylogenetic correlation matrix
    lambda ~ dunif(0,1)
    beta.mat[1:nsp,1:nsp] <- lambda*VCOV[,] + (1-lambda)*ID[,]
    ## convert to precision (or co-precision) matrix
    sigma.beta.sp ~ dunif(0.01,100)
    tau.beta.sp.mat[1:nsp,1:nsp] <- inverse(sigma.beta.sp * beta.mat[,])
    ## draw species-specific slopes
    psi.beta.sp[1:nsp] ~ dmnorm(mu.beta.sp[], tau.beta.sp.mat[,])
    
    for(clone in 1:k) {
      for(sp in 1:nsp) {
        for(site in 1:nsite) {
          for(yr in 1:nyr) {        
            
            logit(psi[site,yr,sp,clone])  <- 
              psi.0[sp] + 
                psi.beta.sp[sp]*env[site]
            Z[site,yr,sp,clone] ~ dbern(psi[site,yr,sp,clone])
          }
        }
      }
    }
  }', fill = TRUE)
  sink()
  
  attach(d$data)
  res <- jags.parallel(data=list('Z', 'nsp', 'nsite', 'nyr', 
                         'env', 'VCOV', 'ID', 'nk'),
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
## **************************************************
