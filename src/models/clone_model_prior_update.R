run.dclones.model.lam.free <- function() {

  ## priors and direct derivatives  *** not to be cloned ***

  ## detectability
  mu.p.0 ~ dnorm(pr[2,1], pr[2,2])
  sigma.p.0 ~ dunif(0,20)
  tau.p.0 <- 1/(sigma.p.0*sigma.p.0)

  ## random effect of site on occupancy
  sigma.psi.site ~ dunif(0, 20)
  tau.psi.site <- 1/(sigma.psi.site*sigma.psi.site)

  ## random effect of year on occupancy
  sigma.psi.year ~ dunif(0,20)
  tau.psi.year <- 1/(sigma.psi.year*sigma.psi.year)

  ## fixed effect of species-level covariate on mean occupancy
  ## psi.0.trait ~ dnorm(0,0.001)
  psi.0.trait <- 0
  ## species-specific slope expectation
  mu.psi.beta ~ dnorm(pr[3,1], pr[3,2])
  ## effect of species-level covariate on slope
  ## beta.trait ~ dnorm(0,0.001)
  beta.trait <- 0
  
  for(sp in 1:nsp) {
    ## construct species-specific occupancy intercepts which include
    ## both a random species effect and fixed effect based on the
    ## species trait value
    psi.0.sp[sp] ~ dnorm(pr[sp+3,1], pr[sp+3,2])
    psi.0[sp] <- psi.0.sp[sp] + psi.0.trait*trait[sp]
    ## species-specific slope expectation terms
    mu.beta.sp[sp] <- mu.psi.beta + beta.trait*trait[sp]
  }
  
  ## lambda allows us to scale the phylogenetic correlation matrix
  lambda.pre ~ dnorm(pr[1,1], pr[1,2])
  logit(lambda) <- lambda.pre

  ## incorporate phylogenetic covariance structure into
  ## species-specific slope expectations.
  beta.mat[1:nsp,1:nsp] <- lambda*VCOV[,] + (1-lambda)*ID[,]
  ## convert to precision (or co-precision) matrix - While lambda
  ## above affects the relative weighting on the diagonal in precision
  ## matrix, sigma.psi.beta below affects the overall magnitudes of the
  ## entries and thus the spread of the species-specific responses.
  sigma.psi.beta ~ dunif(0.01,20)
  tau.beta.sp.mat[1:nsp,1:nsp] <- inverse((sigma.psi.beta^2)*beta.mat[,])

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


run.dclones.model.ou.free <- function() {

  ## priors and direct derivatives  *** not to be cloned ***

  ## detectability
  mu.p.0 ~ dnorm(pr[2,1], pr[2,2])
  sigma.p.0 ~ dunif(0,20)
  tau.p.0 <- 1/(sigma.p.0*sigma.p.0)

  ## random effect of site on occupancy
  sigma.psi.site ~ dunif(0, 20)
  tau.psi.site <- 1/(sigma.psi.site*sigma.psi.site)

  ## random effect of year on occupancy
  sigma.psi.year ~ dunif(0,20)
  tau.psi.year <- 1/(sigma.psi.year*sigma.psi.year)

  ## fixed effect of species-level covariate on mean occupancy
  ## psi.0.trait ~ dnorm(0,0.001)
  psi.0.trait <- 0
  ## species-specific slope expectation
  mu.psi.beta ~ dnorm(pr[3,1], pr[3,2])
  ## effect of species-level covariate on slope
  ## beta.trait ~ dnorm(0,0.001)
  beta.trait <- 0
  
  for(sp in 1:nsp) {
    ## construct species-specific occupancy intercepts which include
    ## both a random species effect and fixed effect based on the
    ## species trait value
    psi.0.sp[sp] ~ dnorm(pr[sp+3,1], pr[sp+3,2])
    psi.0[sp] <- psi.0.sp[sp] + psi.0.trait*trait[sp]
    ## species-specific slope expectation terms
    mu.beta.sp[sp] <- mu.psi.beta + beta.trait*trait[sp]
  }
  
  ## lambda allows us to scale the phylogenetic correlation matrix

  alpha.pre ~ dnorm(pr[1,1], pr[1,2])
  log(alpha) <- alpha.pre

#	alpha~dunif(1,2)
	
  ## incorporate phylogenetic covariance structure into
  ## species-specific slope expectations.
  sigma.psi.beta ~ dunif(0.01,20)
for (i in 1:nsp) {
	for (j in 1:nsp) {
        beta.mat[i,j] <- ((sigma.psi.beta^2)/(2*alpha))*exp(-2*alpha*(1-VCOV[i,j]))*(1-exp(-2*alpha*VCOV[i,j])) 
		}
	}
  ## convert to precision (or co-precision) matrix - While lambda
  ## above affects the relative weighting on the diagonal in precision
  ## matrix, sigma.psi.beta below affects the overall magnitudes of the
  ## entries and thus the spread of the species-specific responses.
  tau.beta.sp.mat[1:nsp,1:nsp] <- inverse(beta.mat[,])

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



run.dclones.model.lam.0 <- function() {

  ## priors and direct derivatives  *** not to be cloned ***

  ## detectability
  mu.p.0 ~ dnorm(pr[1,1], pr[1,2])
  sigma.p.0 ~ dunif(0,20)
  tau.p.0 <- 1/(sigma.p.0*sigma.p.0)

  ## random effect of site on occupancy
  sigma.psi.site ~ dunif(0, 20)
  tau.psi.site <- 1/(sigma.psi.site*sigma.psi.site)

  ## random effect of year on occupancy
  sigma.psi.year ~ dunif(0,20)
  tau.psi.year <- 1/(sigma.psi.year*sigma.psi.year)

  ## fixed effect of species-level covariate on mean occupancy
  ## psi.0.trait ~ dnorm(0,0.001)
  psi.0.trait <- 0
  ## species-specific slope expectation
  mu.psi.beta ~ dnorm(pr[2,1], pr[2,2])
  ## effect of species-level covariate on slope
  ## beta.trait ~ dnorm(0,0.001)
  beta.trait <- 0
  
  for(sp in 1:nsp) {
    ## construct species-specific occupancy intercepts which include
    ## both a random species effect and fixed effect based on the
    ## species trait value
    psi.0.sp[sp] ~ dnorm(pr[sp+2,1], pr[sp+2,2])
    psi.0[sp] <- psi.0.sp[sp] + psi.0.trait*trait[sp]
    ## species-specific slope expectation terms
    mu.beta.sp[sp] <- mu.psi.beta + beta.trait*trait[sp]
  }
  
  ## lambda allows us to scale the phylogenetic correlation matrix
  lambda <- 0

  ## incorporate phylogenetic covariance structure into
  ## species-specific slope expectations.
  
  beta.mat[1:nsp,1:nsp] <- lambda*VCOV[,] + (1-lambda)*ID[,]
  
  ## convert to precision (or co-precision) matrix - While lambda
  ## above affects the relative weighting on the diagonal in precision
  ## matrix, sigma.psi.beta below affects the overall magnitudes of the
  ## entries and thus the spread of the species-specific responses.
  sigma.psi.beta ~ dunif(0.01,20)
  tau.beta.sp.mat[1:nsp,1:nsp] <- inverse((sigma.psi.beta^2) * beta.mat[,])

  
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

run.dclones.model.no.det <- function() {

  ## priors and direct derivatives  *** not to be cloned ***

  ## random effect of site on occupancy
  sigma.psi.site ~ dunif(0, 20)
  tau.psi.site <- 1/(sigma.psi.site*sigma.psi.site)

  ## random effect of year on occupancy
  sigma.psi.year ~ dunif(0,20)
  tau.psi.year <- 1/(sigma.psi.year*sigma.psi.year)

  ## fixed effect of species-level covariate on mean occupancy
  ## psi.0.trait ~ dnorm(0,0.001)
  psi.0.trait <- 0
  ## species-specific slope expectation
  mu.psi.beta ~ dnorm(0, 0.001)
  ## effect of species-level covariate on slope
  ## beta.trait ~ dnorm(0,0.001)
  beta.trait <- 0
  
  for(sp in 1:nsp) {
    ## construct species-specific occupancy intercepts which include
    ## both a random species effect and fixed effect based on the
    ## species trait value
    psi.0.sp[sp] ~ dnorm(0, 0.01)
    psi.0[sp] <- psi.0.sp[sp] + psi.0.trait*trait[sp]
    ## species-specific slope expectation terms
    mu.beta.sp[sp] <- mu.psi.beta + beta.trait*trait[sp]
  }
  
  ## lambda allows us to scale the phylogenetic correlation matrix
  lambda.pre ~ dnorm(0, 0.01)
  logit(lambda) <- lambda.pre

  ## incorporate phylogenetic covariance structure into
  ## species-specific slope expectations.
  beta.mat[1:nsp,1:nsp] <- lambda*VCOV[,] + (1-lambda)*ID[,]
  ## convert to precision (or co-precision) matrix - While lambda
  ## above affects the relative weighting on the diagonal in precision
  ## matrix, sigma.psi.beta below affects the overall magnitudes of the
  ## entries and thus the spread of the species-specific responses.
  sigma.psi.beta ~ dunif(0.01,20)
  tau.beta.sp.mat[1:nsp,1:nsp] <- inverse((sigma.psi.beta^2)*beta.mat[,])

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
      for(site in 1:nsite) {
        for(yr in 1:nyr) {        
          
          logit(psi[site,yr,sp,clone]) <- 
            psi.0[sp] + 
              psi.beta.sp[sp,clone]*env[site] +
                psi.site[site,clone] +
                  psi.year[yr,clone]

          for(rep in 1:nrep[site,yr,sp]) {
            X[site,yr,rep,sp,clone] ~ dbern(psi[site,yr,sp,clone])
          }
        }
      }
    }
  }
}



## specify the parameters to be monitored
get.params <- function() {
	c('mu.p.0',
    'sigma.p.0',
    'sigma.psi.site',
    'sigma.psi.year',
    'psi.0.sp',
    'mu.psi.beta',
    'sigma.psi.beta',
    'lambda') }
 
get.params.ou <- function() {
	c('mu.p.0',
    'sigma.p.0',
    'sigma.psi.site',
    'sigma.psi.year',
    'psi.0.sp',
    'mu.psi.beta',
    'sigma.psi.beta',
    'alpha') }
 
## specify the parameters to be monitored
get.params.no.det <- function() {
	c('sigma.psi.site',
    'sigma.psi.year',
    'psi.0.sp',
    'mu.psi.beta',
    'sigma.psi.beta',
    'lambda') } 
 
    
## function to update priors

## specify how to update priors over multiple iterations of the model
update.priors.lam.free <- function(x) {
  if(missing(x)) {
    np <- length(get.params()) ###This will be a problem for group parameters, like psi.0.sp
    
    all.params <- c(get.params(), rep("psi.0.sp", nsp-1))
    
    all.param.names <- all.params [order(all.params)]
    all.param.names[all.param.names=='psi.0.sp'] <- paste(all.param.names[all.param.names=='psi.0.sp'], 1:nsp)
    
    
    pr.df <- data.frame(first=rep(0, length(all.param.names)), second=rep(0.001,length(all.param.names)), row.names=all.param.names)
    
    return(pr.df)
    #return(cbind(0,0.01))  #Start with only one prior to update. Say mu.psi.beta
  
  } else {
    ncl <- nclones(x)
    if(is.null(ncl)) ncl <- 1

    first.param <- coef(x) ###This will be a problem for group parameters, like psi.0.sp
    second.param <- 1/(dcsd(x)^2)
    
    ## ## Dealing with lambda if it's normal through a link	
    logit.lambda <- mcmcapply(x[,'lambda'], logit)
    logit.lambda <- logit.lambda[! logit.lambda==Inf | logit.lambda==-Inf]
    
    first.param['lambda'] <- mean(logit.lambda)
    second.param['lambda'] <- 1/(sd(logit.lambda)^2)

    ## Dealing with lambda if it's a dunif
    ## first.param ['lambda'] <- mcmcapply(x[,'lambda'], min)
    ## second.param['lambda'] <- mcmcapply(x[,'lambda'], max)

    ## ## Dealing with sigma.beta.sp if it's normal through a link
    ## log.sigma <- mcmcapply(x[,'sigma.psi.site'], log)
    ## log.sigma <- log.sigma[! log.sigma ==Inf | log.sigma ==-Inf]    
    ## first.param['sigma.psi.site'] <- mean(log.sigma)
    ## second.param['sigma.psi.site'] <- 1/(sd(log.sigma)^2)

    ## Dealing with sigma.beta.sp if it's a dunif
    ## first.param ['sigma.beta.sp'] <- mcmcapply(x[,'sigma.beta.sp'], min)
    ## second.param['sigma.beta.sp'] <- mcmcapply(x[,'sigma.beta.sp'], max)
    
    return(cbind(first.param, second.param))
  }
}


## specify how to update priors over multiple iterations of the model
update.priors.ou.free <- function(x) {
  if(missing(x)) {
    np <- length(get.params()) ###This will be a problem for group parameters, like psi.0.sp
    
    all.params <- c(get.params(), rep("psi.0.sp", nsp-1))
    
    all.param.names <- all.params [order(all.params)]
    all.param.names[all.param.names=='psi.0.sp'] <- paste(all.param.names[all.param.names=='psi.0.sp'], 1:nsp)
    
    
    pr.df <- data.frame(first=rep(0, length(all.param.names)), second=rep(0.001,length(all.param.names)), row.names=all.param.names)
    
    pr.df[1,2]<-0.1 ### Setting a narrower prior on alpha, since
    				## it it gets logged in the inner workings
    				## (i.e. alpha.pre ~ dnorm(pr.df[1,1], pr.df[1,2])
    				##  then log(alpha)<-alpha.pre   )
    
    return(pr.df)
    #return(cbind(0,0.01))  #Start with only one prior to update. Say mu.psi.beta
  
  } else {
    ncl <- nclones(x)
    if(is.null(ncl)) ncl <- 1

    first.param <- coef(x) ###This will be a problem for group parameters, like psi.0.sp
    second.param <- 1/(dcsd(x)^2)
    
    ## ## Dealing with alpha if it's normal through a link	
    log.alpha <- mcmcapply(x[,'alpha'], log)
    log.alpha <- log.alpha[! log.alpha==Inf | log.alpha==-Inf]
    
    first.param['alpha'] <- mean(log.alpha)
    second.param['alpha'] <- 1/(sd(log.alpha)^2)

    ## Dealing with lambda if it's a dunif
    ## first.param ['lambda'] <- mcmcapply(x[,'lambda'], min)
    ## second.param['lambda'] <- mcmcapply(x[,'lambda'], max)

    ## ## Dealing with sigma.beta.sp if it's normal through a link
    ## log.sigma <- mcmcapply(x[,'sigma.psi.site'], log)
    ## log.sigma <- log.sigma[! log.sigma ==Inf | log.sigma ==-Inf]    
    ## first.param['sigma.psi.site'] <- mean(log.sigma)
    ## second.param['sigma.psi.site'] <- 1/(sd(log.sigma)^2)

    ## Dealing with sigma.beta.sp if it's a dunif
    ## first.param ['sigma.beta.sp'] <- mcmcapply(x[,'sigma.beta.sp'], min)
    ## second.param['sigma.beta.sp'] <- mcmcapply(x[,'sigma.beta.sp'], max)
    
    return(cbind(first.param, second.param))
  }
}

## specify how to update priors over multiple iterations of the model
update.priors.lam.0 <- function(x) {
  if(missing(x)) {
    rel.params <- get.params()[!get.params()=='lambda']
    all.params <- c(rel.params, rep("psi.0.sp", nsp-1))
    
    all.param.names <- all.params [order(all.params)]
    all.param.names[all.param.names=='psi.0.sp'] <- paste(all.param.names[all.param.names=='psi.0.sp'], 1:nsp)
    
    
    pr.df <- data.frame(first=rep(0, length(all.param.names)), second=rep(0.001,length(all.param.names)), row.names=all.param.names)
    
    return(pr.df)
    #return(cbind(0,0.01))  #Start with only one prior to update. Say mu.psi.beta
  
  } else {
    ncl <- nclones(x)
    if(is.null(ncl)) ncl <- 1

    first.param <- coef(x) ###This will be a problem for group parameters, like psi.0.sp
    second.param <- 1/(dcsd(x)^2)
    
    
    return(cbind(first.param, second.param))
  }
}


## specify how to update priors for no detection version multiple iterations of the model
update.priors.no.det <- function(x) {
  if(missing(x)) {
    rel.params <- get.params.no.det()
    all.params <- c(rel.params, rep("psi.0.sp", nsp-1))
    
    all.param.names <- all.params [order(all.params)]
    all.param.names[all.param.names=='psi.0.sp'] <- paste(all.param.names[all.param.names=='psi.0.sp'], 1:nsp)
    
    
    pr.df <- data.frame(first=rep(0, length(all.param.names)), second=rep(0.001,length(all.param.names)), row.names=all.param.names)
    
    return(pr.df)
    #return(cbind(0,0.01))  #Start with only one prior to update. Say mu.psi.beta
  
  } else {
    ncl <- nclones(x)
    if(is.null(ncl)) ncl <- 1

    first.param <- coef(x) ###This will be a problem for group parameters, like psi.0.sp
    second.param <- 1/(dcsd(x)^2)
    
    
    return(cbind(first.param, second.param))
  }
}

