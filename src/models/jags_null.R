run.R2jags.model <- function(d,
                  nt=10,
                  ni=1100,
                  nb=100,
                  nc=3) {
  
  sink('model.jags')
  cat('model {

    ## detectability
    mu.p.0.prior ~ dunif(0,1)
    mu.p.0 <- logit(mu.p.0.prior)
    ## mu.p.0  ~ dnorm(0,0.001)
    sigma.p.0 ~ dunif(0,20)
    tau.p.0 <- 1/(sigma.p.0*sigma.p.0)

    ## species-specific detectability
    for(sp in 1:nsp) {
      p.0[sp] ~ dnorm(mu.p.0,tau.p.0)
      logit(p[sp]) <- p.0[sp]
    }
    
    ## random effect of site on occupancy
    ## mu.psi.site ~ dnorm(0,0.001)
    sigma.psi.site ~ dunif(0,20)
    tau.psi.site <- 1/(sigma.psi.site*sigma.psi.site)
    for (site in 1:nsite) {
      psi.site[site] ~ dnorm(0, tau.psi.site)
    }

    ## random effect of year on occupancy
    ## mu.psi.year ~ dnorm(0,0.001)
    sigma.psi.year ~ dunif(0,20)
    tau.psi.year <- 1/(sigma.psi.year*sigma.psi.year)
    for (yr in 1:nyr) {
      psi.year[yr] ~ dnorm(0, tau.psi.year)
    }

    ## construct species-specific occupancy intercepts which include
    ## both a random species effect and fixed effect based on the
    ## species trait value
    ##
    ## fixed effect of species-level covariate on mean occupancy
    psi.0.trait ~ dnorm(0,0.001)
    
    ## ## species-specific intercepts are modeled as random effects
    ## mu.psi.0 ~ dnorm(0,0.001)
    ## sigma.psi.0 ~ dunif(0,20)
    ## tau.psi.0 <- 1/(sigma.psi.0*sigma.psi.0)
    ## for(sp in 1:nsp) {
    ##   psi.0.mean[sp] ~ dnorm(mu.psi.0,tau.psi.0)
    ##   psi.0[sp] <- psi.0.mean[sp] + psi.0.trait*trait[sp]
    ## }

    ## species-specific intercepts are modeled as fixed effects
    for(sp in 1:nsp) {
      psi.0.sp.prior[sp]~dunif(0,1)
      psi.0.sp[sp] <- logit(psi.0.sp.prior[sp])
      psi.0[sp] <- psi.0.sp[sp] ## + psi.0.trait*trait[sp] ## No trait for now!
    }
    
    ## create species-specific slopes with regard to the environment
    ## (phylogeny is incorporated here)
    ##
    ## species-specific slope expectation and precision terms
    ## (beta.trait is the effect of species-level covariate on slope)
    mu.psi.beta ~ dnorm(0,0.001)
    beta.trait ~ dnorm(0,0.001)
    for(sp in 1:nsp) {
      mu.beta.sp[sp] <- mu.psi.beta ## + beta.trait*trait[sp] ##No traits for now!
    }
    ## incorporate phylogenetic covariance structure into
    ## species-specific slope expectations.  lambda allows us to scale
    ## the phylogenetic correlation matrix
    lambda.pre <- 0
    
    #logit.lambda <- logit(lambda.pre)
    #logit(lambda) <- logit.lambda
    
    lambda<-0
    
    
    beta.mat[1:nsp,1:nsp] <- lambda*VCOV[,] + (1-lambda)*ID[,]
    ## convert to precision (or co-precision) matrix - While lambda
    ## above affects the relative weighting on the diagonal in
    ## precision matrix, sigma.psi.beta below affects the overall
    ## magnitudes of the entries and thus the spread of the
    ## species-specific responses.
    sigma.psi.beta ~ dunif(0,20)
    tau.beta.sp.mat[1:nsp,1:nsp] <-inverse(sigma.psi.beta^2*beta.mat[,])
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
                       n.thin=nt,
                       n.iter=ni,
                       n.burnin=nb,
                       n.chains=nc,
                       working.directory=NULL)
  detach(d$data)
  
  res
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
	'logit.lambda')