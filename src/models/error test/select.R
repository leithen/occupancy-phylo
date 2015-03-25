run.R2jags.model <- function(d,
                  nt=10,
                  ni=1100,
                  nb=100,
                  nc=3) {
  
  sink('model.jags')
  cat('model {


	mu~dnorm(0,0.01)
	
	for (sp in 1:nsp) {
		mu.beta.sp[sp]<-mu
	}
	
	IND.lambda~dbern(0.5)    #Establish indicator for lambda

    lambda.pre ~ dunif(0.0001, 0.9999) #Establish range that lambda can search in
    lambda.holder~dnorm(MU.lambda, TAU.lambda) #Establish a place holder to maintain model mixing for when IND.lambda is in 0 state
    
    logit.lambda <- logit(lambda.pre)*IND.lambda + (1-IND.lambda)*lambda.holder
    logit(lambda) <- logit.lambda
    
    
    beta.mat[1:nsp,1:nsp] <- lambda*IND.lambda*VCOV[,] + (1-lambda*IND.lambda)*ID[,]
    
    sigma.psi.beta ~ dunif(0,20)
    tau.beta.sp.mat[1:nsp,1:nsp] <- inverse(sigma.psi.beta^2*beta.mat[,])
    ## draw species-specific slopes
    
    
    Y[1:nsp] ~ dmnorm(mu.beta.sp[], tau.beta.sp.mat[,])

  }', fill = TRUE)
  sink()
	
	
	  attach(d$data)
  res <- jags.parallel(data=list('Y', 'nsp', 'VCOV', 'ID', 'MU.lambda', 'TAU.lambda'),
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
  c('logit.lambda',
  	'lambda',
  	'IND.lambda',
    'sigma.psi.beta',
    'mu')



	