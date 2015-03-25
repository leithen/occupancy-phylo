bat.occ.mod <- function() {

  ## detectability
  
  mu.p.0.prior~dunif(0,1)
  mu.p.0  <- logit(mu.p.0.prior)
  
  sigma.p.0 ~ dunif(0,20)
  tau.p.0 <- 1/(sigma.p.0*sigma.p.0)

  ## species-specific detectability
  for(sp in 1:nsp) {
    p.0[sp] ~ dnorm(mu.p.0,tau.p.0)
    logit(p[sp]) <- p.0[sp]
  }
  
  ## random effect of site on occupancy
  sigma.psi.site ~ dunif(0,20)
  tau.psi.site <- 1/(sigma.psi.site*sigma.psi.site)
  for(site in 1:nsite) {
    psi.site[site] ~ dnorm(0, tau.psi.site)
  }

  ## construct species-specific occupancy intercepts which include
  ## both a random species effect and flexibility to have each species
  ## vary by year

  ## random effect of year on occupancy
  sigma.psi.year ~ dunif(0,20)
  tau.psi.year <- 1/(sigma.psi.year*sigma.psi.year)

   ## species-specific intercepts are modeled as fixed effects
   ## But each year intercept is allowed to vary around the mean
   
  for(sp in 1:nsp) {
    psi.0.sp.prior[sp]~dunif(0,1)
  	
  	psi.0.sp[sp]<-logit(psi.0.sp.prior[sp])
  	
  	for (yr in 1:nyr) {
  		psi.0.yr.sp[yr, sp]~ dnorm(psi.0.sp[sp], tau.psi.year)	
 	   	}
    }
  
  
  ## intercept for model on slope

  sigma.beta~dunif(0,20)
  tau.beta<-1/(sigma.beta* sigma.beta)
  
  ##Effect of diets on slope (i.e. habitat affiliation)
    trait.beta~dnorm(0, 0.0001) #Effect of trait on psi.beta.sp
  	mu.beta  ~ dnorm(0, 0.0001) #Intercept of psi.beta.sp



  for(sp in 1:nsp) {  	
  	rand.beta.sp[sp] ~ dnorm(0, tau.beta)
  	psi.beta.sp[sp]<- mu.beta + Trait[sp]*trait.beta 
  						+ rand.beta.sp[sp] 
  }
  
  for(sp in 1:nsp) {
    for(site in 1:nsite) {
      for(yr in 1:nyr) {        
        
        logit(psi[site,yr,sp]) <-
          psi.0.yr.sp[yr, sp] + 
            psi.beta.sp[sp]*env[site] +
              psi.site[site]

        Z[site,yr,sp] ~ dbern(psi[site,yr,sp])
        mu.p[site,yr,sp] <- Z[site,yr,sp]*p[sp]
        for(rep in 1:nrep[site, yr]) {
        	
          X[site,yr,rep,sp] ~ dbern(mu.p[site,yr,sp])
        }
      }
    }
  }
}

get.params<-function() {
	c(
	"psi.0.sp",
	"psi.0.yr.sp",
	"mu.p.0",
	"sigma.p.0",
	"p.0",
	"sigma.psi.site",
	"sigma.psi.year",
	"mu.beta",
	"trait.beta",
	"sigma.beta",	
	"psi.beta.sp",
	"psi.site"
	)
}
