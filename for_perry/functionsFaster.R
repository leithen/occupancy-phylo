expit <- function(x) 1/(1+exp(-x))

## Returns a set of funcions with mu.vec and vcov.mat baked in.
## Exception: the dMuSig function takes a mu vector and a sig2 that
## scales the vcov.mat.
make.mvnorm.funs <- function(mu.vec, vcov.mat) {
  chol.vcov <- chol(vcov.mat)
  inv.chol <- solve(chol.vcov)
  p <- length(mu.vec)
  force(mu.vec)
  logDeterminant <- sum(2*log(diag(chol.vcov)))
  logNormConst <- -0.5*p*log( 2*pi ) - 0.5*logDeterminant
  rd <- function(num.draws) {
    z <- matrix(rnorm(p * num.draws), nrow=p)
    x <- mu.vec + t(chol.vcov) %*% z
    lpx <- logNormConst - 0.5 * colSums(z*z)
    list(x=t(x), lpx=lpx)
  }
  r <- function(num.draws) {
    z <- matrix(rnorm(p * num.draws), nrow=p)
    x <- mu.vec + t(chol.vcov) %*% z
    t(x)
  }
  d <- function(x) {
    v <- if(is.matrix(x)) {
      t(t(x)-mu.vec) %*% inv.chol
    } else {
      (x-mu.vec) %*% inv.chol
    }
    ans <- logNormConst - 0.5 * rowSums(v*v)
    ans
  }
  dMuSig <- function(x, mu, sig2) { ## sig2 is a scalar multipler for cov matrix
    
    v <- if(is.matrix(x)) {
      t(t(x)-mu) %*% inv.chol
    } else {
      (x-mu) %*% inv.chol
    }
    ans <- -0.5*p*log(sig2) + logNormConst - 0.5 * rowSums(v*v)/sig2
    ans
  }
  list(r=r, d=d, rd=rd, dMuSig=dMuSig)
}

## Extracts information needed and generates the mvnorm.funs for the
## specific mu.vec and vcov.mat and the calc.par.LL.faster function
## with all the one-time computations done
setup.jags.LL.faster <- function(dclone.mod, draw.num) {

  model.case <- attr(dclone.mod, 'model.case')
  jags.mod <- attr(dclone.mod, 'jags.mod')
  
  ## set prms
  prms <- jags.mod$data ## this has data and parameters
  
  ## extract raw chains
  chains <- jags.mod$bugs$BUGSoutput$sims.matrix
  ## drop deviance
  chains <- chains[,-which(colnames(chains)=='deviance')]
  
  ## mean vector
  mu.vec <- colMeans(chains)
  ## variance-covariance matrix
  vcov.mat <- cov(chains)	
  
  mvnorm.funs <- make.mvnorm.funs(mu.vec, vcov.mat)

  ## some testing lines
  ## rd <- mvnorm.funs$rd(draw.num)
  ## plot(cov(rd2$x), vcov.mat)
  ## mvnorm.funs$d(rd2$x)[1:10]
  ## dmvnorm(rd2$x,mu.vec, vcov.mat, log=TRUE)[1:10]
  
  par.calc.LL.faster <- 
    make.par.calc.LL.faster(prms=prms,
                            mu.vec=mu.vec,
                            vcov.mat=vcov.mat,
                            model.case=model.case)

  ## list(random.draws=random.draws,
  ##      par.calc.LL.faster=par.calc.LL.faster)
  list(mvnorm.funs=mvnorm.funs,
       par.calc.LL.faster=par.calc.LL.faster)

}

## This is modified to take a setup object, generate the random draws
## here and calc par.calc.LLfaster on the entire block of random.draws
jags.LL.faster <- function(setup,
                           dclone.mod,
                           draw.num,
                           detection=TRUE,
                           return.model.object=TRUE,
                           save.all.LLs = save.all.LLs) {  
  
  ## model.case <- attr(dclone.mod, 'model.case')
  jags.mod <- attr(dclone.mod, 'jags.mod')
  
  ## extract raw chains and drop deviance
  chains <- jags.mod$bugs$BUGSoutput$sims.matrix
  chains <- chains[,-which(colnames(chains)=='deviance')]
  
  mu.vec <- colMeans(chains) # mean vector
  vcov.mat <- cov(chains) # variance-covariance matrix 
  
  rd <- setup$mvnorm.funs$rd(draw.num)
  random.draws <- rd$x
  random.draw.logprobs <- rd$lpx

  if(detection==TRUE)
    vals <- setup$par.calc.LL.faster(1, random.draws=random.draws)
  if(detection==FALSE) {
    browser() ## NEEDS FILLING IN
  }

  ## subtract out the importance weight
  vals <- vals - rd$lpx 

  if(return.model.object==TRUE) {
    attr(dclone.mod, 'lik.vals') <- vals

    mean.vals <- log(mean(exp(vals - mean(vals)))) + mean(vals)
    se.vals <- sqrt(var(vals-mean(vals)))/sqrt(length(vals))
    AIC <- 2*length(coef(dclone.mod)) - 2* mean.vals
    out <- c(mean.vals, se.vals, AIC, draw.num)
    names(out) <- c('LL', 'SE Log Lik.', 'AIC', 'Sample Number')
    out <- round(out, 5)
    attr(dclone.mod, 'stat.sum') <- out

    return(dclone.mod)
  } 
  
  if (return.model.object==FALSE) {
    if (save.all.LLs) {
      return(vals)
    }
    else{
      mean.vals <- log(mean(exp(vals - mean(vals)))) + mean(vals)
      se.vals <- sqrt(var(vals-mean(vals)))/sqrt(length(vals))
      return(c(mean.vals, se.vals))
    }
  }
}

## Does a one-time setup calculations and then returns a function that
## computes from a block of random.draws.
make.par.calc.LL.faster <- function(prms=prms,
                                    mu.vec=mu.vec,
                                    vcov.mat=vcov.mat,
                                    model.case=model.case) {

  ## get important data into easily accessible form
  env <- prms$env
  nsite <- prms $nsite
  nyr <- prms $nyr
  nrep <- max(prms$nrep) ## note: assuming same num of reps
  nsp <- prms $nsp
  X <- prms$X
  
  if(model.case=='phylo.lam' | model.case=='phylo.null') {
    lambda.mat <- prms$lambda*prms$VCOV + (1-prms$lambda)*prms$ID
    lambda <- prms$lambda
  }
  
  ## get mle estimates into easily accessible form
  mu.p.0 <- prms$mu.p.0
  sigma.p.0 <- prms$sigma.p.0
  mu.psi.beta <- prms$mu.psi.beta
  sigma.psi.beta <- prms$sigma.psi.beta
  sigma.psi.site <- prms$sigma.psi.site
  sigma.psi.year <- prms$sigma.psi.year
  psi.0.sp <- prms$psi.0.sp
  psi.0.vec <- psi.0.sp
  
  ## if(model.case=='phylo.ou') {  
  ##   VCOV <- prms$VCOV
  ##   alpha <- prms$alpha
  ## }

  p.0.terms      <- grep('p.0',         names(mu.vec))
  psi.beta.terms <- grep('psi.beta.sp', names(mu.vec))
  psi.site.terms <- grep('psi.site',    names(mu.vec))
  psi.year.terms <- grep('psi.year',    names(mu.vec))

  boolXzero <- as.integer(X==0)
  
  absent.all.reps <- apply(X, c(1,2,4), sum)==0
  Z <- X*0
  for(i in 1:(dim(Z)[3])) {Z[,,i,][absent.all.reps] <- 1}
  boolZzero <- as.vector(Z)
  
  ## Xred<-aperm(X, c(1,2,4,3))
  ## Z.min<-apply(Xred, c(1,2,3), sum)==0
  ## Zred<-array( rep(as.vector(Z.min), nrep), dim=c(nsite, nyr, nsp, nrep))
  ## boolZzero <-aperm(Zred, c(1,2,4,3))
  ## dimnames(boolZzero)<-dimnames(X)
  ## boolZzero<-as.vector(boolZzero)

  
  mvnorm.funs.RR.psi.beta <- make.mvnorm.funs(rep(0, nsp), lambda.mat)
  
  par.calc.LL.faster <- function(chain,
                                 random.draws=random.draws,
                                 use.compiled=F) {

    p.0.mat      <- random.draws[, p.0.terms]
    psi.beta.mat <- random.draws[, psi.beta.terms]
    psi.site.mat <- random.draws[, psi.site.terms]
    psi.year.mat <- random.draws[, psi.year.terms]
    
    if(model.case=='phylo.lam' | model.case=='phylo.null')
      RR.psi.beta <- 
        mvnorm.funs.RR.psi.beta$dMuSig(psi.beta.mat,
                                       rep(mu.psi.beta, nsp),
                                       sigma.psi.beta^2)

    ## if(model.case=='phylo.ou')
    ##   RR.psi.beta <- dmvnorm(psi.beta.vec,
    ##                          mean=rep(mu.psi.beta, nsp),
    ##                          sigma=((sigma.psi.beta^2)/(2*alpha))*exp(-2*alpha*(1-VCOV))*(1-exp(-2*alpha*VCOV)),
    ##                          log=TRUE)
    
    LL.psi.beta.given.theta <- RR.psi.beta

    RR.site <- dnorm(psi.site.mat,
                     mean=0,
                     sd=sigma.psi.site,
                     log=TRUE)
    LL.psi.site.given.theta <- rowSums(RR.site)

    RR.year <- dnorm(psi.year.mat,
                     mean=0,
                     sd=sigma.psi.year,
                     log=TRUE)
    LL.psi.year.given.theta <- rowSums(RR.year) 

    RR.p.0 <- dnorm(p.0.mat,
                    mean=mu.p.0,
                    sd=sigma.p.0,
                    log=TRUE)
    LL.p.given.theta <- rowSums(RR.p.0)
    
    ## Calculates probability of obtaining a certain occupancy status
    ## at a site, at a year, (and if this were actually on
    ## observations rather than occupancy, during a rep), given the
    ## expected psi for the species at the site.
    if(!use.compiled) { 
      LL.X.given.psi.and.p <- numeric(nrow(random.draws))
      for(iRow in 1:nrow(random.draws)) {

        psi.sp.int.vec <-
          rep(psi.site.mat[iRow,],
              ncol(psi.year.mat)*nrep*length(psi.0.vec)) +
                rep(rep(psi.year.mat[iRow,],
                        each=ncol(psi.site.mat)), nrep*length(psi.0.vec)) +
                          rep(psi.0.vec,
                              each=ncol(psi.site.mat)*ncol(psi.year.mat)*nrep) 
        
        ## vector of species-specific environmental slopes
        psi.sp.slope.vec <-
          rep(env, nyr * nrep * ncol(psi.beta.mat)) * ## first part static
            rep(psi.beta.mat[iRow,], each=length(env) * nyr * nrep)
        
        ## occupancy matrix
        psi.mat <- expit(psi.sp.int.vec + psi.sp.slope.vec)
        
        ## Calculate the probability of obtianing the the detection
        ## history given the expected occupancy at a site, and the
        ## estimated detectability of each spp
        ##
        ## As long as species detectability does not depend on site or
        ## year this should work fine. Spp as final layer.
        p.mat <- rep(expit(p.0.mat[iRow,]), each=nsite*nyr*nrep)

        ## compute probabilities of detection histories, given
        ## presence
        X.p.mat <- X
        X.p.mat[X==1] <- p.mat[X==1]
        X.p.mat[X==0] <- 1-p.mat[X==0]
        Z.prod.prob.det <- apply(X.p.mat, c(1,2,4), prod)
        X.psi.mat <- array(psi.mat, dim=dim(X))
        Z.psi.mat <- X.psi.mat[,,1,]
        Z.psi.mat * Z.prod.prob.det + (1-Z.psi.mat)
        QQ <- log(Z.psi.mat)
        ## browser()
        
        ## ## calculate probability of X given psi.mat and p.mat
        ## p.mat[boolXzero==1] <- 1-p.mat[boolXzero==1]
        ## ## QQ <- log(psi.mat*p.mat + (1-psi.mat)*boolXzero)
        ## QQ <- log(psi.mat*p.mat + (1-psi.mat)*boolZzero)
        
        LL.X.given.psi.and.p[iRow] <- sum(QQ)
      }
    } else {
      ## copy of the C function definition for reference:
      ##
      ## void sample_phylo_occupancy_ll(double *expit_p_0_samples,
      ##                                double *psi_beta_samples,
      ##                                double *psi_site_samples, 
      ##                                double *psi_year_samples,
      ##                                double *env,
      ##                                double *psi_0_vec,
      ##                                int *boolXzero,
      ##                                int *pnumSamples,
      ##                                int *pnumSites,
      ##                                int *pnumYears,
      ##                                int *pnumReps,
      ##                                int *pnumSpecies, 
      ##                                double *LL )

      expit.p.0.mat <- expit(p.0.mat)

      LL.X.given.psi.and.p <- .C('sample_phylo_occupancy_ll',
                                 expit.p.0.mat,
                                 psi.beta.mat,
                                 psi.site.mat,
                                 psi.year.mat,
                                 env,
                                 psi.0.vec,
                                 as.integer(boolXzero),
                                 as.integer(boolZzero),
                                 nrow(random.draws),
                                 dim(X)[1],
                                 dim(X)[2],
                                 dim(X)[3],
                                 dim(X)[4],
                                 numeric(nrow(random.draws)),
                                 DUP=FALSE)[[14]]
    }
    
    ## Calculate probability of random effect parameters given the
    ## multi-variate normal we used to draw them from
    ## used
    
    ## calculate the total log likelihood.
    LL.psi.beta.given.theta +
      LL.psi.site.given.theta +
        LL.psi.year.given.theta +
          LL.p.given.theta + 
            LL.X.given.psi.and.p
    
  }
  par.calc.LL.faster
}


anova.pom <- function(object, ...) {
  ancall <- sys.call()
  ancall$verbose <- ancall$test <- NULL
  
  aux <- list(object, ...)
  deparse(substitute(aux[[1]]))
  rt <- length(aux)
  
  statz <- simplify2array(lapply(aux, attr, 'stat.sum'))
  dfModel <- unlist(lapply(lapply(aux, coef), length))
  logLik <- statz[1,]
  AIC <- statz[2,]
  
  aod <- data.frame(df=dfModel, 
                    AIC=AIC, logLik=logLik)

  ddf <- diff(dfModel)
  if(sum(abs(ddf)) > 0) {
    effects <- rep('', rt)
    for (i in 2:rt) {
      if(ddf[i - 1] !=0) {
        effects[i] <- paste(i - 1, i, sep=' vs ')
      }
    }

    pval <- rep(NA, rt - 1)
    ldf <- as.logical(ddf)
    lratio <- -2 * diff(logLik)
    lratio[!ldf] <- NA
    pval[ldf] <- 1 - pchisq(lratio[ldf], df=abs(ddf[ldf]))
    dAIC <- AIC-min(AIC)
    
    
    aod <- data.frame(aod,
                      dAIC=dAIC,
                      Test=effects,
                      L.Ratio=c(NA, lratio),
                      `p-value`=c(NA, pval),
                      check.names=FALSE) 
  }
  
  row.names(aod) <- unlist(lapply(as.list(ancall[-1L]), deparse))
  
  aod
}

imp.samp.plot <- function(object ,...) {
  aux <- list(object, ...)
  ancall <- sys.call()
  ancall$verbose <- ancall$test <- NULL
  
  x.most <-
    max(unlist(lapply(1:length(aux),
                      function(x) length(attr(aux[[x]], 'lik.vals')))))
  y.range <-
    range(unlist(lapply(1:length(aux),
                        function(x) (attr(aux[[x]], 'lik.vals')))))
  colz <- rainbow(length(aux))

  plot(0,0, ylab='Log Likelihoods', xlab='Draw Number', 
       ylim=y.range, xlim=c(1,x.most), type='n')

  ot <- list()
  
  for(i in 1:length(aux)) {
    vals <- attr(aux[[i]], 'lik.vals')
    ot[[i]] <- numeric()
    notch <- length(vals)/1000
    
    for (j in 1:1000) {
      ot[[i]][j] <- log(mean(exp(vals[1:(j*notch)] -
                                 mean(vals[1:(j*notch)])))) +
                                   mean(vals[1:(j*notch)])
    }
    x <- seq(notch, length(vals) , by=notch)

    points(vals, col=alpha(colz[i], 0.2), pch=1)	
  }

  for (i in 1:length(aux)) {
    points(x, ot[[i]], type='l', col='black', lwd=6)
    points(x, ot[[i]], type='l', col=colz[i], lwd=3)
  }

  namez <- unlist(lapply(as.list(ancall[-1L]), deparse))
  legend('bottomleft', namez, lwd=3, col=colz)
}

######################################################################
##  Luke and Leithen's additions to Perry's functions
######################################################################

batch.wrapper <- function(dclone.mod,
                          batch.size,
                          n.batches,
                          save.all.LLs=save.all.LLs) {
  
  draw.num <- batch.size
  
  if(save.all.LLs) {
    LLs <- matrix(nrow=batch.size, ncol=n.batches)
  } else {
    LLs <- numeric()
  }
  
  for (i in 1:n.batches) {
    
    jags.LL.setup <- setup.jags.LL.faster(dclone.mod, draw.num)

    valz <- jags.LL.faster(setup=jags.LL.setup,
                           dclone.mod=dclone.mod,
                           draw.num=draw.num,
                           detection=TRUE,
                           return.model.object=FALSE,
                           save.all.LLs=save.all.LLs)
    if(save.all.LLs) {
      LLs[,i] <- valz
    } else {
      LLs[i] <- valz[1] 
    }
  }
  LLs
}

pom.LL <- function(dclone.mod,
                   n.cores=1,
                   n.batches,
                   batch.size=100000,
                   save.all.LLs=F) {
  
  ## Paritions all draws among available cores. batch.wrapper loops 
  ## on each core for the total number of batches. This results in 
  ## a absolute total number of draws equal to 
  ## n.cores * n.batches * batch.size
  
  ## save.all.LLs argument is included to save each LL draw, rather 
  ## than just a batch wide average. Turning this argument on 
  ## results in massive objects, so probably only worth having 
  ## save.all.LLs = TRUE for debugging.
  
  vals <- mclapply(1:n.cores, function(x)
                   batch.wrapper(dclone.mod=dclone.mod,
                                 batch.size=batch.size, 
                                 n.batches=n.batches,
                                 save.all.LLs=save.all.LLs),
                   mc.cores=n.cores)
  
  vals <- unlist(vals)
  vals <- as.numeric(vals)		
  
  attr(dclone.mod, 'lik.vals') <- vals

  mean.vals <- log(mean(exp(vals - mean(vals)))) + mean(vals)
  AIC <- round(2*length(coef(dclone.mod)) - 2* mean.vals, 3)
  Samples <- n.cores*n.batches*batch.size
  
  if(save.all.LLs) {
    se.vals <- sqrt(var(vals-mean(vals)))/sqrt(length(vals))
    out <- c(mean.vals, AIC, Samples, se.vals)
    names(out) <- c('LL', 'AIC', 'Samples', 'SE.LL')
  } else {
    out <- c(mean.vals, AIC, Samples)
    names(out) <- c('LL', 'AIC', 'Samples')
  }
  
  attr(dclone.mod, 'stat.sum') <- out

  dclone.mod
}
