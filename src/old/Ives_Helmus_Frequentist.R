

setwd("~/Dropbox/ccb_banding/analyses/multi-species-chase/saved/simulated/Frequentist Model")
library(picante)


pglmm.data.fix<-function(samp, tree=NULL, env=NULL, modelflag=c(1,2), nudge=FALSE, VCOV=NULL) {
		
		if (is.null(tree)==FALSE) {
		tree <- prune.sample(samp, tree)
        samp <- samp[, tree$tip.label]
        V <- vcv.phylo(tree, corr = TRUE)
        }
        
        if(is.null(VCOV)==FALSE) {  	
        	V<-VCOV
        	samp <- samp[,colnames(V)]
        	}
        
        species <- colnames(samp)
        preval <- colSums(samp)/sum(sum(samp))
        species <- species[preval > 0]
        V <- V[species, species]
        if(nudge==TRUE) {
        	nslots<-length(V[upper.tri(V)])
        	add<-runif(nslots,min=0,max=0.0001)
			V[upper.tri(V)]<-V[upper.tri(V)]+add
			V[lower.tri(V)]<-t(V)[lower.tri(V)]       	
        	}
        
        #Vcomp <- Vcomp[species, species]
        samp <- samp[, colnames(V)]
        #traits <- as.matrix(traits[species, ])
    
   
 	Y <- samp
    #X <- traits
    nsites <- dim(Y)[1]
    nspp <- dim(Y)[2]
    YY <- t(Y)
    YY <- as.matrix(as.vector(as.matrix(YY)))
    if(is.null(env)==FALSE) {U <- matrix(env, nrow = length(env), ncol = 1)}



    if (modelflag == 1) {
        Vfullspp <- kronecker(diag(nsites), V)
        Vfullsite <- kronecker(diag(nsites), matrix(1, nspp, 
            nspp))
        VV <- list(Vfullspp = Vfullspp, Vfullsite = Vfullsite)
        XX <- kronecker(matrix(1, nsites, 1), diag(nspp))
        return(list(YY = YY, VV = VV, XX = XX))
    }
    if (modelflag == 2) {
        
        u <- scale(U)
        U <- kronecker(u, matrix(1, nspp, 1))
        Vfullspp <- kronecker(matrix(1, nsites, nsites), diag(nspp))
        VfullsppV <- kronecker(matrix(1, nsites, nsites), V)
        VfullUCU <- diag(as.vector(U)) %*% Vfullspp %*% diag(as.vector(U))
        VfullUCUV <- diag(as.vector(U)) %*% VfullsppV %*% diag(as.vector(U))
        Vfullsite <- kronecker(diag(nsites), matrix(1, nspp, 
            nspp))
        VV <- list(VfullUCU = VfullUCU, VfullUCUV = VfullUCUV, 
            Vfullsite = Vfullsite)
        XXspp <- kronecker(matrix(1, nsites, 1), diag(nspp))
        XX <- cbind(U, XXspp)
        YY <- as.vector(t(Y))
        return(list(YY = YY, VV = VV, XX = XX))
    }
  }

pglmm.fit.tracker<-function (dat = NULL, Y = NULL, X = NULL, VV = NULL, sp.init = 0.5, 
    maxit = 25, exitcountermax = 50) 
{
    PGLMM.reml <- function(sp) {
        sp <- abs(Re(sp))
        Cd <- matrix(0, dim(tVV[[1]])[1], dim(tVV[[1]])[2])
        for (i in 1:length(sp)) {
            Cd = Cd + sp[i] * tVV[[i]]
        }
        V <- tinvW + Cd
        invV <- solve(V, diag(x = 1, nrow = dim(V)[1], ncol = dim(V)[2]))
        if (all(eigen(V)$values > 0)) {
            cholV <- chol(V)
            LL = 0.5 * (2 * sum(log(diag(cholV))) + t(tH) %*% 
                invV %*% tH + log(det(t(tX) %*% invV %*% tX)))
        }
        else {
            LL <- 10^10
        }
        LL
    }
    if (!require(corpcor)) {
        stop("The 'corpcor' package is required")
    }
    if (!is.null(dat)) {
        X <- dat$XX
        Y <- dat$YY
        VV <- dat$VV
    }
    is.empty <- function(x) {
        length(x) == 0
    }
    if (any(unlist(lapply(list(X = X, Y = Y, VV = VV), is.empty)))) {
        stop("a data matrix is empty")
    }
    if (any(unlist(lapply(list(X = X, Y = Y, VV = VV), is.na)))) {
        stop("a data matrix is NA")
    }
    n <- dim(X)[1]
    p <- dim(X)[2]
    sp <- matrix(sp.init, length(as.list(VV)), 1)
    B0 <- matrix(mean(Y), p, 1)
    oldB0 <- matrix(10, p, 1)
    counter <- 0
    while (matrix((t(oldB0 - B0) %*% (oldB0 - B0) > 10^-8)) & 
        (counter < 100)) {
        oldB0 <- B0
        counter <- counter + 1
        mu <- exp(X %*% B0)/(1 + exp(X %*% B0))
        W <- as.vector((mu * (1 - mu))^-1)
        invW <- (diag(W))
        Z <- (X %*% B0 + (Y - mu)/(mu * (1 - mu)))
        denom <- (t(X) %*% invW %*% X)
        options(warn = -1)
        if (any(c(is.nan(denom), is.infinite(denom), is.null(denom)))) {
            B0 <- solve((t(X) %*% X), (t(X) %*% Y))
            (counter <- 100)
        }
        else {
            num <- (t(X) %*% invW %*% Z)
            B0 <- pseudoinverse(denom) %*% num
        }
    }
    if (is.nan(B0) | is.infinite(B0)) {
        mu <- matrix(mean(Y), p, 1)
        B0 <- log(mu/(1 - mu))
    }
    B <- B0
    b <- matrix(0, n, 1)
    bet <- (rbind(B, b))
    mu <- (exp(X %*% B)/(1 + exp(X %*% B)))
    XX <- cbind(X, diag(n))
    Cdum <- matrix(0, n, n)
    for (i in 1:length(sp)) {
        Cdum = Cdum + (sp[i] * VV[[i]])
    }
    est <- t(rbind(sp, B))
    oldest <- matrix(10^6, dim(est)[1], dim(est)[2])
    exitcounter <- 0
    while (matrix(((est - oldest) %*% t(est - oldest)) > 10^-4) & 
        exitcounter <= exitcountermax) {
        oldest <- est
        W <- as.vector((mu * (1 - mu))^-1)
        invW <- (diag(W))
        V <- invW + Cdum
        invV <- solve(V, diag(n))
        Z <- X %*% B + b + (Y - mu)/(mu * (1 - mu))
        denom <- (t(X) %*% invV %*% X)
        num <- (t(X) %*% invV %*% Z)
        B <- pseudoinverse(denom) %*% num
        b <- Cdum %*% invV %*% (Z - X %*% B)
        bet <- (rbind(B, b))
        mu <- exp(XX %*% bet)/(1 + exp(XX %*% bet))
        tH <- Z - X %*% B
        tinvW <- diag(as.vector((mu * (1 - mu))^-1))
        tX <- X
        tVV <- VV
        sp <- abs(optim(sp, PGLMM.reml, control = list(maxit = maxit, 
            abstol = 10^-1))$par)
        if (exitcounter == 0) {
            scrn.output <- c(exitcounter, t(sp), t(B[1:4]))
            names(scrn.output) <- c("exitcounter", paste("sigma", 
                1:length(sp)), paste("B", 1:4))
          
            Lrun <- PGLMM.reml(sp)

           tracker<-c(t(sp), t(B) ,Lrun)
			names(tracker)<-c(paste("sigma", 1:length(sp)), paste("B", 1:length(B), sep=""), "LL")
            print(scrn.output)
        }
        else {
            print(c(exitcounter, t(sp), t(B[1:4])))
        	Lrun <- PGLMM.reml(sp)
         	tracker<-rbind(tracker, c(t(sp), t(B) ,Lrun))
        	save(tracker, file="temporary.tracker.rdata")
        }
        Cdum <- matrix(0, n, n)
        for (i in 1:length(sp)) {
            Cdum = Cdum + sp[i] * tVV[[i]]
        }
        est <- t(rbind(sp, B))
        exitcounter <- exitcounter + 1
    }
    flag <- "converged"
    if (exitcounter >= exitcountermax) {
        flag <- "did not converge, try increasing exitcountermax"
    }
    W <- as.vector((mu * (1 - mu))^-1)
    invW <- (diag(W))
    V <- invW + Cdum
    invV <- solve(V, diag(n))
    Z <- X %*% B + b + (Y - mu)/(mu * (1 - mu))
    B <- solve((t(X) %*% invV %*% X), (t(X) %*% invV %*% Z))
    Chi025 <- 5.0238
    Chi05 <- 3.8414
    S95int <- NULL
    S <- sp
    for (iss in 1:length(sp)) {
        Lmin <- PGLMM.reml(S)
        if (sp[iss] > 0.02) {
            Smin <- optimize(PGLMM.reml, c(0, 0.9 * sp[iss]), 
                tol = 1e-04)$minimum
        }
        else {
            Smin <- 0
        }
        if (Smin < 0.01) {
            Smin <- 0
        }
        if (sp[iss] > 0.02) {
            Smax <- optimize(PGLMM.reml, c(1.1 * sp[iss], 10 * 
                sp[iss]), tol = 1e-04)$minimum
        }
        else {
            Smax <- optimize(PGLMM.reml, c(0.1, 5), tol = 1e-04)$minimum
        }
        S95int <- rbind(S95int, c(abs(Smin), abs(Smax)))
    }
    names(sp) <- names(VV)
    colnames(S95int) <- c("0.05", "0.95")
    return(list(B = B, B0 = B0, s = cbind(sp, S95int), LL = Lmin, 
        flag = flag))
}



set.seed(1)
nsp<-64
tree <- balanced.tree(nsp)
dd.model <- make.data(nsp=nsp, lambda=1, tree=tree, mu.p.0=0, sigma.beta.sp=10, det.image=T, psi.image=T, occ.image=T)

str(dd.model)

out<-list()
for (i in 1:5){
	mat3<-apply(dd.model$X, c(1,2,4), sum)
	mat2<-apply(dd.model$X, c(1,4), sum)
	mat2[mat2>0]<-1

	mat3.1<-mat3[,i,]
	mat3.1[mat3.1>0]<-1


	tree$tip.label<-dimnames(dd.model$X)$species

	pglmmDATA1<-pglmm.data.fix(samp=mat3.1, tree=tree, env=dd.model$landcover, modelflag=2, nudge=FALSE)

	str(pglmmDATA1)

	out[[i]]<-pglmm.fit.tracker(dat= pglmmDATA1)
	print(paste(date(), "completed", i))
}


load('~/Dropbox/ccb_banding/analyses/multi-species-chase/saved/simulated/Frequentist Model/temporary.tracker.rdata', verbose=T)

tracker

par(mar=c(5,5,1,1))
par(mfrow=c(3,3))
X<-40
for (i in 1:9) {
plot(tracker[,X-1+i], type="l", ylab=colnames(tracker)[X-1+i])
}


grabber<-function(list,item) {list[item]}

matrix(unlist(lapply(out, grabber, item='B')), ncol=5)


