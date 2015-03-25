# This program contains the calculations needed to draw the profile likelihood for 
# the parameter "lambda" for the Stochastic Beverton-Holt (SBH) model with Poisson    
# Sampling error, described in the paper

# "Confidence Intervals, model selection and Likelihood ratio tests for data cloning". 
#  J.M. Ponciano, Mark L. Taper, Brian Dennis and Subhash R. Lele.  Ecology 90:356-362.

# The number of parameters in the SBH model is 3: lambda, b and sigma (the process error 
# variance). The data used in this file is the first of three replicated time series for 
# Paramecium aurelia, recorded by G.F. Gause in his 1934 book: "The struggle for existence".
# In the following code, WinBugs is called from R using the package "BRugs".
# WinBugs used in conjunction with BRugs in R is convenient and easy to code.
# The following small "txt" files with WinBugs code need to be saved in the same directory 
# than this program (their contents for the SBH model are shown at the end of this file):
#
# BevHoltPoisson1ts.txt:  Contains the "model" code in WinBugs needed to estimate the initial 
#                         ML estimates of all the parameters in the model:
#
# BevHoltPoisson1ts_Inits.txt:  Contains the parameter values used to initialize the chain for 
#                               the initial ML estimation. 
#
# LprofBevHoltPoisson1rep:  Contains the "model" code in WinBugs needed to estimate the restricted 
#                           ML estimates of the parameters 'b' and 'sigma', while profiling over
#                           values of the parameter lambda. 
#
# BevHoltSSM_Inits.txt:   Contains the parameter values used to initialize the chain for the 
#                         restricted maximizations
#
# hofXgYBevHoltPoisson1rep.txt:  Contains the  "model" code in WinBugs needed to sample from the 
#                                conditional P(X|Y,theta.hat) where theta.hat are the initial ML estimates
#
# 
# The calculations proceed in 5 sequential steps (Only the order of steps 2 and 3 can be interchanged):
# 
# 1. Finding the ML estimates for ALL the parameters in the models.  To do that, increase the number of clones  
#    and compare results. Here, the number of clones (k) was increased five times by a two-fold amount.   
#    The users should only proceed with the other calculations after confirming that convergence of the MCMC 
#    chains has been achieved. Failing to do so will return erroneous profiles. The number of MCMC iterations
#    needed in order to achieve convergence might be different from what is shown here depending on your system.
#
# 2. Find the restricted ML estimates for 'b' and 'sigma' for an array of 'lambda' values.
#
# 3. Take many (here, about 46000) samples from the posterior P(X|Y,theta.hat) to later compute the Monte Carlo 
#    Likelihood Ratio approximation.
#
# 4. Calculate the Monte Carlo likelihood ratio's denominator for every single sampled X from step 3.
#
# 5. Calculate the Monte Carlo likelihood ratio;s numerator for every single sampled X from step 3 and for all
#    the lambda values in the same array used in step 2. Once that is done, the entire profile is calculated
#    
#
#
# The only thing the users need to change from this code are all the path statements. 
# A sample of the output is interleaved with the code as commented text in some places.
# All the programs were written by Jose Miguel Ponciano in June 2007 and any concerns or questions
# should be addressed to him at ponciano@cimat.mx.


setwd("~/Dropbox/ccb_banding/analyses/multi-species-chase/src/MCLR")
library(BRugs)
library("MASS")


#############################################################################################
#######--------------------------------         STEP 1
#######---------------------Doing a profile likelihood using just 1 replicate 
#######--------------------- of the time series of P. aurelia:

kvec <- c(20,40,80,160,240);
for(i in 1:5){

	K <- kvec[i];

	Paurelia1 <- c(17,29,39,63,185,258,267,392,510,570,650,560,575,650,550,480,520,500);
	qp1       <- 18;
	tskreps1<- matrix(rep(Paurelia1,K),nrow=K,ncol=qp1,byrow=T);
	dataforBugs <- list(K=K,Y1=tskreps1);
	bugsData(data=dataforBugs,
	fileName=file.path("C:\\PONCIANO\\GauseDataCloning\\LikelihoodRatio","Paurelia10k1rep.txt"),digits=4);


	#######---------------- Stochastic Beverton-Holt model with Lognormal observation error---------------####
	modelCheck("BevHoltPoisson1ts.txt");         # check the model file
	modelData("Paurelia10k1rep.txt");            # read data file
	modelCompile(numChains=1);                   # compile model with 1 chain
	modelInits("BevHoltPoisson1ts_Inits.txt");   # read initial parameter values data file
	modelGenInits();                             # get initial guesses of the X process
	
	samplesSet(c("lam1","b1","sig1"));           # parms. that should be monitored
	samplesSetBeg(5000);
	samplesSetEnd(10000);
	modelUpdate(10000);          
	#all.nodes <- samplesHistory("*",plot=FALSE)
	#len.chains <- length(all.nodes[[1]]);
	#BevertonHolt.history1rep <-matrix(0,ncol=3,nrow=len.chains);
	#for(j in 1:3){BevertonHolt.history1rep[,j] <- all.nodes[[j]]};
	Sum.chainsBH.1rep   <-samplesStats("*",beg=5000,end=10000);	# Get summarized results -print them to take a look-
	MLES.BH.1rep        <- Sum.chainsBH.1rep[,1];
	SDMLEs.BH.1rep      <- sqrt(K)*(Sum.chainsBH.1rep[,2]);
	print(paste("K=",K));
	print(rbind(MLES.BH.1rep,SDMLEs.BH.1rep,Sum.chainsBH.1rep[,2]));
}

#K = 20
#
#10000 updates took 84 s
#                       [,1]      [,2]       [,3]
#MLES.BH.1rep   0.0021420000 2.1760000 0.12030000
#SDMLEs.BH.1rep 0.0002743655 0.1127873 0.03121551
#               0.0000613500 0.0252200 0.00698000
#
#K = 40
#
#10000 updates took 163 s
#                       [,1]      [,2]       [,3]
#MLES.BH.1rep   0.0021370000 2.1740000 0.11930000
#SDMLEs.BH.1rep 0.0002794189 0.1167513 0.03247027
#               0.0000441800 0.0184600 0.00513400
#
#K = 80
#10000 updates took 347 s
#                       [,1]      [,2]       [,3]
#MLES.BH.1rep   0.0021320000 2.1730000 0.11880000
#SDMLEs.BH.1rep 0.0002841595 0.1177961 0.03117079
#               0.0000317700 0.0131700 0.00348500
#
#K= 160
#                       [,1]      [,2]       [,3]
#MLES.BH.1rep   0.0021330000 2.1730000 0.11840000
#SDMLEs.BH.1rep 0.0002666433 0.1117169 0.03111681
#               0.0000210800 0.0088320 0.00246000
#
#K= 240
#
#10000 updates took 1177 s
#                       [,1]      [,2]       [,3]
#MLES.BH.1rep   0.0021320000 2.1730000 0.11860000
#SDMLEs.BH.1rep 0.0002752917 0.1154924 0.03186691
#               0.0000177700 0.0074550 0.00205700


MLES.BH.1rep <-  c(0.0021320000, 2.1730000, 0.11860000)


#######--------------------------------         STEP 2
#######------- Find the restricted ML estimates for 'b' and 'sigma' for an array of 'lambda' values.
K       <- 20;
lam.vec <- seq(1.81,2.80,by=0.01);
len.lam <- length(lam.vec);
prof.mles1rep <- matrix(0,nrow=len.lam,ncol=2)

for(i in 1:len.lam){

	dataforBugs <- list(K=K,lam=lam.vec[i],Y1=tskreps1);
	bugsData(data=dataforBugs,
	fileName=file.path("C:\\PONCIANO\\GauseDataCloning\\LikelihoodRatio","Paureliaprof10k1rep.txt"),digits=4)


	#######---------------- Stochastic Beverton-Holt model with Poisson observation error---------------####
	modelCheck("LprofBevHoltPoisson1rep.txt");     # check the model file
	modelData("Paureliaprof10k1rep.txt");         # read data file
	modelCompile(numChains=1);            # compile model with 1 chain
	modelInits("BevHoltSSM_Inits.txt");   # read initial parameter values data file
	modelGenInits();                      # get initial guesses of the X process
	
	samplesSet(c("b","sig"));       # parms. lam,b,sig, tau should be monitored
	samplesSetBeg(5000);
	samplesSetEnd(10000);
	modelUpdate(10000);          
	Sum.chainsBH   <- samplesStats("*",beg=5000,end=10000);	# Get summarized results -print them to take a look-
	MLES.BH        <- Sum.chainsBH[,1];
	prof.mles1rep[i,]  <- MLES.BH;
	#print(Sum.chainsBH)
}


#write.matrix(prof.mles1rep,file="profmlesPaureliaBH1rep.txt")


#######--------------------------------         STEP 3
#####----- After getting the MLes, sample from h(X|Y), tilting at the mles


dataforBugs <- list(lam=MLES.BH.1rep[2],b=MLES.BH.1rep[1],sig=MLES.BH.1rep[3],Y1=Paurelia1);
bugsData(data=dataforBugs,
fileName=file.path("C:\\PONCIANO\\GauseDataCloning\\LikelihoodRatio\\","PofXgYdatBHPoiss1rep.txt"),digits=5)
lenXs     <- 19; #(18+ 1 latent variable)

modelCheck("hofXgYBevHoltPoisson1rep.txt");      # check the model file
modelData("PofXgYdatBHPoiss1rep.txt");         # read data file
modelCompile(numChains=1);            # compile model with 1 chain
modelGenInits();                      # get initial guesses of the X process

samplesSet("X1"); # parms. a,b,sig, tau should be monitored
samplesSetBeg(40000);
samplesSetEnd(1940000);
samplesSetThin(40);
modelUpdate(1900000);                 
Sum.chains  <-samplesStats("*")	# Get summarized results          	
print(Sum.chains)
all.nodes <- samplesHistory("*",plot=FALSE)
len.chains <- length(all.nodes[[1]]);
out.name <- "C:\\PONCIANO\\GauseDataCloning\\LikelihoodRatio\\XgYBevHoltHistory1rep.txt";
out.history <-matrix(0,ncol=lenXs,nrow=len.chains);
for(j in 1:lenXs ){out.history[,j] <- all.nodes[[j]]};
write.matrix(out.history,file=out.name)

#> print(Sum.chains)
#        mean      sd  MC_error val2.5pc median val97.5pc start sample
#X1[1]  1.537 0.10190 0.0004999    1.338  1.537     1.738 40001  46500
#X1[2]  2.377 0.11420 0.0005451    2.152  2.377     2.599 40001  46500
#X1[3]  3.116 0.10790 0.0004998    2.903  3.116     3.325 40001  46500
#X1[4]  3.734 0.09531 0.0004823    3.545  3.735     3.919 40001  46500
#X1[5]  4.352 0.07992 0.0003958    4.194  4.352     4.507 40001  46500
#X1[6]  5.134 0.06355 0.0002939    5.009  5.135     5.257 40001  46500
#X1[7]  5.530 0.05495 0.0002332    5.421  5.531     5.637 40001  46500
#X1[8]  5.649 0.05175 0.0002556    5.547  5.649     5.749 40001  46500
#X1[9]  5.974 0.04564 0.0002223    5.884  5.975     6.062 40001  46500
#X1[10] 6.226 0.04110 0.0001758    6.146  6.227     6.307 40001  46500
#X1[11] 6.343 0.03923 0.0001864    6.266  6.344     6.420 40001  46500
#X1[12] 6.460 0.03752 0.0001763    6.385  6.460     6.532 40001  46500
#X1[13] 6.334 0.03932 0.0001828    6.257  6.335     6.411 40001  46500
#X1[14] 6.356 0.03891 0.0001898    6.279  6.357     6.432 40001  46500
#X1[15] 6.459 0.03760 0.0001709    6.384  6.459     6.532 40001  46500
#X1[16] 6.310 0.03969 0.0001863    6.233  6.311     6.388 40001  46500
#X1[17] 6.190 0.04169 0.0001980    6.107  6.190     6.271 40001  46500
#X1[18] 6.250 0.04070 0.0001880    6.169  6.250     6.328 40001  46500
#X1[19] 6.222 0.04193 0.0002071    6.140  6.222     6.303 40001  46500


#######--------------------------------         STEP 4
#######---------------- Then compute the denominator from E Thompson's ratio

#### Calculating EThompson's MC Like ratio: first the denominator:

B         <- len.chains;
b.max     <- MLES.BH.1rep[1];
lam.max   <- MLES.BH.1rep[2];
sig.max   <- MLES.BH.1rep[3];
X1.star1  <- out.history[,1:19];
Y1        <- Paurelia1;
logpxy.hat<- rep(0,B);
lenY      <- length(Y1);

#### Computing the denominator of EThompson's Monte Carlo Likelihood Ratio
for(i in 1:B){
	
	Xloglikes.hat      <- rep(0,lenY+1);
	
	for(j in 2:(lenY+1)){
			Xloglikes.hat[j] <- dlnorm(
			exp(X1.star1[i,j]),
			meanlog = X1.star1[i,(j-1)] + log(lam.max) - log(1+b.max*exp(X1.star1[i,(j-1)])) , sdlog=sig.max,log=T);
	}
	X.logL  <- sum(Xloglikes.hat);
	Y.logL <- sum(dpois(Y1,lambda=exp(X1.star1[i,(2:19)]),log=T)); # Since the sampling model is Poisson and  
								       # there are no sampling model extra parameters, this line cancels 
								       # with the one on the numerator (see next part), but I include it here to 
								       # recall users that in general it does not cancel out!!
	logpxy.hat[i] <- Y.logL + X.logL; # log(P(y|x)P(X|theta.hat))  
}


#######--------------------------------         STEP 5
#######------ Finally, loop over different values of 'lambda', and calculate EThompson's ratio.

##### Now a double loop to calculate the numerator for each value of 'lam' from the profile
##### and the corresponding mles of 'b' and 'sigma' at those lambda values:

Like.profile1rep <- rep(0,len.lam);
half.indexes <- seq(1,len.lam,2);
for(h in half.indexes){

	lam.hat       <- lam.vec[h];
	b.hat         <- prof.mles1rep[h,1];
	sig.hat       <- prof.mles1rep[h,2];
	logratio.Prof <- rep(0,B);

	for(i in 1:B){
	
		Xloglikes.star      <- rep(0,lenY+1);
	
		for(j in 2:(lenY+1)){
			Xloglikes.star[j] <- dlnorm(exp(X1.star1[i,j]),meanlog=X1.star1[i,(j-1)]+ log(lam.hat)
			-log(1+b.hat*exp(X1.star1[i,(j-1)])),sdlog=sig.hat,log=T);
		}
		X.logL  <- sum(Xloglikes.star);
		Y.logL <- sum(dpois(Y1,lambda=exp(X1.star1[i,(2:19)]),log=T));
	
		#print(c(Y.logL+X.logL, -logpxy.hat[i]))
		logratio.Prof[i] <- Y.logL + X.logL - logpxy.hat[i]; # log(P(y|x)P(X|theta.hat))
	}

	Like.profile1rep[h] <- mean(exp(logratio.Prof));
	
}

plot(lam.vec[half.indexes],Like.profile1rep[half.indexes],type="l",col="red", lwd=2, 
ylab="Likelihood score relative to the maximum", xlab=expression(lambda), 
main="Profile likelihood", cex=1.5,cex.lab=1.5);


###  OK: re-tilting with better parameter values:

max.index <- which(Like.profile1rep==max(Like.profile1rep),arr.ind=T);
lam.max   <- lam.vec[max.index];
b.max     <- prof.mles1rep[max.index,1];
sig.max   <- prof.mles1rep[max.index,2];


print(c(b.max,lam.max,sig.max))
#> print(c(b.max,lam.max,sig.max))
#[1] 0.002047 2.130000 0.120300
#> 
#> MLES.BH.1rep
#[1] 0.002132 2.173000 0.118600

dataforBugs <- list(lam=lam.max,b=b.max,sig=sig.max,Y1=Paurelia1);
bugsData(data=dataforBugs,
fileName=file.path("C:\\PONCIANO\\GauseDataCloning\\LikelihoodRatio\\","PofXgYdatBHPoiss1rep.txt"),digits=5)
lenXs     <- 19; #(18+ 1 latent variable)

modelCheck("hofXgYBevHoltPoisson1rep.txt");      # check the model file
modelData("PofXgYdatBHPoiss1rep.txt");         # read data file
modelCompile(numChains=1);            # compile model with 1 chain
modelGenInits();                      # get initial guesses of the X process

samplesSet("X1"); # parms. a,b,sig, tau should be monitored
samplesSetBeg(40000);
samplesSetEnd(1940000);
samplesSetThin(40);
modelUpdate(1900000);                 
Sum.chains  <-samplesStats("*")	# Get summarized results          	
print(Sum.chains)
all.nodes <- samplesHistory("*",plot=FALSE)
len.chains <- length(all.nodes[[1]]);
out.name <- "C:\\PONCIANO\\GauseDataCloning\\LikelihoodRatio\\XgYBevHoltHistory1rep2ndtour.txt";
out.history <-matrix(0,ncol=lenXs,nrow=len.chains);
for(j in 1:lenXs ){out.history[,j] <- all.nodes[[j]]};
#write.matrix(out.history,file=out.name);



B         <- len.chains;
X1.star1  <- out.history[,1:19];
Y1        <- Paurelia1;
logpxy.hat2ndtour<- rep(0,B);
lenY      <- length(Y1);

#### Computing the denominator of EThompson's Monte Carlo Likelihood Ratio
for(i in 1:B){
	
	Xloglikes.hat      <- rep(0,lenY+1);
	
	for(j in 2:(lenY+1)){
			Xloglikes.hat[j] <- dlnorm(exp(X1.star1[i,j]),meanlog=X1.star1[i,(j-1)]+ log(lam.max)
			-log(1+b.max*exp(X1.star1[i,(j-1)])),sdlog=sig.max,log=T);
	}
	X.logL  <- sum(Xloglikes.hat);
	Y.logL <- sum(dpois(Y1,lambda=exp(X1.star1[i,(2:19)]),log=T));
	
	logpxy.hat2ndtour[i] <- Y.logL + X.logL; # log(P(y|x)P(X|theta.hat))
}


####------ Finally, loop over different values of 'lambda', and calculate EThompson's ratio 
####------ using previous result

##### Now a double loop to calculate the numerator for each value of 'lam' from the profile
##### and the corresponding mles of 'b' and 'sigma' at those lambda values:

Like.profile1rep2ndtour <- rep(0,len.lam);
half.indexes <- seq(1,len.lam,2);
for(h in half.indexes){

	lam.hat       <- lam.vec[h];
	b.hat         <- prof.mles1rep[h,1];
	sig.hat       <- prof.mles1rep[h,2];
	logratio.Prof <- rep(0,B);

	for(i in 1:B){
	
		Xloglikes.star      <- rep(0,lenY+1);
	
		for(j in 2:(lenY+1)){
			Xloglikes.star[j] <- dlnorm(exp(X1.star1[i,j]),meanlog=X1.star1[i,(j-1)]+ log(lam.hat)
			-log(1+b.hat*exp(X1.star1[i,(j-1)])),sdlog=sig.hat,log=T);
		}
		X.logL  <- sum(Xloglikes.star);
		Y.logL <- sum(dpois(Y1,lambda=exp(X1.star1[i,(2:19)]),log=T));
	
		#print(c(Y.logL+X.logL, -logpxy.hat[i]))
		logratio.Prof[i] <- Y.logL + X.logL - logpxy.hat2ndtour[i]; # log(P(y|x)P(X|theta.hat))
	}

	Like.profile1rep2ndtour[h] <- mean(exp(logratio.Prof));
	
}

plot(lam.vec[half.indexes],Like.profile1rep2ndtour[half.indexes],type="l",col="red", 
lwd=2, ylab="Likelihood score relative to the maximum", xlab=expression(lambda), 
main="Profile likelihood", cex=1.5,cex.lab=1.5);

#### ---- END OF LIKELIHOOD PROFILE CALCULATION PROGRAM-----

#### ----  IN WHAT FOLLOWS THE CONTENTS OF THE VARIOUS "BUGS" PROGRAM FILES CALLED ABOVE ARE SHOWN
#### ----     please save each of these 5 code sections as separate 'txt' files

# 1. BevHoltPoisson1ts.txt:  Contains the "model" code in WinBugs needed to estimate the initial 
#                         ML estimates of all the parameters in the model:
#### BevHoltPoisson, a single time series
model{
for(k in 1:K){
	# Updating the state:
	X1[k,1]~dnorm(meano1,varE1);
	
	for(i in 2:19){
		X1[k,i]~dnorm(meanX1[k,i],varE1);
		meanX1[k,i] <- X1[k,(i-1)]+ log(lamm11+1)-log(1+b1*exp(X1[k,(i-1)]));	
	}

	# Updating the observations
	for(i in 1:18){
		Y1[k,i]~dpois(mu1[k,i]);
		log(mu1[k,i]) <- X1[k,(i+1)];

	}
}
# Priors on model parameters
lamm11~dlnorm(0,1);
b1~dlnorm(-1,1); 
sig1~dunif(0,1);
lam1 <- lamm11 + 1;
meano1 <- log(2)  + log(lamm11+1)-log(1+b1*2);
varE1        <- 1/pow(sig1,2);

}

# 2. BevHoltPoisson1ts_Inits.txt:  Contains the parameter values used to initialize the chain for 
#                               the initial ML estimation.

list(lamm11=1.2,b1=0.13,sig1=0.174)

# 3. LprofBevHoltPoisson1rep:  Contains the "model" code in WinBugs needed to estimate the restricted 
#                           ML estimates of the parameters 'b' and 'sigma', while profiling over
#                           values of the parameter lambda. 


# 4. BevHoltSSM_Inits.txt:   Contains the parameter values used to initialize the chain for the 
#                         restricted maximizations
#
list(b=0.13,sig=0.174)


# 5. hofXgYBevHoltPoisson1rep.txt:  Contains the  "model" code in WinBugs needed to sample from the 
#                                conditional P(X|Y,theta.hat) where theta.hat are the initial ML estimates
model{

# Updating the observations
for(i in 1:18){
	Y1[i]~dpois(mu1[i]);
	log(mu1[i]) <- X1[(i+1)];

}

# Priors on the X's
X1[1]~dnorm(meano,varE);
	
for(i in 2:19){
	X1[i]~dnorm(meanX1[i],varE);
	meanX1[i] <- X1[(i-1)]+ log(lam)-log(1+b*exp(X1[(i-1)]));	
		
}

# Constants definitions
meano <- log(2)  + log(lam)-log(1+(b*2));
varE        <- 1/pow(sig,2);
}

