
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


compile.bf.modruns<-function(directory) {

all.files<-list.files(directory)

output<-matrix(NA, nrow=length(all.files), ncol=10)

colnames(output)<-c("lambda.case", "nsp.case", "lambda.KNOWN", "lambda.KNOWN.P", "lambda.EST", "lambda.EST.mode", "lambda.EST.dens", "mod.weight", "bf", "rhat")

simz<-character()

for (i in 1:length(all.files)) {

	quer<-all.files[i]
	quer.final<-paste(directory,"/",quer, sep="")
	load(quer.final)

	lambda.case<-as.numeric(strsplit(quer, "_")[[1]][4])
	nsp.case<-as.numeric(strsplit(quer, "_")[[1]][6])
	sim.case<-strsplit(strsplit(quer, "_")[[1]][8], "\\.")[[1]][1]

	logit.lambda <- res.full$BUGSoutput$sims.matrix[,'logit.lambda']

	lambda.KNOWN<-data.frame(phylosig(dd.simulated$inputs$tree, dd.simulated$inputs	$beta.sp, 'lambda', test=T))[1]
	lambda.KNOWN.P<-data.frame(phylosig(dd.simulated$inputs$tree, dd.simulated$inputs	$beta.sp, 'lambda', test=T))[4]

	
	lambda.EST<-expit(res.full$BUGSoutput$summary['logit.lambda',c('mean', '2.5%', '97.5%')])[1]
	mod.weight<-tryCatch( bf.test(res.select)[1,1], error=function(e) {NA})
	bf<-tryCatch( bf.test(res.select)[1,2], error=function(e) {NA})	
	rhat<-tryCatch( bf.test(res.select)[1,5], error=function(e) {NA})

	lambda.EST.mode<-Mode(
		as.numeric(as.character(cut(expit(logit.lambda), breaks=seq(0,1, by=0.05), labels=seq(0,1, by=0.05)[-1]-0.025)))
	)

	denz<-density(expit(logit.lambda))
	lambda.EST.dens<-denz$x[denz$y==max(denz$y)]



	output[i,]<-unlist(c(lambda.case, nsp.case, lambda.KNOWN, lambda.KNOWN.P ,lambda.EST, lambda.EST.mode, lambda.EST.dens ,mod.weight, bf, rhat))
	simz[i]<-sim.case
}

final.df<-data.frame(sim.case=simz, output)
final.df$lbf<-2*log(final.df$bf)

return(final.df)
}


bf.test<-function(res.select) {

### Select
bf.pre <- res.select$BUGSoutput$summary['IND.lambda','mean'] #Pull out indicator posterior

bf <- bf.pre/(1-bf.pre) ### Calculate Bayes Factor (How many times more likely is the full model than the null?)

pseudoLR <- 2*log(bf) ### Convert Bayes Factor to something that looks like a likelihood ratio

pseudoP<-pchisq(pseudoLR, df=1, lower.tail=F) ### Do a completely silly 'LR' test. This doesn't actually make any sense, but helps me visualize how 'significant' the parameter is.


#traceplot(res.select)


rhat<-res.select$BUGSoutput$summary['IND.lambda','Rhat']

model.weight<- res.select$BUGSoutput$summary['IND.lambda','mean'] 



data.frame(Model.Weights=c(model.weight),BayesFactors=c(bf), LR=c(pseudoLR), P=c(pseudoP), Rhat=c(rhat), row.names=c("Select"))
}

### Old function to compare 'RJ' verusus single variable 'select' means of model selection
bf.comparison<-function(res.RJ, res.select) {
# RJ

bf.pre <- res.RJ$BUGSoutput$summary['IND','mean'] #Pull out indicator posterior

bf <- (1-bf.pre)/bf.pre ### Calculate Bayes Factor (How many times more likely is the full model than the null?)
bf
bf.RJ<-bf


pseudoLR <- 2*log(bf) ### Convert Bayes Factor to something that looks like a likelihood ratio
pseudoLR
pseudoLR.RJ<-pseudoLR


pseudoP<-pchisq(pseudoLR, df=1, lower.tail=F) ### Do a completely silly 'LR' test. This doesn't actually make any sense, but helps me visualize how 'significant' the parameter is.
pseudoP.RJ<-pseudoP

### Select
bf.pre <- res.select$BUGSoutput$summary['IND.lambda','mean'] #Pull out indicator posterior

bf <- bf.pre/(1-bf.pre) ### Calculate Bayes Factor (How many times more likely is the full model than the null?)
bf.sel<-bf

pseudoLR <- 2*log(bf) ### Convert Bayes Factor to something that looks like a likelihood ratio
pseudoLR
pseudoLR.sel<-pseudoLR

pseudoP<-pchisq(pseudoLR, df=1, lower.tail=F) ### Do a completely silly 'LR' test. This doesn't actually make any sense, but helps me visualize how 'significant' the parameter is.
pseudoP.sel<-pseudoP


#traceplot(res.select)


rhat.RJ<-res.RJ$BUGSoutput$summary['IND','Rhat']
rhat.sel<-res.select$BUGSoutput$summary['IND.lambda','Rhat']

model.weight.RJ<- 1- res.RJ$BUGSoutput$summary['IND','mean'] 
model.weight.sel<- res.select$BUGSoutput$summary['IND.lambda','mean'] 



data.frame(Model.Weights=c(model.weight.RJ, model.weight.sel),BayesFactors=c(bf.RJ, bf.sel), LR=c(pseudoLR.RJ, pseudoLR.sel), P=c(pseudoP.RJ, pseudoP.sel), Rhat=c(rhat.RJ, rhat.sel), row.names=c("RJ", "Select"))
}
