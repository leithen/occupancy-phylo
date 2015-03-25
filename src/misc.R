## create a balanced tree with branch lengths
balanced.tree <- function(nsp) {
  tree<-stree(n=nsp, tip.label=1:nsp , type='balanced')
  tree$edge.length <- node.depth(tree)[tree$edge[,2]]
  tree
}

## Quick test function to see how good ML is at estimating the true
## value of lambda depending on different tree topologies and
## different nsps
phy.sig.tester <- function(tree, lambda) {

  nsp <- length(tree$tip.label)
  vcov <- vcv(tree, corr=T)

  final.vcov <- vcov*lambda+(1-lambda)*diag(nsp) #calculate vcov based on lambda

  beta.spp <- mvrnorm(n=1, mu=rep(0, nsp), Sigma=final.vcov) #Simulate


  ## Quality model fitting method. 
  fC <- fitContinuous(tree, beta.spp, model="lambda")


  ## Other estimate of phylogenetic signal.
  ps <- phylosig(tree, beta.spp, method="lambda", test=F)

  return(data.frame(fitContinuous=fC$opt$lambda, phylosig=ps$lambda))

}



### Extracts model coeficients, and optionally creates a plot
pom.coefs<-function(dd.simulated, res, phy=NULL, parameter.group, plot=T, color.gradient=c("red", "blue", "orange"), tree.label.cex=0.7, cex.tip.labs=1 , english=F,...) {


if (!is.null(phy)) {
		tree<-phy
}

else {
	tree<-dd.simulated$inputs$tree
}

model.summary<-res$BUGSoutput$summary[,c('mean', '2.5%', '97.5%','Rhat') ]

psi.beta.sp.ID<-grep(parameter.group, row.names(model.summary))

coefs<-model.summary[psi.beta.sp.ID,]

coefs

######################################################
###Add Phylogeny
######################################################
require(ape)

#Load tree

if (!is.null(phy)) {	
		nm<-read.csv("~/Dropbox/ccb_banding/data/original/NameMatch_Master_WorkingVersion.csv")

		spp<-dimnames(dd.simulated$X)$species
		jetz<-nm[match(spp, nm$Chase) , 'Jetz']
		row.names(coefs)<-jetz
		
		noshows<-tree$tip.label[!tree$tip.label %in% jetz]
		tree<-drop.tip(tree, noshows)
	}
	
else {row.names(coefs)<-as.character(tree$tip.label)}


tip.labels.in.order<-tree$tip.label[tree$edge[,2][tree$edge[,2] %in% 1:dim(coefs)[1]]] #Have to reorder that slopes final plots in same order as tree

coefs.final <-coefs[match( tip.labels.in.order, row.names(coefs)),]

if (plot==TRUE) {
############################################################
###With Tip Labels
############################################################
par(mar=c(6,0,2,0))
#layout(matrix(c(1,2,3), nrow=1), widths=c(1,1,.2))
layout(matrix(c(1,2), nrow=1), widths=c(1, 1))


r.co<-range(coefs.final[,1])


colz.tree<-colorize(coefs.final[,1], colors= color.gradient, min=-max(abs(r.co)), max=max(abs(r.co)))

if(english){
	traits<-read.csv("~/Dropbox/ccb_banding/data/original/Traits/Bird_Traits_WorkingVersion.csv")
	spp.codes.tree<-nm$Jim[match(tree$tip.label, nm$Jetz)]
	tree$tip.label <-as.character(traits$english_name [match(spp.codes.tree, traits$species_code)])
}

plot(tree, show.tip.label=T, cex=tree.label.cex, label.offset=2)
tiplabels(pch=22, bg=colz.tree, cex=cex.tip.labs)

colz.plot<-colorize(coefs.final[,1], colors=color.gradient)

plot(y=1:dim(coefs.final)[1], x= coefs.final[,1], bty='n', yaxt='n', ylab="", pch=21, bg= colz.plot, ...)
arrows(y0=1:dim(coefs.final)[1], x0= coefs.final[,2], x1= coefs.final[,3], length=0)
abline(v=0, lty=2)
}

return(coefs.final)
}


## ## Try both of these
## tree <- balanced.tree(32)
## ## -or-
## tree <- rcoal(32)


## ## This function spits out two estimates of lambda, conducted using
## ## different commonly used ML functions in R. 

## phy.sig.tester(tree, lambda=0.7) 
## ## Performs terribly with a balanced tree. Less terribly with a random
## ## coalescent. But very innaccurate when nsp is small.
