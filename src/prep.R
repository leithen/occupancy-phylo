prep.data  <-  function(mat.path, case, threshold, drop.migrants) {
  browser()
  ## load matrices
  load(mat.path)
  mm <- make.mat(mat, site.list(case), threshold=threshold, nzero=0)
  mat <- mm$mat
  dates <- mm$dates

  ## make 4D matrix
  mats <- lapply(1:dim(mat)[2], function(x) mat[,x,])
  mats.split <- split(mats, dimnames(mat)$date)
  yr.table <- table(dimnames(mat)$date)
  X <- array(NA, dim=c(dim(mat)[1], length(yr.table),
                   max(yr.table), dim(mat)[3]))
  dimnames(X) <- list(site=dimnames(mat)$site,
                      year=unique(dimnames(mat)$date),
                      rep=1:max(yr.table),
                      species=dimnames(mat)$species)
  
  null.mat <- matrix(NA, dim(mat)[1], dim(mat)[3],
                     dimnames=dimnames(X)[c('site', 'species')])
  f <- function(i) {
    missing <- max(yr.table)-yr.table[i]
    if(missing==0) return(mats.split[[i]])
    c(mats.split[[i]], lapply(1:missing, function(x) null.mat))
  }
  tmp <- lapply(seq_along(yr.table), f)
  for(i in 1:length(yr.table))
    for(j in 1:max(yr.table))
      X[,i,j,] <- tmp[[i]][[j]]
  
  ## only keep sites with some positive number of reps
  no.reps <- apply(mat, 1, function(x) sum(x>=0, na.rm=TRUE))==0
  X <- X[!no.reps,,,,drop=FALSE]

  date.mat <- array(NA, dim=dim(X)[1:3], dimnames=dimnames(X)[1:3])
  for(i in seq_along(dates)) {
    year <- as.numeric(format(as.Date(dates[[i]], format='%Y-%m-%d'),
                              format = '%Y'))
    lengths <- rle(as.vector(year))$lengths
    ind <- cbind(rep(i, sum(lengths)),
                 match(year, dimnames(X)$year),
                 as.vector(unlist(sapply(lengths, seq_len))))
    date.mat[ind] <- strptime(dates[[i]], '%Y-%m-%d')$yday+1 
  }
  
  date.mat <- date.mat-min(date.mat, na.rm=TRUE)
  date.mat <- date.mat/max(date.mat, na.rm=TRUE)

  ## --- drop migrants ---
  if(drop.migrants==TRUE) {
    dd.migratory  <- 
      read.csv('data/original/chase/covariates/Migrant List.csv',
               as.is=TRUE)
    ## keep 'non-migrants'
    X <- X[,,,!dimnames(X)$species %in% dd.migratory$Latitudinal.migrants]
  }

  ## --- site covariates ---
  path <- 'data/original/chase/covariates'
  dd.site  <- 
    read.csv(file.path(path, 'site_specific_spatial_covariates.csv'),
             row.names=NULL, as.is=TRUE)
  drop.site <- is.na(match(dimnames(X)$site, dd.site$Site.Code))
  
  X <- X[!drop.site,,,]
  date.mat <- date.mat[!drop.site,,]

  cover.cols  <- 
    sapply(1:20*50, function(x)
           sprintf('%s.%d.%s',
                   'Percent.forest.cover.within',
                   x, 'meter.radius.of.center'))

  landcover <- dd.site[match(dimnames(X)$site, dd.site$Site.Code),
                       cover.cols]
  colnames(landcover) <- sapply(1:20*50,
                                function(x) sprintf('fc.%d', x))
  
  ## --- species covariates ---
  dd.traits <-
    read.csv('data/original/Traits/Bird_Traits_WorkingVersion_31-07-2014.csv')
  drop.sp <- is.na(match(dimnames(X)$species, dd.traits$ChaseNames))
  
  if(any(drop.sp)) X <- X[,,,!drop.sp]

  ## subset and re-order trait data
  dd.traits <- dd.traits[dd.traits$ChaseNames %in%
                         dimnames(X)$species,] 
  dd.traits <- dd.traits[match(dd.traits$ChaseNames,
                               dimnames(X)$species),] 

  ## cts and ordinal
  body.mass <- log(dd.traits$body_mass_male)
  clutch.size <- log(dd.traits$Clutch.Size)
  diet.breadth <- dd.traits$DietBreadth
  breeding.season.length <- dd.traits$Breeding_season_length..mo.
  ## binary
  solitary <- dd.traits$solitary
  aerial <- dd.traits$aerial
  berry.fruit <-
    dd.traits$berry.fruit.plucking_.from_trees.bushes_above_ground.__rarely
  seeds <- dd.traits$seeds
  
  ## add names to species-level covariates
  names(body.mass) <- dimnames(X)$species
  names(clutch.size) <- dimnames(X)$species
  names(diet.breadth) <- dimnames(X)$species
  names(breeding.season.length) <- dimnames(X)$species
  names(solitary) <- dimnames(X)$species
  names(aerial) <- dimnames(X)$species
  names(berry.fruit) <- dimnames(X)$species
  names(seeds) <- dimnames(X)$species
  
  ## --- other ---
  year <- standardize(as.numeric(dimnames(date.mat)$year))

  ## specify the initial values
  Z <- apply(X, c(1,2,4),
             function(x) (sum(x,na.rm=TRUE)>0)*1)
  Z[apply(X, c(1,2,4), function(x) !any(!is.na(x)))] <- NA

  ## create site x species matrix indicating site presence
  site.presence <- (apply(X, c(1,4), sum)>0)*1

  ## create site x species matrix indicating fraction of times a
  ## species was present
  frac.presence.num <- apply(X, c(1,4), sum, na.rm=TRUE)
  frac.presence.denom <- apply(X, c(1,4), function(x) sum(!is.na(x)))
  frac.presence <- frac.presence.num/frac.presence.denom
  pres.yr.1 <- (apply(X[,1,,], c(1,3), sum, na.rm=TRUE)>0)*1
  frac.presence[pres.yr.1==1] <- 1

  ## load and sync up VCOV matrix
  load('data/original/chase/ChaseVCOV.rdata')
  ind <- match(dimnames(X)$species,
               dimnames(ChaseVCOV)[[1]])
  keep <- dimnames(X)$species[!is.na(ind)]
  
  ## subset to birds we have phylogenetic data for, and rescale so
  ## root of phylogeny is 0, and round to five digits
  VCOV <- ChaseVCOV[keep,keep]
  VCOV <- round((VCOV-min(VCOV))/(1-min(VCOV)), 5)
  
  X <- X[,,,keep]
  Z <- Z[,,keep]
  site.presence <- site.presence[,keep]
  frac.presence <- frac.presence[,keep]
  body.mass <- body.mass[keep] 
  clutch.size <- clutch.size[keep]
  
  list(X=X,
       Z=Z,
       VCOV=VCOV,
       ID=diag(1, dim(X)[4]),
       VV.mu=matrix(0,ncol=dim(X)[4], nrow=2),
       nsite=dim(X)[1],
       nyr=dim(X)[2],
       nrep=apply(X, c(1,2,4),
         function(x) sum(x>=0,na.rm=TRUE)),
       nsp=dim(X)[4],
       date.mat=date.mat,
       site.presence=site.presence,
       frac.presence=frac.presence,
       landcover=landcover,
       body.mass=body.mass,
       clutch.size=clutch.size,
       diet.breadth=diet.breadth,
       breeding.season.length=breeding.season.length,
       solitary=solitary,
       aerial=aerial,
       berry.fruit=berry.fruit,
       seeds=seeds)
}

## below function is only to be used when working with actual data
prep.for.analysis <- function(dd,
                              model.case,
                              trait='no.covariates') {

  if(!trait %in% c('no.covariates')) {
    ## drop species with missing trait values
    sp.to.drop <- id((sapply(trait, function(x) which(is.na(dd[[x]])))))
    if(length(sp.to.drop)>0) {
      dd$X <- dd$X[,,,-sp.to.drop]
      dd$Z <- dd$Z[,,-sp.to.drop]
      dd$nsp <- dd$nsp-length(sp.to.drop)
      dd$frac.presence <- dd$frac.presence[-sp.to.drop,]
      dd$site.presence <- dd$site.presence[-sp.to.drop,]

      dd$VCOV  <- dd$VCOV[-sp.to.drop,-sp.to.drop]
      dd$ID    <- dd$ID[-sp.to.drop,-sp.to.drop]
      dd$VV.mu <- dd$VV.mu[,-sp.to.drop]
      for(i in trait)
        dd[[i]] <- dd[[i]][-sp.to.drop]
    }
  }
  dd$trait <- dd[[trait]]
  
  ## update naming scheme
  dd$env <- as.numeric(dd$landcover[,'fc.100'])
  
  ## env for categorical cases
  if(model.case=='categorical' |
     length(grep("categorical", model.case))==1 |
     model.case=='categorical.indicator' |
     model.case=='categorical.fixed' |
     model.case=='categorical.fixed.BM' |
     model.case=='categorical.fixed.WN' |
     model.case=='categorical.fixed.LAM' |
     model.case=='categorical.phyloVCOV' |
     model.case=='categorical.phyloVCOV.2' |
     model.case=='categorical.nophylosig.2' |
     model.case=='categorical.phyloVCOV.4') {
    dd$env <- (dd$env>=0.5)*1+1
    agr.presence <- (apply(dd$X[dd$env==1,,,], 4, sum)>0)*1
    for.presence <- (apply(dd$X[dd$env==2,,,], 4, sum)>0)*1
    dd$site.type.presence <- rbind(agr.presence, for.presence)
  }

  dd
}

analyse <- function(dd) {

  Z <- dd$Z
  my.inits <- function() {
    list(Z=Z)
  }

  my.params <- get.params()

  dd <- list(data=dd, inits=my.inits, params=get.params())
  
  bugs <- run.R2jags.model(dd)
  list(data=dd,
       bugs=bugs,
       summary=bugs$BUGSoutput$summary)
}
