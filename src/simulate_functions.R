## create data-sets with estimated lambda near the actual lambda
create.trees <- function(save.dir,
                         param,
                         values,
                         nsp) {
  save.path <-
    file.path('analyses/occupancy-phylo/saved/simulated',
              param, save.dir)
  if(!file.exists(save.path)) dir.create(save.path)

  create.data <- function(val) {
    create <- function() {
      tree <- balanced.tree(nsp)
      if(param=='lambda') {
        dd.model <- make.data(nsp=nsp,
                              lambda=val,
                              tree=tree)
        dd.model$inputs$lambda.estimated <- 
          phylosig(dd.model$inputs$tree,
                   dd.model$inputs$beta.spp,
                   method="lambda", test=F)$lambda
      }
      if(param=='sigma.psi.beta') {
        dd.model <- make.data(nsp=nsp,
                              sigma.psi.beta=val,
                              tree=tree)
        ## dd.model$inputs$sigma.psi.beta.estimated <- xxx
      }
      dd.model
    }
    dd.model <- create()

    if(param=='lambda')
      while(abs(dd.model$inputs$lambda.estimated - val)>0.05)
        dd.model <- create()

    if(param=='sigma.psi.beta')
      while(abs(dd.model$inputs$sigma.psi.beta - val)>0.05)
        dd.model <- create()

    to.print <- sprintf('actual lambda = %s, estimated lambda = %s\n',
                        dd.model$inputs$lambda,
                        dd.model$inputs$lambda.estimated)
    cat(to.print)

    fn <- sprintf('%s/%2.3f.RData', save.path, val)
    save(dd.model, file=fn)
  }

  sapply(values, create.data)
}

## run and save for various values of parameter 'param'
vary.param.dclone <- function(save.dir,
                              param,
                              values,
                              n.iter,
                              n.adapt,
                              n.update,
                              n.clones,
                              nsp,
                              cores, ...) {

  ## save directory (and create if the latter does
  ## not exist) 
  path <- file.path('analyses/occupancy-phylo/saved/simulated',
                    param, save.dir)
  if(!file.exists(path)) dir.create(path)

  ## function to run analysis across reps
  run.analysis <- function(x) {
    load.path <-
      file.path(path, 'trees')
    load(sprintf('%s/%2.3f.RData', load.path, x))

    save.path <- file.path(path, sprintf('%2.3f', x))
    if(!file.exists(save.path)) dir.create(save.path)
    
    for(attempt in 1:10) {

      cat(sprintf('running %s = %2.3f, attempt %d\n',
                  param, x, attempt))
    
      model.out <- run.dc.fit(dd.model=dd.model,
                              nc=n.clones,
                              n.iter=n.iter,
                              n.adapt=n.adapt,
                              n.update=n.update,
                              num.cores=1)

      res <- list(dd.model=dd.model, model.out=model.out)

      ## save output
      fn <- sprintf('%s/%d.RData', save.path, attempt)
      save(res, file=fn)
      
      if(max(summary(model.out)[[1]][,'R hat']) < 1.2) break
    }
    NULL
  }
  ## run analysis across all values
  if(cores==1)
    res <- lapply(values, run.analysis)
  if(cores>1)
    res <- mclapply(values, run.analysis,
                    mc.cores=cores,
                    mc.preschedule=FALSE)
  res
}

## ## run and save for various values of parameter 'param'
## vary.param.dclone <- function(save.dir,
##                               param,
##                               values,
##                               n.iter,
##                               n.adapt,
##                               n.update,
##                               nsp,
##                               cores, ...) {

##   ## save directory (and create if the latter does
##   ## not exist)  
##   save.path <-
##     file.path('analyses/occupancy-phylo/saved/simulated',
##               param, save.dir)
##   if(!file.exists(save.path)) dir.create(save.path)

##   ## function to run analysis across reps
##   run.analysis <- function(x) {

##     cat(sprintf('running %s = %2.3f\n', param, x))
    
##     tree <- balanced.tree(nsp)

##     if(param=='lambda')
##       dd.model <- make.data(nsp=nsp,
##                             lambda=x,
##                             tree=tree, ...)
##     if(param=='sigma.psi.beta')
##       dd.model <- make.data(nsp=nsp,
##                             sigma.psi.beta=x,
##                             tree=tree, ...)
    
##     model.out <- run.dc.fit(dd.model=dd.model,
##                             nc=c(20),
##                             n.iter=n.iter,
##                             n.adapt=n.adapt,
##                             n.update=n.update,
##                             num.cores=1)

##     res <- list(dd.model=dd.model, model.out=model.out)

##     ## save output
##     fn <- sprintf('%s/%2.3f.RData', save.path, x)
##     save(res, file=fn)
##     NULL
##   }
##   ## run analysis across all values
##   if(cores==1)
##     res <- lapply(values, run.analysis)
##   if(cores>1)
##     res <- mclapply(values, run.analysis,
##                     mc.cores=cores,
##                     mc.preschedule=FALSE)
##   res
## }

## run and save for various values of parameter 'param'
vary.param.R2jags <- function(save.dir,
                              param,
                              values,
                              nsp,
                              nreps, ...) {

  ## load JAGS model
  model.case <- 'IvesHelmus_p_mII_cts_trait_lambda_rand'
  source(file.path('analyses/occupancy-phylo/src/models',
                   sprintf('%s.R', model.case)))

  ## specify save directory (and creat it if it does not exist)
  save.path <-
    file.path('analyses/occupancy-phylo/saved/simulated',
              param, save.dir)
  if(!file.exists(save.path)) dir.create(save.path)

  ## function to run analysis across reps
  run.analysis <- function(x) {
    for(rep in 1:nreps) {
      save.folder <- sprintf('%s/%2.3f', save.path, x)
      if(!file.exists(save.folder)) dir.create(save.folder)
      cat(sprintf('running %s = %2.3f, rep %d\n', param, x, rep))
      tree <- balanced.tree(nsp)

      if(param=='lambda')
        dd.model <- make.data(nsp=nsp, lambda=x,        tree=tree, ...)
      if(param=='sigma.psi.beta')
        dd.model <- make.data(nsp=nsp, sigma.psi.beta=x, tree=tree, ...)
      
      model.out <- analyse(dd.model, model.case, trait='sim.trait')
      res <- list(dd.model=dd.model, model.out=model.out)

      ## save output
      fn <- sprintf('%d.RData', rep)
      save(res, file=file.path(save.folder, fn))
    }
    NULL
  }
  ## run analysis across all values
  lapply(values, run.analysis)
}
