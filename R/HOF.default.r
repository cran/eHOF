HOF.default <- function(
		occ, 
		grad, 
		M = max(occ), 
		y.name, 
		family=binomial, 
		lim=100, 
		bootstrap=100, 
		test = c('AICc', 'BIC', 'AIC','Dev'), 
		...)  {
  if(any(c('data.frame', 'matrix','list') %in% class(occ))) stop('Occurrence data for HOF.default must be a vector.')
  x.name <- deparse(substitute(grad))
  if (missing(y.name)) y.name <- deparse(substitute(occ))
  if(any(is.na(occ))) stop('NA in occurrence vector is not allowed!')
  if(!is.numeric(grad)) print('Gradient must be a numeric vector')
#  if(is.null(bootstrap)) print('If you want to ensure model stability use bootstraps!') else 
  if(!is.null(bootstrap)) if(bootstrap == 0) stop('If you do not want to bootstrap your data use "bootstrap=NULL"!')
  out <- HOF.model(occ, grad, M, y.name, x.name, family=family, lim=lim,...)
  
  IC.weights <- function(x, test = 'AICc') {
	  penal <- sapply(x$models, function(x) length(x$par))
	  dev <- deviance(x)
	  ll <- logLik(x)
	  AICc <- -2 * ll + 2 * penal + 2 * penal *(penal + 1)/(x$nobs - penal - 1) 
	  d.AICc <- AICc - min(AICc, na.rm=TRUE)
	  AICc.W <- round(exp(-0.5*AICc)/ sum(exp(-0.5*AICc), na.rm=TRUE),4)
	  return(AICc.W)
  }
  
  if(!is.null(bootstrap)) {
    test <- match.arg(test)
    modeltypes <- character(length=bootstrap)
    mods <- vector('list', length=bootstrap)
	weights <-  matrix(nrow=bootstrap, ncol=7); colnames(weights) <- eHOF.modelnames
    pb <- txtProgressBar (min = 0, max = bootstrap, char = '.',  width = 45, style = 3)
    for(i in 1:bootstrap) {
      take <- sample(length(grad), replace=TRUE)
      mods[[i]] <- HOF.model(occ[take], grad[take], M=M, y.name, x.name, family=family, lim=lim,...)
      modeltypes[i] <- pick.model(mods[[i]], quiet=TRUE, ...)
	  weights[i,] <- IC.weights(mods[[i]])
      setTxtProgressBar(pb, bootstrap - (bootstrap - i))
    }
    close (pb) ## Close progress bar
    out$call <- match.call()
    out$bootstraptest <- test
    out$bootstrapmodels <- modeltypes
	out$ICweights <- weights
  }
  out
}

