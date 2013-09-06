"print.HOF.list" <-  function (
		x, 
		test = 'AICc',  
		selectMethod = 'bootselect', 
		...) {
	if(is.null(x[[1]]$bootstrapmodels)) {
		cat('Warning: No bootstrap results in object. Using selectMethod "raw" instead.\n')
		selectMethod <- 'raw'
	}
    cat("Deviances:\n")
    printCoefmat(sapply(x, deviance), na.print="", ...)
    cat(paste("\nSuggested best models (",test, ", ", selectMethod, "):", sep=''))
    tmp <- sapply(x, pick.model, test=test, quiet = TRUE, selectMethod = selectMethod, ...)
    names(tmp) <- names(x)
    cat('\n')
    print(noquote(tmp))
#    if(!is.null(x[[1]]$bootstrapmodels)) {
#      bootpick <- sapply(x, function(x) x$bootstrapmodels)
#      cat('\nPercentage of equal model types during bootstrapping\n')
#      perc <- vector('integer', length=ncol(bootpick))
#      names(perc) <- names(x)
#      for(i in 1:ncol(bootpick)) perc[i] <- round(sum(bootpick[,i]==tmp[i])/nrow(bootpick)*100,0)
#      print(perc)
#      }
    invisible(x)
}
