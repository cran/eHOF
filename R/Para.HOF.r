"Para.HOF" <- function (
		resp, 
		model, 
		newdata = NULL, 
		...)
{
    if (missing(model)) model <- pick.model(resp, gam=FALSE, ...)
    if (is.null(newdata)) x <- scale01(seq(resp$range[1], resp$range[2], length.out=1000)) else x <- scale01(newdata, ...)
    M <- resp$M
    opt <- Para.opt(resp, model=model, ...)
    border <- Para.niche(resp, model=model, top=opt$top, opt=opt$opt, mini=opt$mini, pess=opt$pess, ...)
    slope <- Para.deriv(resp, x, p=resp$models[[model]]$par, model)
    infl <- Para.deriv(resp, x, p=resp$models[[model]]$par, model, optima=opt$opt, pessima=opt$pess, type='inflection')
    max.sl <- max(abs(slope))
    out <- list(species = resp$y.name, abund.sum = sum(resp$y/M), range = resp$range, model = model, para = resp$models[[model]]$par, M = M, mini = opt$mini, pess= opt$pess, top = opt$top, opt = opt$opt, max.slope=max.sl,  inflection=infl, expect = opt$expect) 
    out$centralBorder <- border$centralBorder
    out$outerBorder <- border$outerBorder
    out$raw.mean <-  mean(resp$x[resp$y>0])
    class(out) <- c("Para.HOF")
    out
}


