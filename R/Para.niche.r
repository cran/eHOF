"Para.niche" <- function (resp, model, top, opt, pess, central = exp(-0.5), outer = exp(-2), newdata = NULL, ...) 
{
 if (is.null(newdata)) x <- resp$x else x <- newdata
 M <- resp$M
 if (missing(model)) 
    model <- pick.model(resp, gam=FALSE, ...)
 ranx <- diff(resp$range)
 HOFfun <- function(resp, x, y, M) abs(y - predict(resp, new = x, M = M, model = model))

      if (model == "I") {
          outer.low <- resp$range[1]
          outer.high <- resp$range[2]
          central.low <- resp$range[1]
          central.high <- resp$range[2]
          orient <- NA
          }
      if (model == "II") {
          outer <- c(opt, optimize(HOFfun, y = top*eval(outer), resp$range, resp = resp, maximum = FALSE)$minimum)
          outer.low  <- min(outer)
          outer.high <- max(outer)
          central <- c(opt, optimize(HOFfun, y = top*eval(central), resp$range, resp = resp, maximum = FALSE)$minimum)
          central.low  <- min(central)
          central.high <- max(central)
          if(resp$models[[model]]$par[2] < 0) orient = 'increase' else orient = 'decrease'
      }
      if (model == "III") {
          top1 <- resp$models$III$fitted[which(resp$x == resp$range[1])][1]
          top2 <- resp$models$III$fitted[which(resp$x == resp$range[2])][1]
          if (top1 > top2) 
              outer.low <- central.low <- resp$range[1] else 
              outer.high <- central.high <- resp$range[2]

          tmp <- optimize(HOFfun, y = top*eval(outer), resp$range, resp = resp, maximum = FALSE)
          if (top1 > top2) 
              outer.high <- tmp$minimum else 
              outer.low <- tmp$minimum
          tmp <- optimize(HOFfun, y = top*eval(central), resp$range, resp = resp, maximum = FALSE)
          if (top1 > top2) 
              central.high <- tmp$minimum else 
              central.low <- tmp$minimum
          if(resp$models[[model]]$par[2] < 0) orient = 'increase' else orient = 'decrease'
      }
      if (model == "IV") {
          tmp <- optimize(HOFfun, y = top*eval(outer), resp$range, resp = resp, maximum = FALSE)
	  tmp2 <-  opt - diff(c(opt,tmp$minimum))
	  outer.high <- max(tmp$minimum, tmp2)
	  outer.low <- min(tmp$minimum, tmp2)
          tmp <- optimize(HOFfun, y = top*eval(central), resp$range, resp = resp, maximum = FALSE)
	  tmp2 <-  opt - diff(c(opt,tmp$minimum))
	  central.high <- max(tmp$minimum, tmp2)
	  central.low <- min(tmp$minimum, tmp2)
          orient = NA
      }
      if (model == "V") {
 	  x <- seq(min(resp$x)* (-sign(min(resp$x))*10),opt, length.out=10)
          tmp <- optimize(HOFfun, y = top*eval(outer), x, resp = resp, maximum = FALSE)
	  outer.low <- tmp$minimum
	  x <- resp$x[resp$x > opt]
          tmp <- optimize(HOFfun, y = top*eval(outer), c(opt,resp$range[2]), resp = resp, maximum = FALSE)
	  outer.high <- tmp$minimum
          tmp <- optimize(HOFfun, y = top*eval(central), c(resp$range[1], opt), resp = resp, maximum = FALSE)
	  central.low <- tmp$minimum
	  x <- resp$x[resp$x > opt]
          tmp <- optimize(HOFfun, y = top*eval(central), c(opt,resp$range[2]), resp = resp, maximum = FALSE)
	  central.high <- tmp$minimum
          orient <- if((opt - outer.low) < (outer.high - opt)) 'skewed.left' else 'skewed.right'
      }
    if(model %in% c('VI','VII')) {
          central.low <- optimize(HOFfun, y = top['top1']*eval(central), c(resp$range[1],opt['opt1']), resp = resp, maximum = FALSE)$minimum
          central.high <- optimize(HOFfun, y = top['top1']*eval(central), c(opt['opt1'], pess), resp = resp, maximum = FALSE)$minimum
          central2.low <- optimize(HOFfun, y = top['top2']*eval(central), c(pess, opt['opt2']), resp = resp, maximum = FALSE)['minimum']
          central2.high <- optimize(HOFfun, y = top['top2']*eval(central), c(opt['opt2'],resp$range[2]), resp = resp, maximum = FALSE)['minimum']

	  centralBorder <- c(central.low=central.low[[1]], central.high=central.high[[1]], central2.low=central2.low[[1]], central2.high=central2.high[[1]])

          outer.low <- optimize(HOFfun, y = top['top1']*eval(outer), c(resp$range[1],opt['opt1']), resp = resp, maximum = FALSE)$minimum
          outer.high <- optimize(HOFfun, y = top['top1']*eval(outer), c(opt['opt1'], pess), resp = resp, maximum = FALSE)$minimum
          outer2.low <- optimize(HOFfun, y = top['top2']*eval(outer), c(pess, opt['opt2']), resp = resp, maximum = FALSE)['minimum']
          outer2.high <- optimize(HOFfun, y = top['top2']*eval(outer), c(opt['opt2'],resp$range[2]), resp = resp, maximum = FALSE)['minimum']

	  outerBorder <- c(outer.low=outer.low[[1]], outer.high=outer.high[[1]], outer2.low=outer2.low[[1]], outer2.high=outer2.high[[1]])

          orient = NA
	  relfreq.outer <- NA
    } else {
      indx <- outer.low <= resp$x & resp$x <= outer.high
      relfreq.outer <- sum(resp$y[indx] > 0) / sum(indx, na.rm=TRUE)
      outerBorder <- c(outer.low=outer.low, outer.high=outer.high)
      centralBorder <- c(central.low=central.low, central.high=central.high)
    }
  border <- list(centralBorder=centralBorder, outerBorder=outerBorder, orient=orient, relfreq.outer=relfreq.outer)
  border
}

