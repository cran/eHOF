Para.opt <- function (
		resp, 
		model=NULL, 
		punctual = FALSE, 
		newdata = NULL, 
		...) {
  if(is.null(model)) model <- pick.model(resp, gam=FALSE, ...)
  M <- resp$M
#  if (missing(newdata)) x <- seq(resp$range[1], resp$range[2], length.out=10000) else x <- newdata
  x <-  if (missing(newdata)) seq(min(resp$range),max(resp$range),length.out=10000)	else newdata

  HOFfun <- function(x, y, resp) abs(y - predict(resp, newdata = x, M = M, model = model))

  if (model == "I") {
      opt <- NA  # resp$range
      if (missing(newdata)) top <- fitted(resp, 'I')[1] else  top <- predict(resp, 'I', newdata)[1]
      mini <- top
      pess <- NA
  } 

  if (model == "II") {
      tmp <- optimize(HOFfun, resp$range, y=0, resp = resp, maximum = TRUE)
      opt <- as.numeric(tmp$maximum)
      top <- as.numeric(tmp$objective)
      tmp <- optimize(HOFfun, resp$range, y=0, resp = resp, maximum = FALSE)
      pess <- as.numeric(tmp$minimum)
      mini <- as.numeric(tmp$objective)
  }

  if (model == "III") {
        top <- optimize(HOFfun, interval=resp$range, y = 0, resp = resp, maximum = TRUE)$objective
	## optimize does not always find an adequate solution
        opt <- optimize(HOFfun, resp$range, y = top * 9/10, resp = resp, maximum = FALSE)$minimum
        if(resp$models$III$par['a'] < 0) 
		opt <- if(punctual) opt - (opt - min(resp$range))/2	else
		c(opt.min = resp$range[1], opt.max = opt)
	if(resp$models$III$par['a'] >= 0) 
		opt <- if(punctual) opt + (max(resp$range) - opt)/2	else
		c(opt.min = opt, opt.max = resp$range[2])
        pess <- optimize(HOFfun, y = 0, resp$range, resp = resp, maximum = FALSE)
        mini <- pess$objective
        if (predict(resp, newdata = min(resp$range), M = M, model = model) > predict(resp, newdata = max(resp$range), M = M, model = model))  pess <- c(pess$minimum, max(resp$range)) else pess <- c(min(resp$range), pess$minimum)
    }

   if (model == "IV") {
      p <- coef(resp, model)
      ranx <- diff(resp$range)
      minx <- resp$range[1]           
      opt <- (p[3] - p[1])/p[2]/2
      opt <- as.numeric(ranx * opt + minx)
      top <- as.numeric(predict(resp, newdata = opt, M = M, model = "IV"))
      tmp <- optimize(HOFfun, resp$range, y=0, resp = resp, maximum = FALSE)
      mini <- as.numeric(tmp$objective)
      pess <- as.numeric(tmp$minimum)
  }
  if (model == "V") {
      p <- coef(resp, model)
      if (p[2] * p[4] >= 0) {
          tmp <- optimize(HOFfun, resp$range, y=0, resp = resp, maximum = TRUE)
          opt <- as.numeric(tmp$maximum)
          top <- as.numeric(tmp$objective)
          pess <- as.numeric(tmp$minimum)
          if (top < 16 * .Machine$double.eps) {
              tmp <- seq(resp$range[1], resp$range[2], len = 31)
              ytmp <- predict(resp, newdata = tmp, model = "V")
              tmp <- tmp[ytmp > 16 * .Machine$double.eps]
              tmp <- optimize(HOFfun, range(tmp), resp = resp, y=0, maximum = TRUE)
              opt <- as.numeric(tmp$maximum)
              top <- as.numeric(tmp$objective)
              if (tmp$obj <= 0) 
                opt <- top <- NA
          }
      }
      else opt <- top <- NA
      mini <- as.numeric(optimize(HOFfun, resp$range, y=0, resp = resp, maximum = FALSE)[["objective"]])
      pess <- as.numeric(tmp$minimum)
  }
  
  if (model == "VI") {
      max1 <- optimize(HOFfun, resp$range, resp = resp, y=0, maximum = TRUE)
  	  min2 <- optimize(HOFfun, lower=max1$maximum, upper=resp$range[2], resp = resp, y=0, maximum = FALSE)
#  if(round(max1$objective,2) < round(max(predict(resp)),2))  
#  	  min <- optimize(HOFfun, lower=max1$maximum, upper=resp$range[2], resp = resp, y=max1$objective, maximum = FALSE)
#      max2 <- optimize(HOFfun, c(min1$maximum, resp$range[[2]]), resp = resp, y=min1$objective, maximum = FALSE)
#       mini <- min(predict(resp, newdata=x))
      mini <- min2$objective
      max2 <- optimize(HOFfun, c(min2[[1]], resp$range[[2]]), resp = resp, y=mini, maximum = TRUE)
      top <- c(top1=max1$objective, top2 =max2$objective)
      pess <- min2$minimum
      opt <- c(opt1 = max1$maximum, opt2 = max2$maximum)
      new <- seq(resp$range[1], pess, length.out = 5000)
      pm <- predict.HOF(resp, newdata = new, model = model)
      expect1 <- sum(pm * new)/sum(pm)
      new <- seq(pess, resp$range[2], length.out = 5000)
      pm <- predict.HOF(resp, newdata = new, model = model)
      expect2 <- sum(pm * new)/sum(pm) 
      expect <- c(expect1, expect2)
   }

  
  if (model == 'VII') {
      max1 <- optimize(HOFfun, resp$range, resp = resp, y=0, maximum = TRUE)
  	  max2 <- optimize(HOFfun, lower=max1$maximum, upper=resp$range[2], resp = resp, y=max1$objective, maximum = FALSE)
      tmp <- optimize(HOFfun, c(max1$maximum, max2[[1]]), resp = resp, y=0, maximum = FALSE)
      mini <- as.numeric(tmp$objective)
      top <- c(top1=min(max1$objective, max2$objective), top2 =max(max1$objective, max2$objective))
      pess <- as.numeric(tmp$minimum)
      opt <- c(opt1 = min(max1$maximum, max2$minimum), opt2 = max(max1$maximum, max2$minimum))
      pess <- as.numeric(tmp$minimum)
      new <- seq(resp$range[1], pess, length.out = 5000)
      pm <- predict.HOF(resp, newdata = new, model = model)
      expect1 <- sum(pm * new)/sum(pm)
      new <- seq(pess, resp$range[2], length.out = 5000)
      pm <- predict.HOF(resp, newdata = new, model = model)
      expect2 <- sum(pm * new)/sum(pm) 
      expect <- c(expect1, expect2)
  }
  if(model %in% c('I', 'II', 'III', 'IV', 'V')) {      
# 		new <- seq(resp$range[1], resp$range[2], length.out = 10000)
		pm <- predict.HOF(resp, newdata = x, model = model)
		expect <- sum(pm * x)/sum(pm)
	}

  list(opt = opt, top = top, pess = pess, mini = mini, expect = expect)
}
