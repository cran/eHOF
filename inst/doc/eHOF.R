
## ----1,eval=FALSE--------------------------------------------------------
## vignette("eHOF")


## ----prep, echo=FALSE----------------------------------------------------
options(warn=1)


## ----load, messages=FALSE------------------------------------------------
library(eHOF)


## ----4site.echo, results='hide', message=FALSE---------------------------
library(vegdata)
db <- 'elbaue'
site <- tv.site(db)
veg <- tv.veg(db, tax=FALSE, spcnames='numbers')
obs <- tv.obs(db)
#taxa <- tax(unique(obs$TaxonUsageID), verbose=TRUE)
#write.csv2(taxa, file='taxa.csv')
taxa <- read.delim('taxonnames.csv')
names(veg) <- sub('.0', '', names(veg), fixed=TRUE)
for(i in 1:ncol(veg)) names(veg)[i] <- as.character(taxa$LETTERCODE[match(as.numeric(names(veg)[i]), taxa$TaxonUsageID)])


## ----5veg.2, results='hide'----------------------------------------------
veg.sqrt <- tv.veg(db, cover.transform='sqrt', tax=FALSE, spcnames='numbers')
names(veg.sqrt) <- sub('.0', '', names(veg.sqrt), fixed=TRUE)
names(veg.sqrt) <- taxa$LETTERCODE[match(names(veg.sqrt), taxa$TaxonUsageID)]


## ----5veg.3, results='hide'----------------------------------------------
veg.pa <- tv.veg(db, cover.transform='pa', tax=FALSE, spcnames='numbers')
names(veg.pa) <- sub('.0', '', names(veg.pa), fixed=TRUE)
names(veg.pa) <- taxa$LETTERCODE[match(names(veg.pa), taxa$TaxonUsageID)]


## ----modeltypes, warning=FALSE, results='hide'---------------------------
data(acre)
sel <- c('ELYMREP', 'VEROPES', 'CONSREG', 'DESUSOP', 'VEROARV', 'ARTE#VU', 'ACHIMIL')
mo <- HOF(acre[match(sel, names(acre))], acre.env$PH_KCL, M=1, bootstrap=NULL)
par(mar=c(4,4,1,1)+.1)
autolayout(7)
par(mar=c(4,4,1,1)+.1)
for(i in 1:7) plot(mo[[i]], model = eHOF.modelnames[i], marginal ='n')


## ----percentageVEG, results='hide', warning=FALSE------------------------
mods <- HOF(veg, site$MGL, M=100, family = poisson, bootstrap = NULL)


## ----printPerc-----------------------------------------------------------
mods


## ----7sqrt, results='hide', warning=FALSE--------------------------------
mods.sq <- HOF(veg.sqrt, site$MGL, M=10, family= poisson, freq.limit=10, bootstrap=NULL)
plot(mods.sq)


## ----7pres-abs, results='hide', warning=FALSE----------------------------
mods.pa <- HOF(veg.pa, site$MGL, M=1, bootstrap=NULL)
plot(mods.pa)


## ----paraplot, warning=FALSE, out.width='.5\\textwidth'------------------
plot(mods.pa[['RANCREP']], para=TRUE, onlybest=FALSE)


## ----sqrt-ELYMREP, results='hide', warning=FALSE, out.width='.5\\textwidth'----
mod.ELYM.sqrt <- HOF(veg.sqrt$ELYMREP, site$MGL, M=10, family=poisson, bootstrap = 10)
plot(mod.ELYM.sqrt, marginal='point', para=TRUE, onlybest=FALSE, newdata=seq(min(mod.ELYM.sqrt$range), max(mod.ELYM.sqrt$range),length.out=10000) )


## ----11, results='hide', warning=FALSE, fig.height=4, fig.width=4--------
mods.pa <- HOF(veg.pa, site$MGL, M=1, freq.limit=7, family=poisson, bootstrap=NULL)
p <- Para(mods.pa)
p4 <- vector('list', length(mods.pa))
for(i in 1:length(mods.pa)) p4[[i]] <- Para(mods.pa[[i]], modeltypes='IV')
#p4
ind <- !sapply(p, function(x) x$model) %in% c('III','VI', 'VII')
p.opt <- sapply(p, function(x) x$opt)
p4.opt <- sapply(p4, function(x) x$opt)
plot(p.opt[ind], p4.opt[ind], xlab='mean groundwater level (MGL)', ylab='MGL only with model IV')
for(i in 1:length(p.opt[!ind]))
    lines(unlist(p.opt[!ind][i]), rep(p4.opt[!ind][i], 2), col='darkgreen')
points(c(-270), p4.opt[is.na(p.opt)], pch='+')
lines(c(-600,0),c(-600,0), lty=2)


