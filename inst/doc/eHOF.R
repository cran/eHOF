## ----prep, echo=FALSE, results='hide'-----------------------------------------
suppressPackageStartupMessages(library(vegdata))
tmp <- tempdir()
options(tv_home = tmp)
dir.create(file.path(tmp, 'Species'))
dir.create(file.path(tmp, 'Popup'))
dir.create(file.path(tmp, 'Data'))
file.copy(from = file.path(path.package("vegdata"), 'tvdata', 'Popup'), to = tmp, recursive = TRUE)
file.copy(from = file.path(path.package("vegdata"), 'tvdata', 'Species'), to = tmp, recursive = TRUE)
file.copy(from = file.path(path.package("vegdata"), 'tvdata', 'Data'), to = tmp, recursive = TRUE)

## ----load, messages=FALSE-----------------------------------------------------
library(eHOF)

## ----4site.echo, results='hide', message=FALSE--------------------------------
library(vegdata)
db <- 'elbaue'
site <- tv.site(db)
veg <- tv.veg(db, tax=FALSE, spcnames='numbers')
obs <- tv.obs(db)
# taxa <- tax(unique(obs$TaxonUsageID), verbose=TRUE)
# write.csv2(taxa, file='taxonnames.csv')
taxa <- read.csv('taxonnames.csv')
names(veg) <- sub('.0', '', names(veg), fixed=TRUE)
for(i in 1:ncol(veg)) names(veg)[i] <- as.character(taxa$LETTERCODE[match(as.numeric(names(veg)[i]), taxa$TaxonUsageID)])

## ----5veg.2, results='hide'---------------------------------------------------
veg.sqrt <- tv.veg(db, cover.transform='sqrt', tax=FALSE, spcnames='numbers')
names(veg.sqrt) <- sub('.0', '', names(veg.sqrt), fixed=TRUE)
names(veg.sqrt) <- taxa$LETTERCODE[match(names(veg.sqrt), taxa$TaxonUsageID)]

## ----5veg.3, results='hide'---------------------------------------------------
veg.pa <- tv.veg(db, cover.transform='pa', tax=FALSE, spcnames='numbers')
names(veg.pa) <- sub('.0', '', names(veg.pa), fixed=TRUE)
names(veg.pa) <- taxa$LETTERCODE[match(names(veg.pa), taxa$TaxonUsageID)]

## ----modeltypes, warning=FALSE, results='hide'--------------------------------
data(acre)
sel <- c('ELYMREP', 'VEROPES', 'CONSREG', 'DESUSOP', 'VEROARV', 'ARTE#VU', 'ACHIMIL')
mo <- HOF(acre[match(sel, names(acre))], acre.env$PH_KCL, M=1, bootstrap=NULL)
par(mar=c(4,4,1,1)+.1)
autolayout(7)
par(mar=c(4,4,1,1)+.1)
for(i in 1:7) plot(mo[[i]], model = eHOF.modelnames[i], marginal ='n')

## ----percentageVEG, results='hide', warning=FALSE-----------------------------
mods <- HOF(veg, site$MGL, M=100, family = poisson, bootstrap = NULL)

## ----printPerc----------------------------------------------------------------
mods

## ----7sqrt, results='hide', warning=FALSE-------------------------------------
mods.sq <- HOF(veg.sqrt, site$MGL, M=10, family= poisson, freq.limit=10, bootstrap=NULL)
plot(mods.sq)

## ----7pres-abs, results='hide', warning=FALSE---------------------------------
mods.pa <- HOF(veg.pa, site$MGL, M=1, bootstrap=NULL)
plot(mods.pa)

## ----paraplot, warning=FALSE, out.width='.5\\textwidth'-----------------------
plot(mods.pa[['RANCREP']], para=TRUE, onlybest=FALSE)

## ----sqrt-ELYMREP, results='hide', warning=FALSE, out.width='.5\\textwidth'----
mod.ELYM.sqrt <- HOF(veg.sqrt$ELYMREP, site$MGL, M=10, family=poisson, bootstrap = 10)
plot(mod.ELYM.sqrt, marginal='point', para=TRUE, onlybest=FALSE, newdata=seq(min(mod.ELYM.sqrt$range), max(mod.ELYM.sqrt$range),length.out=10000) )
Para(mod.ELYM.sqrt)

