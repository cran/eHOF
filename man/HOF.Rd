\encoding{UTF-8}
\name{HOF}
\alias{HOF}
\alias{HOF.default}
\alias{HOF.data.frame}
\alias{fitted.HOF}
\alias{pick.model}
\alias{pick.model.HOF}
\alias{pick.model.HOF.list}
\alias{predict.HOF}
\alias{print.HOF}
\alias{print.HOF.list}

\title{Hierarchical logistic regression}

\description{
 Fit seven hierarchical logistic regression models and select the most appropriate model 
 by information criteria and a bootstrap approach to guarantee model stability. 
 The first five shapes are known as Huisman-Olff-Fresco (HOF) models in ecology (Huisman et al. 1993). 
 Additionally the package provides two bimodal shapes.
}

\usage{
	\method{HOF}{default}(occ, grad, M = max(occ), y.name, family = binomial, lim = 100, 
  bootstrap = 100,   test = c('AICc', 'BIC', 'AIC','Dev'), ...)
	\method{HOF}{data.frame}(veg, grad, M, freq.limit = 10, ...)
	\method{pick.model}{HOF}(object, level = 0.95, test = c("AICc", "BIC", "AIC", "Dev"), 
  modeltypes, penal = 'df', gam = FALSE, selectMethod = c('bootselect','IC.weight',
  'pick.model'), ...)
	\method{fitted}{HOF}(object, model, ...)
	\method{predict}{HOF}(object, model, newdata, ...)
}

\arguments{
  \item{occ}{species occurrences (=response) vector.}
  \item{veg}{vegetation data frame. Either as matrix with species in columns and plots in rows, or in Turboveg format.}
  \item{grad}{gradient data vector.}
  \item{M}{maximum attainable value in the HOF model, similar to binomial denominator.}
  \item{y.name}{name of the species.}
  \item{family}{error distribution. Alternatives are \code{binomial}, \code{poisson} and \code{gaussian}.}
  \item{lim}{limiting criterion for minimization function, see details.}
  \item{bootstrap}{number of bootstrap samplings to check model robustness, see details.}
  \item{freq.limit}{lowest frequency of species analysed.}
  \item{object}{HOF model object, results from function HOF. Either for a single species or a list for several species.}
  \item{level}{probability for model selection (1-P) for F or Chisq tests.}
  \item{test}{information criterion for model selection. Alternatives are \code{"AICc"}, \code{"BIC"}, \code{"AIC"} or just \code{"Dev"iance}. }
  \item{modeltypes}{vector of model types, when only a subset of the six implemented modeltypes should be considered.}
  \item{penal}{penalty term for model types, default is the number of model parameter.}
  \item{gam}{calculate AIC of GAM model and compare it with HOF models.}
  \item{selectMethod}{choose the model type selection method in case of divergence of the most frequent bootstrapped model 
  		or the one with highest information criterion weight from the primary chosen type.}
  \item{newdata}{vector of gradient values for prediction.}
  \item{model}{manually specify the HOF model type to be used, otherwise it is selected through function \code{\link{pick.model}}.}
  \item{\dots}{further arguments passed to or from other methods}
}

\details{
  The minimization function changed compared to package gravy (Oksanen 2002) from nlm to nlminb to be able to set a limit 
  for estimated parameters (default= -100 to 100). The old models III and V have been often too sharp, 
  \code{lim=Inf} will produce results similar to gravy.
  Function \code{\link{pick.model}} finds the most adequate modeltype according to the chosen Information Criterion (AICc is default).
  Function \code{fitted} returns the fitted values for the used \code{grad}ient, and \code{predict} for any values in \code{newdata}.
  
  To improve and check model stability a bootstrapping mechanism is implemented in function HOF. If the initially chosen model type is different from the most frequent one, the latter will be chosen by default.
 }

\value{
  \code{HOF.fit} returns an object of class \code{"HOF"} which contains the parameters of the call, 
  the fitting results for every model type and a vector of chosen model types from bootstrapping.
 }

\seealso{
  \code{\link{plot.HOF}} provides advanced plotting schemes for HOF models. 
  \code{\link{Para}} derives model parameters like optimum, niche (width), slope etc.
 }

\references{
Jansen, F. & Oksanen, J. (2013) What shape are species responses along ecological gradients? 
- Huisman-Olf-Fresco models revisited. Journal of Vegetation Science, DOI: 10.1111/jvs.12050

Oksanen, J. & Minchin, P.R. (2002). Continuum theory revisited: what shape are species responses along ecological gradients? 
\emph{Ecological  Modelling} 157, 119-129.

Huisman, J., Olff, H. & Fresco, L.F.M. (1993). A hierarchical set of models for species response analysis. 
\emph{Journal of Vegetation Science} 4, 37-46.
 }

\author{Florian Jansen, Jari Oksanen}

\examples{
	data(acre)
	sel <- c('MATRREC', 'RUMEACT', 'SILENOC', 'APHAARV', 'MYOSARV', 'DESUSOP', 'ARTE#VU')
	mo <- HOF(acre[,match(sel, names(acre))], acre.env$PH_KCL, M=1, bootstrap=NULL)
	mo
}

\keyword{ models }
\keyword{ regression }
\keyword{ nonlinear }
