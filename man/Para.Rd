\encoding{UTF-8}
\name{Para}
\alias{Para}
\alias{Para.HOF}
\alias{Para.HOF.list}

\title{Curve parameters of eHOF models}

\description{
  Derive common shape parameters from the different model types. 
  Calculate a set of parameters (see values below) from eHOF models.
}

\usage{
\method{Para}{HOF}(resp, model, newdata = NULL, ...)
\method{Para}{HOF.list}(resp, ...)
}

\arguments{
 \item{resp}{response model results, see \code{\link{HOF}}.}
 \item{model}{response model type. If not specified, the output of \code{\link{pick.model}} will be used.}
 \item{newdata}{new gradient vector.}
 \item{\dots}{further arguments passed to or from other methods, e.g. for \code{\link{pick.model}}.}
}

\value{
	\item{species}{Name or abbreviat of the species.}
	\item{abund sum}{Abundance sum, i.e. sum of all response values divided by M.}
	\item{range}{Range of x values.}
	\item{model}{Model type, if not specified the result of \code{\link{pick.model}}.}
	\item{para}{Model parameters (a to d).}
	\item{M}{Maximum response value, specified in the HOF function call.}
	\item{mini}{Location of the minimum, i.e. the gradient value, where the response is lowest, for model VI and VII the lowest response between the two optima.}
	\item{pess}{Lowest estimated response value.}
	\item{top}{Highest estimated response value(s).}
	\item{opt}{Location of the optimum, i.e. the gradient value, where the species response is highest. NA for model I and an optimum interval for model type III.}
	\item{expect}{ Expectancy value, i.e. average x value under the model curve).}
	\item{max slope}{Highest slope, i.e. maximum of the first derivation of the curve.}
	\item{centralBorder}{Following Heegard, the central borders are calculated as the gradient values, where the response reaches \code{"exp(-1/2)"} of the top.}
	\item{outerBorder}{Following Heegard, the outer borders of the species niche are calculated as the gradient values, where the response reaches \code{exp(-2)} of the top.}
	\item{raw mean}{Average of measured x values.}
}

\details{
For models VI and VII \code{Para} will give you two expectancy values for the ranges left and right of the pessimum between the model optima. If you want to have the overall expectancy value, use: 
\code{
  gradient <- seq(min(Para(resp)$range), max(Para(resp)$range), length.out=10000)
  weighted.mean(gradient, predict(resp, newdata=gradient))
}
}

\references{
Heegard, E. (2002) The outer border and central border for species-environmental relationships estimated by non-parametric generalised additive models. Ecological Modelling, 157, 131-139.
Damgaard, C. (2006) Modelling ecological presence-absence data along an environmental gradient: threshold levels of the environment. Environ Ecol Stat 13:229-236.
}

\author{Florian Jansen}

\keyword{models}

