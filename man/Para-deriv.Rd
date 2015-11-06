\encoding{UTF-8}
\name{Para_deriv}
\alias{Para_deriv}

\title{Curve derivatives of eHOF models}

\description{
  Slope and inflection points for the eHOF model types.
}

\usage{
Para_deriv(resp, newdata = NULL, model = NULL, p, type = c("slope", "inflection"), ...)
}

\arguments{
 \item{resp}{response model results, see \code{\link{HOF}}.}
 \item{newdata}{new gradient vector, if NULL a vector with 10000 values within gradient range will be used}
 \item{model}{response model type. If not specified, the output of \code{\link{pick.model}} will be used.}
 \item{p}{Vector of model parameters. If not specified, resp$models[[model]]$par will be used.}
 \item{type}{Output type: slope or inflection points of model curves}
 \item{\dots}{further arguments passed to or from other methods, e.g. for \code{\link{pick.model}}.}
}

\value{
  For slope a vector of slopes for all newdata points. Inflection points are the extremes of the second derivative.
}

\references{
Peper, J., Jansen, F., Pietzsch, D., Manthey, M. 2011. Patterns of plant species turnover along grazing gradients. Journal of Vegetation Science 22: 457–466.
}

\author{Florian Jansen}

\keyword{models}
