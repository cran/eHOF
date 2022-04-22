\encoding{UTF-8}
\name{betaresp}
\alias{betaresp}
\alias{print.betaresp}
\alias{plot.betaresp}
\alias{betaresponse}

\title{Beta response models}

\description{
 Generating beta response models. Implemented to compare the hierarchical logistic regression models 
 with beta response models.
}

\usage{
 betaresp(x, p1, p2, alpha, gamma, hi)
}

\arguments{
 \item{x}{Gradient values}
 \item{p1, p2}{response endpoints}
 \item{alpha, gamma}{shape parameters of the response}
 \item{hi}{maximum height of response}
}


\references{
Jansen, F. & Oksanen, J. (2013) How to model species responses along ecological gradients? - 
 Huisman-Olff-Fresco models revisited, Journal of Vegetation Science doi:10.1111/jvs.12050

Oksanen, J. & Minchin, P.R. (2002). Continuum theory revisited: what shape are species responses along ecological gradients? \emph{Ecological Modelling} 157, 119-129.
 }
 
\author{Jari Oksanen}

\examples{
x <- seq(0,3, len=201)
m <- betaresp(x, 2, 0.8, 0.1, 0.5, 0.8)
plot(m)
}

\keyword{ models }
\keyword{ regression }
\keyword{ nonlinear }

