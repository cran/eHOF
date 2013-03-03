\encoding{UTF-8}
\name{acre}
\alias{acre}
\alias{acre.env}
\docType{data}

\title{Arable fields in North-Eastern Germany.}

\usage{
data(acre)
}

\description{
Species and Environmental Data from arable fields in North-Eastern Germany. The \code{acre} data frame has 100 selected sites (rows) and 6 species (columns). 
}
\format{
  Two data frames, one data frame called \code{acer} with 100 vegetation plots (rows). Species are named with 7 letter shortletters according to the taxonomic reference list GermanSL. Plant performance is given as presence-absence. The second data frame, named \code{acre.env} consists of plot information about date of survey, survey area, total plant cover and measured soil pH. }

\details{ Soil pH was measured in 1mol KCl solution.}

\references{
Manthey, M. (2003). Vegetations\"okologie der \"Acker und Ackerbrachen Mecklenburg-Vorpommerns. 209 p., Berlin [u.a.]: Cramer.}


\examples{
data(acre)
hist(acre.env$PH_KCL)
}

\keyword{datasets}

