\encoding{UTF-8}
\name{autolayout}
\alias{autolayout}

\title{Automatically Specifying Complex Plot Arrangements}

\description{
 Automatic version of function \code{\link{layout}} for given number of subplots. 
 Defined matrix schemes for \code{\link{layout}} for up to 32 subfigures.
}

\usage{
autolayout(N, byrow = TRUE, ...)
}

\arguments{
  \item{N}{number of subfigures}
  \item{byrow}{see \code{\link{matrix}}}
  \item{...}{additional arguments for \code{\link{layout}} or \code{\link{matrix}}.}
}

\value{Returns the number of figures, N, see \code{\link{layout}}.}

\author{
 Florian Jansen
}

\seealso{
 \code{\link{layout}}; \code{\link{matrix}}
}

\examples{
for(i in sample(1:32, 1)) {
  autolayout(i, byrow=TRUE)
  layout.show(i) 
}
}

\keyword{misc}
