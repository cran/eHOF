\encoding{UTF-8}
\name{autolayout}
\alias{autolayout}

\title{Automatically Specifying Complex Plot Arrangements}

\description{
 Automatic version of function \code{\link{layout}} for given number of subplots. 
 Defined matrix schemes for \code{\link{layout}} for up to 32 subfigures.
}

\usage{
autolayout(N, byrow = TRUE, widths, heights, respect)
}

\arguments{
  \item{N}{number of subfigures}
  \item{byrow}{see \code{\link{matrix}}}
  \item{widths}{a vector of values for the widths of columns on the device. Relative widths are specified with numeric values. Absolute widths (in centimetres) are specified with the lcm() function (see examples).}
  \item{heights}{a vector of values for the heights of rows on the device. Relative and absolute heights can be specified, see widths above.}
  \item{respect}{either a logical value or a matrix object. If the latter, then it must have the same dimensions as mat and each value in the matrix must be either 0 or 1.}
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
