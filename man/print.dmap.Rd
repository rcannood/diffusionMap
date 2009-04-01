\name{print.dmap}
\alias{print.dmap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Extract diffusion map eigen-multipliers and epsilon}
\description{
  A generic function that prints the eigen-multipliers of a dmap object
  and the epsilon used in its construction.
}
\usage{
print.dmap(x,\dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ diffusion map object from diffuse() }
  \item{\dots}{ further arguments}
}
\author{ Joseph Richards \email{jwrichar@stat.cmu.edu} }
\seealso{ \code{\link{diffuse}} }

\examples{
data(annulus)
dmap = diffuse(dist(annulus)) # compute diffusion map
print.dmap(dmap)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
\keyword{ nonparametric }% __ONLY ONE__ keyword per line

