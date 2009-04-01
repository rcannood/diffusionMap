\name{plot.dmap}
\alias{plot.dmap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot diffusion map }
\description{
  A generic function that plots the diffusion map stored in a dmap object
}
\usage{
plot.dmap(x,\dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ diffusion map object from diffuse() }
  \item{\dots}{ other options used for plot, currently inactive}
}
\author{ Joseph Richards \email{jwrichar@stat.cmu.edu} }
\seealso{ \code{\link{diffuse}} }

\examples{
data(annulus)
dmap = diffuse(dist(annulus)) # compute diffusion map
plot.dmap(dmap)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
\keyword{ nonparametric }% __ONLY ONE__ keyword per line

