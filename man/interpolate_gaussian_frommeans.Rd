\name{interpolate_gaussian_frommeans}
\alias{interpolate_gaussian_frommeans}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Interpolation for means with Gaussian noise
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Interpolates using a `Gaussian' algorithm when the assessed points are means over time. This function is called by \code{interpolate_betw_elic_times} if \code{method='gaussian.from.means'}.
}
\usage{
interpolate_gaussian_frommeans(z, nperiods, maxchange)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z}{
%%     ~~Describe \code{z} here~~
a single row of \code{y}
}
  \item{nperiods}{
%%     ~~Describe \code{nperiods} here~~
a vector indicating how many observations each assessed time point is to be converted into. See \code{interpolate_betw_elic_times}. 
}
  \item{maxchange}{
%%     ~~Describe \code{maxchange} here~~
gives the estimated maximum possible change between adjacent time points (e.g. year-to-year) in the constructed sequence. See \code{interpolate_betw_elic_times}. 
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
a vector containing the full generated time series (both interpolated and any assessed points) 
}
\references{
%% ~put references to the literature/web site here ~
}
\author{Ian Durbach
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
z = c(3,6,7,18)
nperiods = c(5,10,20)
interpolate_gaussian_frommeans(z, nperiods,maxchange=3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
