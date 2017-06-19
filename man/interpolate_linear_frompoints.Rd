\name{interpolate_linear_frompoints}
\alias{interpolate_linear_frompoints}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Linear interpolation between time-points
%%  ~~function to do ... ~~
}
\description{Linear interpolation when the assessed points are point-in-time values. This function is called by \code{interpolate_betw_elic_times} if \code{method='linear.from.points'}.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
interpolate_linear_frompoints(z, nperiods)
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
interpolate_linear_frompoints(z, nperiods)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
