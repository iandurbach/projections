\name{interpolate_rwalk_frompoints}
\alias{interpolate_rwalk_frompoints}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Interpolation between time-points using a random walk
%%  ~~function to do ... ~~
}
\description{Interpolates using a "random walk" algorithm when the assessed points are point-in-time values. This function is called by \code{interpolate_betw_elic_times} if \code{method='rwalk.from.points'}.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
interpolate_rwalk_frompoints(z, nperiods, random_step_all)
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
  \item{random_step_all}{
%%     ~~Describe \code{random_step_all} here~~
a reasonable step size for the random walk, as computed by \code{compute_random_step} with parameter \code{alpha}. See \code{compute_random_step} and \code{interpolate_betw_elic_times} for details.
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
interpolate_rwalk_frompoints(z, nperiods,random_step_all=c(1,1,1))}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
