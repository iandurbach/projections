\name{interpolate_betw_elic_times}
\alias{interpolate_betw_elic_times}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generate interpolated sequences
%%  ~~function to do ... ~~
}
\description{Takes a set of samples at E assessed time points (as output from \code{simulate_at_elic_times} for example) and fills in a user-specified number of values between each of the assessed time points. The `filling in' algorithm differs depending on whether point-in-time values or means were assessed at each of the assessed time points, and on the desired properties of the generated sequences (linear interpolation, Gaussian noise, random walk). 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
interpolate_betw_elic_times(y, nperiods, method, 
maxchange = NULL, alpha = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
%%     ~~Describe \code{y} here~~
an (N x E) matrix or data frame where N is the number of simulations and E is the number of assessed time points. Each column should contain the samples from an assessed probability distribution.
}
  \item{nperiods}{
%%     ~~Describe \code{nperiods} here~~
a vector indicating how many observations each assessed time point is to be converted into. This depends on whether the assessed time points are point-in-time values or means. If method is \code{linear.from.points} or \code{rwalk.from.points} then \code{nperiods} contains the desired number of points to be inserted \emph{between} each assessed time point (i.e. excluding the assessed time-points, which will be inserted in the appropriate place). The total number of generated points will be the sum of the elements of \code{nperiods} plus the number of assessed time-points E. If method is \code{linear.from.means} or \code{gaussian.from.means} then \code{nperiods} contains the desired number of points in each time window. The total number of generated points will be the sum of the elements of \code{nperiods}.
}
  \item{method}{
%%     ~~Describe \code{method} here~~
An interpolation method, one of the following options: \cr\cr
  \code{linear.from.points}: Linear interpolation when the assessed points are point-in-time values. \cr\cr
  \code{rwalk.from.points}: A "random walk" algorithm when the assessed points are point-in-time values. \cr\cr
  \code{linear.from.means}: Linear interpolation when the assessed points are point-in-time values. \cr\cr
  \code{gaussian.from.means}: A `Gaussian' algorithm when the assessed points are means over time.
}
  \item{maxchange}{
%%     ~~Describe \code{maxchange} here~~
If \code{method='gaussian.from.means'} then \code{maxchange} gives the estimated maximum possible change between adjacent time points (e.g. year-to-year) in the constructed sequence. Defaults to \code{NULL}, do not specify if using a method other than \code{gaussian.from.means}. 
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
If \code{method='rwalk.from.points'} then \code{alpha} is a multiplier controlling the size of the random steps between adjacent time points (e.g. year-to-year) in the constructed sequence. Precisely, the step size is \code{alpha * Q} where Q is the largest absolute change observed between adjacent time points in \code{y}, divided by the number of periods to be generated between those time points. Defaults to \code{NULL}, do not specify if using a method other than \code{rwalk.from.points}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
The random walk algorithm should be used with caution as it can give errors or strange results. It is sensitive to the parameter "alpha" (which controls the size of the year-to-year changes). The basic idea is to fix the end points in each time period, and then generate the points in-between as a random walk, subject to some drift and constraints an period-to-period change. If alpha is specified too small, it will not be possible for the algorithm to `reach' the intended destination, giving an error. If this occurs try increasing the value of alpha. Reasonable values seem to be between 2 and 10. Larger values of alpha will give very large period-to-period change. Can also give strange results depicting a "pinching" at some of the end-points. Use with caution!
}
\value{
a (T x N) matrix or data frame, where T is the total number of time points used and N is the number of simulated samples. Each column contains a single projection.
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
# create some data
x_min = c(2.4,2.5,1.5,1.5)
x_mode1 = c(2.6,3,3,3)
x_mode2 = x_mode1
x_max <- c(2.8,3.5,6,6)
correls = matrix(0.3,4,4)
diag(correls) = 1

# construct correlated samples at assessed time points
xe = simulate_at_elic_times(nsims=200,Sigma=correls,x_min,x_mode1,x_mode2,
x_max,mean_rev=c(0,0,1,0))
# interpolate between assessed time points
xt_lfp = interpolate_betw_elic_times(y=xe,nperiods=c(5,14,14),method='linear.from.points')
xt_lfm = interpolate_betw_elic_times(y=xe,nperiods=c(7,15,15),method='linear.from.means')

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
