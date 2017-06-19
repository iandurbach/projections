\name{compute_random_step}
\alias{compute_random_step}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Computes step size for random walk
%%  ~~function to do ... ~~
}
\description{Called by \code{interpolate_rwalk_frompoints} to compute a reasonable step size for the random walk. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
compute_random_step(y, nperiods, alpha)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
%%     ~~Describe \code{y} here~~
an (N x E) matrix or data frame where N is the number of simulations and E is the number of assessed time points. Each column should contain the samples from an assessed probability distribution.
}
  \item{nperiods}{
%%     ~~Describe \code{nperiods} here~~
a vector indicating the desired number of points to be inserted \emph{between} each assessed time point (i.e. excluding the assessed time-points in \code{y}, which will be inserted in the appropriate place).
}
  \item{alpha}{a multiplier controlling the size of the random steps between adjacent time points (e.g. year-to-year) in the constructed sequence. Precisely, the step size is \code{alpha * Q} where Q is the largest absolute change observed between adjacent time points in \code{y}, divided by the number of periods to be generated between those time points.
%%     ~~Describe \code{alpha} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
This function is called internally by \code{interpolate_betw_elic_times} if \code{method='rwalk.from.points'}.
}
\value{
a (1 x E) vector containing a (hopefully) reasonable random step size for each time period.
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
# compute random steps
compute_random_step(xe,nperiods=c(5,14,14),alpha=5)
# would normally be called as part of interpolate_betw_elic_times function
xt_rfp = interpolate_betw_elic_times(y=xe,nperiods=c(5,14,14),
method='rwalk.from.points',alpha=5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
