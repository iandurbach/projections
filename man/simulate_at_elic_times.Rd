\name{simulate_at_elic_times}
\alias{simulate_at_elic_times}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generate samples at assessed time-points
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Generates possibly correlated (or otherwise dependent) sequences at elicited time-points. Assumes that, at each of E assessed time points, parameters of a trapezoidal probability distribution have been assessed. Inter-temporal relationships between assessed time-points are represented by correlations or as a statement of "mean reversion", where this means that extreme values in one time period must be followed by less extreme values in the next.
}
\usage{
simulate_at_elic_times(nsims = 1000, Sigma, x_min, x_mode1, x_mode2, 
x_max, mean_rev = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nsims}{
%%     ~~Describe \code{nsims} here~~
number of samples to generate
}
  \item{Sigma}{
%%     ~~Describe \code{Sigma} here~~
a (E x E) matrix of correlations, containing any correlations between assessed time-points
}
  \item{x_min}{
%%     ~~Describe \code{x_min} here~~
a (1 x E) vector containing the "minimum" parameter of each assessed trapezoidal probability distribution.
}
  \item{x_mode1}{
%%     ~~Describe \code{x_mode1} here~~
a (1 x E) vector containing the "lower mode" parameter of each assessed trapezoidal probability distribution. This is the smallest value at which the trapezoidal distribution achieves its maximum probability.
}
  \item{x_mode2}{
%%     ~~Describe \code{x_mode2} here~~
a (1 x E) vector containing the "upper mode" parameter of each assessed trapezoidal probability distribution. This is the largest value at which the trapezoidal distribution achieves its maximum probability.
}
  \item{x_max}{
%%     ~~Describe \code{x_max} here~~
a (1 x E) vector containing the "maximum" parameter of each assessed trapezoidal probability distribution.
}
  \item{mean_rev}{
%%     ~~Describe \code{mean_rev} here~~
a (1 x E) vector containing a binary indicator of whether or not mean reversion applies in the following time period. If a 1 appears at assessed time point \code{i}, then extreme values in time period \code{i} will be followed by less extreme values in time period \code{i+1}. See details.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Mean reversion is implemented by matching extreme ranks in period \code{i} with intermediate ranks in period \code{i+1}. This means that the correlations specified between the two relevant time periods in \code{Sigma} will no longer apply -- the mean reversion supercedes the correlation. The rank reversal implementation of mean reversio still leaves some role for correlation to play though. For example, we can say that both extremely high and low values at \code{i} should be followed by intermediate values at \code{i+1}, but that the high values should appear above the low values (positive correlation) or below them (negative correlation). The same applies for intermediate values: values that were above-average at \code{i} will have a larger chance of being extremely high at \code{i+1} than values that were below-average. Thus the user can specify desired "mean reversion correlations" in the appropriate place in the \code{Sigma} i.e. in elements \code{(i,i+1)} and \code{(i+1,i)}.
}
\value{
%%  ~Describe the value returned
a (N x E) matrix where N is the number of simulations and E is the number of assessed time points. Each column contains the samples from an assessed probability distribution, after accounting for inter-temporal relationships between the assessed time points.
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
xe = simulate_at_elic_times(nsims=200,Sigma=correls,x_min,x_mode1,x_mode2,x_max)
head(xe)
# with mean reversion between t3 and t4
xe = simulate_at_elic_times(nsims=200,Sigma=correls,x_min,x_mode1,x_mode2,
x_max,mean_rev=c(0,0,1,0))
head(xe)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
