\name{plot_trajectories}
\alias{plot_trajectories}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plots simulated trajectories
%%  ~~function to do ... ~~
}
\description{Plots all simulated trajectories with user-specified quantiles, using \code{ggplot2}.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot_trajectories(X, xlabs = NULL, myprobs = c(0.025, 0.1, 0.5, 0.9, 0.975), 
xtitle = "", ytitle = "", basesize = 12)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
a (T x N) matrix or data frame, where T is the total number of time points used and N is the number of simulated samples. Each column contains a single projection.
}
  \item{xlabs}{
%%     ~~Describe \code{xlabs} here~~
a (1 x T) vector of labels for tick marks on the x-axis, usually indicators of time (e.g. years). Must contain the same number of time-points as \code{X}. If not specified defaults to 1:T.
}
  \item{myprobs}{
%%     ~~Describe \code{myprobs} here~~
desired quantiles to plot. Defaults to \code{c(0.025, 0.1, 0.5, 0.9, 0.975)}. The two outer quantiles (2.5, 97.5) are shown in red; the two inner quantiles (10, 90) are shown in blue; the central quantile (50) is shown in green. Must be of length 5.
}
  \item{xtitle}{
%%     ~~Describe \code{xtitle} here~~
title for x-axis.
}
  \item{ytitle}{
%%     ~~Describe \code{ytitle} here~~
title for y-axis.
}
  \item{basesize}{
%%     ~~Describe \code{ytitle} here~~
base size for plot text.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
a \code{ggplot2} line plot
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
xt = interpolate_betw_elic_times(y=xe,nperiods=c(5,14,14),method='linear.from.points')
# plot generated sequences
plot_trajectories(X=xt,xlabs=2014:2050)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
