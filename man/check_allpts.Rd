\name{check_allpts}
\alias{check_allpts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Checks for generated sequences 
%%  ~~function to do ... ~~
}
\description{Takes the output of \code{interpolate_between_elic_times} and creates a 
number of summaries of the simulated samples: correlations, marginal distributions, and sample trajectories.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
check_allpts(y, nperiods)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{a (T x N) matrix or data frame, where T is the total number of time points used and N is the number of simulated samples. Each column should contain a single projection.
%%     ~~Describe \code{y} here~~
}
  \item{nperiods}{a vector indicating how the T time points should be aportioned into time `windows'. It indicates the number of time points in each of the observation periods for which summaries are to be produced.
%%     ~~Describe \code{nperiods} here~~
}
}
\details{Returns summary statistics and plots from simulated samples (all time points).
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{histograms_means }{Histogram of the mean values in each time window}
  \item{histograms_points }{Histogram of the values in each time window}
  \item{sample_trajectories }{A line plot of a small number of trajectories}
  \item{all_trajectories }{A line plot of all trajectories}
%% ...
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
# check generated sequences
check2 = check_allpts(xt,nperiods=c(7,15,15))
# view output
check2$histogram_points
check2$all_trajectories
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
