\name{check_keypts}
\alias{check_keypts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Checks for generated samples at assessed time-points
%%  ~~function to do ... ~~
}
\description{Takes the output of \code{simulate_at_elic_times} and creates a 
number of summaries of the simulated samples: correlations, marginal distributions, and sample trajectories.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
check_keypts(y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{an (N x E) matrix or data frame where N is the number of simulations and E is the number of assessed time points. Each column should contain the samples from an assessed probability distribution.
%%     ~~Describe \code{y} here~~
}
}
\details{Returns summary statistics and plots from simulated samples (assessed time points).
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{correlations }{Matrix of correlations between time periods}
\item{histograms }{Histogram of values at each assessed time point}
\item{trajectories }{Line plot showing a few random trajectories}
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
# check samples at assessed time points
check1 = check_keypts(xe)
# view output
check1$correlations
check1$histograms
check1$trajectories

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
