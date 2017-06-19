\name{projections-package}
\alias{projections-package}
\alias{projections}
\docType{package}
\title{Tools for generating and visualizing samples from elicited probability distributions
}
\description{
Accompanies the working papers "Using expert elicitation to characterize input uncertainty in E3 forecasting models" and "Probabilistic projections of baseline greenhouse gas emissions in South Africa to 2050". Contains tools for generating temporal projections from elicited (or otherwise constructed) trapezoidal probability distributions. Generally these distributions will be time-ordered, with some inter-temporal relationships specified as correlations or mean reversion. The package has two main goals: (1) generating samples at the E assessed time-points in such a way that inter-temporal information is included; (2) generating interpolated values for time-points occuring between the assessed time-points, with the idea that it will often not be possible to elicit information at all time points. Some functions for checking and visualization of the generated sequences are also provided.
}
\details{
\tabular{ll}{
Package: \tab projections\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-11-04\cr
License: \tab GPL-2\cr
}
The package would usually employ the functions below, in the order given. It is assumed that at each of E assessed time points, parameters of a trapezoidal probability distribution have been assessed in some way. \cr\cr
  \code{simulate_at_elic_times}: Generate samples from elicited distributions, including inter-temporal information. \cr\cr
  \code{check_keypts}: Run checks on samples to see they match desired distributional properties. \cr\cr
  \code{interpolate_betw_elic_times}: Generate trajectories from each sample, by interpolating between assessed time-points. \cr\cr
  \code{check_allpts}: Run checks on trajectories to see they match desired properties. \cr\cr
  \code{plot_trajectories}: Visualize simulated trajectories. \cr\cr
}
\author{
Ian Durbach

Maintainer: Ian Durbach <indurbach@gmail.com>
}
\references{
Durbach, I, Merven, B, McCall, B. (2015). Using expert elicitation to characterize input uncertainty in E3 forecasting models. MAPS working paper. \cr\cr
Merven, B, Durbach, I, McCall, B. (2015). Probabilistic projections of baseline greenhouse gas emissions in South Africa to 2050. MAPS working paper.
}
\keyword{ package }
\examples{
# create some data
x_min = c(2.4,2.5,1.5,1.5)
x_mode1 = c(2.6,3,3,3)
x_mode2 = x_mode1
x_max <- c(2.8,3.5,6,6)
correls = matrix(0.3,4,4)
diag(correls) = 1

# construct correlated samples at assessed time points
xe = simulate_at_elic_times(nsims=200,Sigma=correls,x_min,x_mode1,
x_mode2,x_max,mean_rev=c(0,0,1,0))

# check samples at assessed time points
check1 = check_keypts(xe)
# view output
check1$correlations
check1$histograms
check1$trajectories

# interpolate between assessed time points
xt = interpolate_betw_elic_times(y=xe,nperiods=c(5,14,14),method='linear.from.points')

# check generated sequences
check2 = check_allpts(xt,nperiods=c(7,15,15))
# view output
check2$histogram_points
check2$all_trajectories

# plot generated sequences
plot_trajectories(X=xt,xlabs=2014:2050)
}
