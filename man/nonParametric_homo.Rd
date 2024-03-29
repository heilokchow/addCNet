% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nonParamteric_homo.R
\name{nonParametric_homo}
\alias{nonParametric_homo}
\title{Non Parametric Estimate for homophily effects only}
\usage{
nonParametric_homo(data, zij, n, p, tz = 1, h1 = 0.05, test = 0, pvpval = 0)
}
\arguments{
\item{data}{The whole dataset containing pairwise communications, the first column
is the person with outgoing communication and the second column is the person with
incoming communication. The last column is the communication time point (should be
in numeric format).}

\item{zij}{The transfromed non-homogeneous effect (array with dimension (n,n,p,t)). The
transformation can be done using zProjection() function.}

\item{n}{number of individuals in the whole network.}

\item{p}{number of non-homogeneuous effect.}

\item{tz}{Time points for non-homofily effects (ranging from (0,1)), the length of this
vector should equal to the 4th dimension of zij which reflect the non-homofily effect
on or before those specific time points.}

\item{h1}{Bandwidth for calculate \eqn{\alpha(t)} and \eqn{\beta(t)}, default is 0.05.}
}
\value{
\item{nonhomo_coefficients}{List nonhomogeneous effect from scaled time range
[0,1]. The parameters are estimated at (0.01, 0.02, ..., 1)}
\item{nonhomo_sd}{List nonhomogeneous parameter's standard deviation for cumulative
effect from scaled time range [0,1]. The parameters are estimated at (0.01, 0.02, ..., 1)}
}
\description{
This function allows you to do non parametric estimation for additive network model's homophily effects (reduced model)
}
