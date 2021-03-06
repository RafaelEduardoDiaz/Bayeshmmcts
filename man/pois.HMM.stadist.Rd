% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BayesHMM_F.R
\name{pois.HMM.stadist}
\alias{pois.HMM.stadist}
\title{Stationary distribution}
\usage{
pois.HMM.stadist(mod)
}
\arguments{
\item{mod}{a Poisson hidden Markov model.}
}
\value{
A stationary distribution vector.
}
\description{
Computing stationary distribution using a transition probability matrix.
}
\examples{

set.seed(2019)
o <- c(rpois(n = 100, lambda = 15),rpois(n = 50, lambda = 25))
mod_2states <- pois.HMM.mle(o = o, m = 2, lambda0 = c(20 ,25),
                            A0 = matrix(c(0.9, 0.1, 0.1, 0.9), 2, 2, byrow = TRUE),
                            pi0 = c(0.5,0.5), stationary = TRUE)
pois.HMM.stadist(mod = mod_2states)

}
