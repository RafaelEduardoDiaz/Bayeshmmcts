% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BayesHMM_F.R
\name{pois.HMM.pseudo_residuals}
\alias{pois.HMM.pseudo_residuals}
\title{Pseudo-residuals}
\usage{
pois.HMM.pseudo_residuals(o, mod)
}
\arguments{
\item{o}{observed series.}

\item{mod}{a Poisson hidden Markov model.}
}
\value{
A matrix with lower, mid and upper normal pseudo-residuals.
}
\description{
Find ordinary normal pseudo-residuals, using Conditional probabilities
from \code{pois.HMM.conditional} function.
}
\details{
The form of the output is this. For each time \eqn{t} from 1 to \eqn{n}, the function provides the lower and upper normal pseudo-residuals, \eqn{z_t^{-}} and \eqn{z_t^{+}}, and the mid-pseudo-residual  \eqn{z_t^m = \Phi^{-1}((u_t^{-} + u_t^{+})/2)}.
}
\examples{

set.seed(2019)
o <- c(rpois(n = 100, lambda = 15),rpois(n = 50, lambda = 25))
mod_2states <- pois.HMM.mle(o = o, m = 2, lambda0 = c(20 ,25),
                            A0 = matrix(c(0.9, 0.1, 0.1, 0.9), 2, 2, byrow = TRUE),
                            pi0 = c(0.5,0.5), stationary = TRUE)
res <- pois.HMM.pseudo_residuals(o = o, mod = mod_2states)
print(res)
}
