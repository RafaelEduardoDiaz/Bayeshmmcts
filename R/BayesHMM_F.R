#' Transforming natural parameters to working
#'
#' The purpose of the transformation of the natural parameters (Poisson means,
#' transition probabilities and, if appropriate, the initial distribution of the
#' Markov chain) to working parameters is simply to convert
#' a constrained optimization problem to an unconstrained one.
#'
#' @param m the number of states.
#' @param lambda the estimates of \eqn{\lambda}.
#' @param A the estimates of \eqn{A}.
#' @param pi the estimates of \eqn{\pi}.
#' @param stationary indicate if the Markov chain is stationary, by default is \code{TRUE}.
#' @return parvect parameters vector of working.
#' @examples
#'
#' m0 = 2
#' lambda0 = c(15, 25)
#' A0 = matrix(c(0.9, 0.1, 0.1, 0.9), 2, 2, byrow = TRUE)
#' pi0 = c(0.5,0.5)
#' parvect_work <- pois.HMM.pn2pw(m = m0, lambda = lambda0, A = A0, pi = pi0, stationary = TRUE)
#' print(parvect_work)
#'
#' @export
pois.HMM.pn2pw <- function(m, lambda, A, pi = NULL, stationary = TRUE){
  tlambda <- log(lambda)
  if(m == 1) return(tlambda)
  foo <- log(A / diag(A))
  tA <- as.vector(foo[!diag(m)])
  if(stationary) {tpi <- NULL}
  else{tpi <- log(pi[-1] / pi[1])}
  parvect <- c(tlambda, tA, tpi)
  return(parvect)
}

#' Transforming working parameters to natural
#'
#' Transform the estimates of the working parameters to estimates of
#' the natural parameters.
#'
#' @param m the number of states.
#' @param parvect parameters vector of working from \code{pois.HMM.pn2pw} function.
#' @param stationary indicate if the Markov chain is stationary, by default is \code{TRUE}.
#' @return A list of parameters
#' \itemize{
#'   \item lambda - estimates values of \eqn{\lambda}.
#'   \item A - estimates values of \eqn{A}.
#'   \item pi - estimates values of \eqn{\pi}.
#' }
#' @examples
#'
#' m0 = 2
#' lambda0 = c(15, 25)
#' A0 = matrix(c(0.9, 0.1, 0.1, 0.9), 2, 2, byrow = TRUE)
#' pi0 = c(0.5,0.5)
#' parvect_work <- pois.HMM.pn2pw(m = m0, lambda = lambda0, A = A0, pi = pi0, stationary = TRUE)
#' pois.HMM.pw2pn(m = 2, parvect = parvect_work)
#'
#' @export
pois.HMM.pw2pn <- function(m, parvect, stationary = TRUE ){
  lambda <- exp(parvect[1:m])
  A <- diag(m)
  if(m == 1) return(list(lambda = lambda , A = A , pi = 1))
  A[!A] <- exp(parvect[(m + 1):(m * m)])
  A <- A / apply(A ,1, sum)
  if(stationary){pi <- solve(t(diag(m) - A +1), rep(1,m))}
  else{foo <-c(1, exp(parvect[(m*m +1) :(m * m + m - 1)]))
  pi <- foo / sum(foo)}
  return(list(lambda = lambda, A = A, pi = pi))
}

#' Computing minus the log-likelihood from the working parameters
#'
#' Compute minus the log-likelihood of a Poisson HMM for given values of
#' the working parameters.
#'
#' @param parvect parameters vector of working from \code{pois.HMM.pn2pw} function.
#' @param o a non-negative vector with of count data time serie.
#' @param m the number of states.
#' @param stationary indicate if the Markov chain is stationary, by default is \code{TRUE}.
#' @return mllk Minus the log-likelihood.
#' @importFrom stats dpois
#' @examples
#'
#' m0 = 2
#' lambda0 = c(15, 25)
#' A0 = matrix(c(0.9, 0.1, 0.1, 0.9), 2, 2, byrow = TRUE)
#' pi0 = c(0.5,0.5)
#' parvect_work <- pois.HMM.pn2pw(m = m0, lambda = lambda0, A = A0, pi = pi0, stationary = TRUE)
#' o <- c(rpois(n = 100, lambda = 15),rpois(n = 50, lambda = 25))
#' pois.HMM.mllk(parvect = parvect_work, o = o, m = 2, stationary = TRUE)
#'
#' @export
pois.HMM.mllk <- function(parvect, o, m, stationary = TRUE){
  if(m == 1) return(-sum(dpois(o, exp(parvect), log = TRUE)))
  n <- length(o)
  pn <- pois.HMM.pw2pn(m, parvect, stationary = stationary )
  foo <- pn$pi * dpois(o[1] , pn$lambda)
  sumfoo <- sum(foo)
  lscale <- log(sumfoo)
  foo <- foo / sumfoo
  for(i in 2:n){
    if(!is.na(o[i])){P <- dpois(o[i], pn$lambda)}
    else{P <- rep(1,m)}
    foo <- foo %*% pn$A * P
    sumfoo <- sum(foo)
    lscale <- lscale + log(sumfoo)
    foo <- foo / sumfoo
  }
  mllk <- -lscale
  return(mllk)
}

#' Computing the MLEs, given starting values for the natural parameters
#'
#' Estimate the (natural) parameters of the model by using numerical minimization
#' of minus the log-likelihood. For more details see the function `nlm{stats}`.
#'
#' @param o a non-negative vector with of count data time serie.
#' @param m the number of states.
#' @param lambda0 starting value of \eqn{\lambda}, the vector of state-dependent means in a Poisson HMM.
#' @param A0 starting value of \eqn{A}, the transition probability matrix.
#' @param pi0 starting value of \eqn{\pi}, the stationary distribution of the Markov chain.
#' @param stationary indicate if the Markov chain is stationary, by default is \code{TRUE}.
#' @param ... Arguments passed to `stats::nlm` (e.g. gradtol, steptol, iterlim).
#' @return A list of parameters
#' \itemize{
#'   \item m - the number of states.
#'   \item lambda - estimates values of \eqn{\lambda}.
#'   \item A - estimates values of \eqn{A}.
#'   \item pi - estimates values of \eqn{\pi}.
#'   \item code - an integer indicating why the optimization process terminated. For more details see \code{nlm{stats}}.
#'   \item mllk - Minus the log-likelihood.
#'   \item AIC - The Akaike information criterion (AIC).
#'   \item BIC - the Bayesian information criterion (BIC)
#' }
#' @importFrom stats nlm
#' @examples
#'
#' set.seed(2019)
#' o <- c(rpois(n = 100, lambda = 15),rpois(n = 50, lambda = 25))
#' acf(o)
#' mod_2states <- pois.HMM.mle(o = o, m = 2, lambda0 = c(20 ,25),
#'                             A0 = matrix(c(0.9, 0.1, 0.1, 0.9), 2, 2, byrow = TRUE),
#'                             pi0 = c(0.5,0.5), stationary = TRUE)
#' print(mod_2states)
#'
#' @export
pois.HMM.mle <- function(o, m, lambda0, A0, pi0 = NULL, stationary = TRUE,...){
  parvect0 <- pois.HMM.pn2pw(m, lambda0, A0, pi0, stationary = stationary)
  mod <- nlm(pois.HMM.mllk, parvect0, o = o, m = m, stationary = stationary)
  pn <- pois.HMM.pw2pn(m = m, mod$estimate, stationary = stationary)
  mllk <- -mod$minimum
  p <- length(parvect0)
  AIC <- -2*mllk + 2*p
  n <- sum(!is.na(o))
  BIC <- -2*mllk + p*log(n)
  list(m = m, lambda = pn$lambda, A = pn$A, pi = pn$pi, code = mod$code,
       mllk = mllk, AIC = AIC, BIC = BIC)
}

#' Generating a sample
#'
#' This function generates a realization, of length ns, of the PHMM mod.
#'
#' @param ns size of sample.
#' @param mod a Poisson hidden Markov model, from \code{pois.HMM.mle} function.
#' @return x a random values from Poisson hidden Markov model.
#' @importFrom stats rpois
#' @examples
#'
#' set.seed(2019)
#' o <- c(rpois(n = 100, lambda = 15),rpois(n = 50, lambda = 25))
#' mod_2states <- pois.HMM.mle(o = o, m = 2, lambda0 = c(20 ,25),
#'                             A0 = matrix(c(0.9, 0.1, 0.1, 0.9), 2, 2, byrow = TRUE),
#'                             pi0 = c(0.5,0.5), stationary = TRUE)
#' pois.HMM.generate_sample(ns = 10, mod = mod_2states)
#'
#' @export
pois.HMM.generate_sample <- function(ns , mod){
  mvect <- 1:mod$m
  state <- numeric(ns)
  state [1] <- sample(mvect, 1, prob = mod$pi)
  for(i in 2: ns) state[i] <- sample(mvect, 1, prob = mod$A[state[i - 1], ])
  o <- rpois(ns, lambda = mod$lambda[state])
  return(o)
}

#' Global decoding by the Viterbi algorithm
#'
#' Given the model mod and observed series x, this function performs global decoding.
#'
#' @param o observed series.
#' @param mod a Poisson hidden Markov model.
#' @return iv result of Viterbi algorithm.
#' @examples
#'
#' set.seed(2019)
#' o <- c(rpois(n = 100, lambda = 15),rpois(n = 50, lambda = 25))
#' mod_2states <- pois.HMM.mle(o = o, m = 2, lambda0 = c(20 ,25),
#'                             A0 = matrix(c(0.9, 0.1, 0.1, 0.9), 2, 2, byrow = TRUE),
#'                             pi0 = c(0.5,0.5), stationary = TRUE)
#' pois.HMM.viterbi(o = o, mod = mod_2states)
#'
#' @export
pois.HMM.viterbi <- function(o, mod){
  n <- length(o)
  oi <- matrix(0, n, mod$m)
  foo <- mod$pi * dpois(o[1], mod$lambda)
  oi[1 ,] <- foo / sum(foo)
  for(i in 2:n){
    foo <- apply(oi[i - 1, ] * mod$A, 2, max) * dpois(o[i], mod$lambda)
    oi[i, ] <- foo / sum(foo)
  }
  iv <- numeric(n)
  iv[n] <- which.max(oi[n, ])
  for(i in (n - 1):1) iv[i] <- which.max(mod$A[ ,iv[i + 1]] * oi[i, ])
  return(iv)
}

#' Computing log(forward probabilities)
#'
#' Given data o and model mod, this function uses the recursion \eqn{\alpha_{t+1} = \alpha_t A P(o_{t+1})}
#' to find all the vectors of forward probabilities, in logarithmic form.
#'
#' @param o observed series.
#' @param mod a Poisson hidden Markov model.
#' @return A matrix (lalpha) is returned.
#' @importFrom stats dpois
#' @examples
#'
#' set.seed(2019)
#' o <- c(rpois(n = 100, lambda = 15),rpois(n = 50, lambda = 25))
#' mod_2states <- pois.HMM.mle(o = o, m = 2, lambda0 = c(20 ,25),
#'                             A0 = matrix(c(0.9, 0.1, 0.1, 0.9), 2, 2, byrow = TRUE),
#'                             pi0 = c(0.5,0.5), stationary = TRUE)
#' pois.HMM.lforward(o = o, mod = mod_2states)
#'
#' @export
pois.HMM.lforward <- function(o, mod){
  n <- length(o)
  lalpha <- matrix(NA, n, mod$m)
  foo <- mod$pi * dpois(o[1], mod$lambda)
  sumfoo <- sum(foo)
  lscale <- log(sumfoo)
  foo <- foo / sumfoo
  lalpha[1, ] <- lscale + log(foo)
  for(i in 2:n){
    foo <- foo %*% mod$A * dpois(o[i], mod$lambda)
    sumfoo <- sum(foo)
    lscale <- lscale + log(sumfoo)
    foo <- foo / sumfoo
    lalpha[i, ] <- log(foo) + lscale
  }
  return(lalpha)
}

#' Computing log(backward probabilities)
#'
#' Similarly a log(forward probabilities), this function uses the recursion
#' \eqn{\beta'_t = A P(o_{t+1} \beta'_{t+1})} to find all the vectors of
#' backward probabilities, in logarithmic form.
#'
#' @param o observed series.
#' @param mod a Poisson hidden Markov model.
#' @return A matrix (lalpha) is returned.
#' @importFrom stats dpois
#' @examples
#'
#' set.seed(2019)
#' o <- c(rpois(n = 100, lambda = 15),rpois(n = 50, lambda = 25))
#' mod_2states <- pois.HMM.mle(o = o, m = 2, lambda0 = c(20 ,25),
#'                             A0 = matrix(c(0.9, 0.1, 0.1, 0.9), 2, 2, byrow = TRUE),
#'                             pi0 = c(0.5,0.5), stationary = TRUE)
#' pois.HMM.lbackward(o = o, mod = mod_2states)
#'
#' @export
pois.HMM.lbackward <- function(o, mod){
  n <- length(o)
  m <- mod$m
  lbeta <- matrix(NA, n, m)
  lbeta[n, ] <- rep(0, m)
  foo <- rep(1/m, m)
  lscale <- log(m)
  for(i in (n - 1):1){
    foo <- mod$A %*% (dpois(o[i + 1] , mod$lambda) * foo)
    lbeta[i, ] <- log(foo) + lscale
    sumfoo <- sum(foo)
    foo <- foo / sumfoo
    lscale <- lscale + log(sumfoo)
  }
  return(lbeta)
}

#' Conditional probabilities
#'
#' Conditional probability that observation at time \eqn{t} equals \code{oc},
#' given all observations other than that at time \eqn{t}.
#'
#' @param oc A input specifies the range of o-values for which these probabilities are required.
#' @param o observed series.
#' @param mod a Poisson hidden Markov model.
#' @return the result is a matrix with conditional probabilities.
#' @importFrom stats dpois
#' @examples
#'
#' set.seed(2019)
#' o <- c(rpois(n = 100, lambda = 15),rpois(n = 50, lambda = 25))
#' mod_2states <- pois.HMM.mle(o = o, m = 2, lambda0 = c(20 ,25),
#'                             A0 = matrix(c(0.9, 0.1, 0.1, 0.9), 2, 2, byrow = TRUE),
#'                             pi0 = c(0.5,0.5), stationary = TRUE)
#' pois.HMM.conditional(oc = 10:15, o = o, mod = mod_2states)
#'
#' @export
pois.HMM.conditional <- function(oc, o, mod){
  n <- length(o)
  m <- mod$m
  noc <- length(oc)
  doc <- matrix(NA, nrow = noc, ncol = n)
  Po <- matrix (NA, nrow = m, ncol = noc)
  for(j in 1:noc) Po[, j] <- dpois(oc[j], mod$lambda)
  la <- t(pois.HMM.lforward(o, mod))
  lb <- t(pois.HMM.lbackward(o, mod))
  la <- cbind(log(mod$pi), la)
  lafact <- apply(la, 2, max)
  lbfact <- apply(lb, 2, max)
  for(i in 1:n){
    foo <- (exp(la[, i] - lafact[i]) %*% mod$A) * exp(lb[, i] - lbfact[i])
    foo <- foo / sum(foo)
    doc[, i] <- foo %*% Po
  }
  return(doc)
}

#' Pseudo-residuals
#'
#' Find ordinary normal pseudo-residuals, using Conditional probabilities
#' from \code{pois.HMM.conditional} function.
#'
#' @param o observed series.
#' @param mod a Poisson hidden Markov model.
#' @return A matrix with lower, mid and upper normal pseudo-residuals.
#' @importFrom stats qnorm
#' @details The form of the output is this. For each time \eqn{t} from 1 to \eqn{n}, the function provides the lower and upper normal pseudo-residuals, \eqn{z_t^{-}} and \eqn{z_t^{+}}, and the mid-pseudo-residual  \eqn{z_t^m = \Phi^{-1}((u_t^{-} + u_t^{+})/2)}.
#' @examples
#'
#' set.seed(2019)
#' o <- c(rpois(n = 100, lambda = 15),rpois(n = 50, lambda = 25))
#' mod_2states <- pois.HMM.mle(o = o, m = 2, lambda0 = c(20 ,25),
#'                             A0 = matrix(c(0.9, 0.1, 0.1, 0.9), 2, 2, byrow = TRUE),
#'                             pi0 = c(0.5,0.5), stationary = TRUE)
#' res <- pois.HMM.pseudo_residuals(o = o, mod = mod_2states)
#' print(res)
#' @export
pois.HMM.pseudo_residuals <- function(o, mod){
  n <- length(o)
  cdists <- pois.HMM.conditional(oc = 0:max(o), o, mod)
  cumdists <- rbind(rep(0, n), apply(cdists, 2,cumsum))
  ulo <- uhi <- rep(NA ,n)
  for(i in 1:n){
    ulo[i] <- cumdists[o[i] + 1 ,i]
    uhi[i] <- cumdists[o[i] + 2 ,i]
  }
  umi <- 0.5 *(ulo + uhi)
  npsr <- qnorm(rbind(ulo, umi, uhi))
  return(data.frame(t(npsr)))
}

#' Plotting pseudo-residuals
#'
#' Plot of the ordinary normal mid-pseudo-residuals
#'
#' @param residual is a ordinary normal pseudo-residuals.
#' @param ... further arguments and graphical parameters passed to plot.histogram and thence to title and axis (if plot = TRUE).
#' @return four plots (index, qq plot, histogram and ACF).
#' @importFrom graphics par
#' @importFrom graphics abline
#' @importFrom stats qqnorm
#' @importFrom stats qqline
#' @importFrom graphics hist
#' @importFrom graphics curve
#' @importFrom stats dnorm
#' @importFrom stats sd
#' @importFrom stats acf
#' @details The top left shows index plot of the normal pseudo-residuals, with horizontal lines at 0, \eqn{\pm 1.96}, \eqn{\pm 2.58}. The top right shows quantile-quantile plots of the normal pseudo-residuals, with the theoretical quantiles on the horizontal axis. The bottom left show histograms of the normal pseudo-residuals. The last plot shows the autocorrelation functions of the normal pseudo-residuals.
#' @examples
#'
#' set.seed(2019)
#' o <- c(rpois(n = 100, lambda = 15),rpois(n = 50, lambda = 25))
#' mod_2states <- pois.HMM.mle(o = o, m = 2, lambda0 = c(20 ,25),
#'                             A0 = matrix(c(0.9, 0.1, 0.1, 0.9), 2, 2, byrow = TRUE),
#'                             pi0 = c(0.5,0.5), stationary = TRUE)
#' res <- pois.HMM.pseudo_residuals(o = o, mod = mod_2states)
#' pois.HMM.plot.residuals(res)
#'
#' @export
pois.HMM.plot.residuals <- function(residual,...){
  par(mfrow = c(2,2), mar = c(3,3,2,2), mgp = c(2,1,0))
  # plot 1
  plot(x = seq_along(residual$umi), y = residual$umi, pch = 22, col = "black", bg = "grey", ylim = c(-4,4), ylab = "", xlab = "")
  abline(h = 0, col = "grey", lty = 2)
  abline(h = c(-1.96,1.96), col = "grey", lty = 2)
  abline(h = c(-2.58,2.58), col = "grey", lty = 2)
  # plot 2
  qqnorm(residual$umi, main = "")
  qqline(residual$umi, lwd = 2, col = "grey")
  # plot 3
  x <- NULL
  hist(residual$umi, prob = TRUE,col = "grey", border = "white", main = "", xlab = "",...)
  curve(dnorm(x, mean=mean(residual$umi), sd=sd(residual$umi)), add=TRUE, lwd = 2)
  # plot 4
  acf(residual$umi, main = "")
}

#' State probabilities
#'
#' Computing state probabilities using log(forward) and log(backward) probabilities.
#'
#' @param o observed series.
#' @param mod a Poisson hidden Markov model.
#' @param digits integer indicating the number of decimal places (round).
#' @return A matrix with state probabilities.
#' @details Here we compute probabilities \eqn{Pr(C_t = i | O^{(T)} = o^{(T)})}, for \eqn{t \in \lbrace 1,2,...,T \rbrace}.
#' @examples
#'
#' set.seed(2019)
#' o <- c(rpois(n = 100, lambda = 15),rpois(n = 50, lambda = 25))
#' mod_2states <- pois.HMM.mle(o = o, m = 2, lambda0 = c(20 ,25),
#'                             A0 = matrix(c(0.9, 0.1, 0.1, 0.9), 2, 2, byrow = TRUE),
#'                             pi0 = c(0.5,0.5), stationary = TRUE)
#' pois.HMM.state_probs(o = o, mod = mod_2states)
#'
#' @export
pois.HMM.state_probs <- function(o, mod, digits = 7){
  n <- length(o)
  la <- pois.HMM.lforward(o, mod)
  lb <- pois.HMM.lbackward(o, mod)
  c <- max(la[n, ])
  llk <- c + log(sum(exp(la[n, ] - c)))
  stateprobs <- matrix(NA, nrow = n, ncol = mod$m)
  for(i in 1:n) stateprobs[i, ] <- exp(la[i, ] + lb[i, ] - llk)
  stateprobs <- as.data.frame(stateprobs)
  colnames(stateprobs) <- paste0("State",seq_len(mod$m))
  return(round(stateprobs, digits = digits))
}

#' State prediction
#'
#' State prediction using a PHMM.
#'
#' @param h Number of periods for forecasting.
#' @param o observed series.
#' @param mod a Poisson hidden Markov model.
#' @return The state output is a matrix even if h=1.
#' @details Compute probabilities using equation \eqn{Pr(C_{T+h} = i | O^{(T)} = o^{(T)})}, for a range of values \eqn{h \in \lbrace 1,2,... \rbrace}.
#' @examples
#'
#' set.seed(2019)
#' o <- c(rpois(n = 100, lambda = 15),rpois(n = 50, lambda = 25))
#' mod_2states <- pois.HMM.mle(o = o, m = 2, lambda0 = c(20 ,25),
#'                             A0 = matrix(c(0.9, 0.1, 0.1, 0.9), 2, 2, byrow = TRUE),
#'                             pi0 = c(0.5,0.5), stationary = TRUE)
#' pois.HMM.state_prediction(h = 5, o = o, mod = mod_2states)
#'
#' @export
pois.HMM.state_prediction <- function(h = 1, o, mod){
  n <- length(o)
  la <- t(pois.HMM.lforward(o, mod))
  c <- max(la[, n])
  llk <- c + log(sum(exp(la[, n] - c)))
  statepreds <- matrix(NA, ncol = h, nrow = mod$m)
  foo <- exp(la[, n] - llk)
  for(i in 1:h){
    foo <- foo %*% mod$A
    statepreds[, i] <- foo
  }
  return(t(statepreds))
}

#' Local decoding
#'
#' This approach determines the most likely sequence of (hidden) state separately for each \eqn{t}.
#'
#' @param o observed series.
#' @param mod a Poisson hidden Markov model.
#' @return A vector with the most likely states.
#' @details The local decoding in this approach consists in determines the most likely state separately for each \eqn{t} by maximizing the conditional probability.
#' @examples
#'
#' set.seed(2019)
#' o <- c(rpois(n = 100, lambda = 15),rpois(n = 50, lambda = 25))
#' mod_2states <- pois.HMM.mle(o = o, m = 2, lambda0 = c(20 ,25),
#'                             A0 = matrix(c(0.9, 0.1, 0.1, 0.9), 2, 2, byrow = TRUE),
#'                             pi0 = c(0.5,0.5), stationary = TRUE)
#' pois.HMM.local_decoding(o = o, mod = mod_2states)
#'
#' @export
pois.HMM.local_decoding <- function(o, mod){
  n <- length(o)
  stateprobs <- pois.HMM.state_probs(o, mod)
  ild <- rep(NA, n)
  for(i in 1:n) ild[i] <- which.max(stateprobs[i, ])
  ild
}

#' Forecast probabilities
#'
#' Computing forecast probabilities using a PHMM.
#'
#' @param of The range of \code{o}-values for which these probabilities are required is specified.
#' @param h The range of times for which they are required.
#' @param o observed series.
#' @param mod a Poisson hidden Markov model.
#' @return output is a matrix with the forecast probabilities.
#' @importFrom stats dpois
#' @examples
#'
#' set.seed(2019)
#' o <- c(rpois(n = 100, lambda = 15),rpois(n = 50, lambda = 25))
#' mod_2states <- pois.HMM.mle(o = o, m = 2, lambda0 = c(20 ,25),
#'                             A0 = matrix(c(0.9, 0.1, 0.1, 0.9), 2, 2, byrow = TRUE),
#'                             pi0 = c(0.5,0.5), stationary = TRUE)
#' forecasts <- pois.HMM.forecast(of = 0:50, h = 1, o, mod_2states)
#' par(mfrow = c(1, 1), las = 1)
#' plot(0:50, forecasts, type = "h", main = paste("forecast distribution for", length(o) + 1),
#' xlim = c(0, max(0:50)), ylim = c(0, 0.12), xlab = "count", ylab = "probability", lwd = 3)
#'
#' @export
pois.HMM.forecast <- function(of, h=1, o, mod){
  n <- length(o)
  nof <- length(of)
  dof <- matrix(0, nrow = h, ncol = nof)
  foo <- mod$pi * dpois(o[1], mod$lambda)
  sumfoo <- sum(foo)
  lscale <- log(sumfoo)
  foo <- foo / sumfoo
  for(i in 2:n){
    foo <- foo %*% mod$A * dpois(o[i], mod$lambda)
    sumfoo <- sum(foo)
    lscale <- lscale + log(sumfoo)
    foo <- foo / sumfoo
  }
  for(i in 1:h){
    foo <- foo %*% mod$A
    for(j in 1:mod$m) dof[i, ] <- dof[i, ] + foo[j] * dpois(of, mod$lambda[j])
  }
  return(t(dof))
}

#' Stationary distribution
#'
#' Computing stationary distribution using a transition probability matrix.
#'
#' @param mod a Poisson hidden Markov model.
#' @return A stationary distribution vector.
#' @examples
#'
#' set.seed(2019)
#' o <- c(rpois(n = 100, lambda = 15),rpois(n = 50, lambda = 25))
#' mod_2states <- pois.HMM.mle(o = o, m = 2, lambda0 = c(20 ,25),
#'                             A0 = matrix(c(0.9, 0.1, 0.1, 0.9), 2, 2, byrow = TRUE),
#'                             pi0 = c(0.5,0.5), stationary = TRUE)
#' pois.HMM.stadist(mod = mod_2states)
#'
#' @export
pois.HMM.stadist <- function(mod){
  A <- mod$A
  m <- ncol(A)
  pi <- rep(1, m) %*% solve(diag(m) - A + matrix(1, m, m))
  return(c(pi))
}

#' Poisson HMM moments
#'
#' Computing moment from a Poisson hidden Markov model, with stationary Markov chain.
#'
#' @param mod a Poisson hidden Markov model.
#' @param lag.max Maximum lag for correlation function \eqn{\rho}.
#' @return A list of moments and rho
#' \itemize{
#'   \item mu - estimate value of first moment \eqn{\mu}.
#'   \item sigma - estimate value of variance \eqn{\sigma}.
#'   \item rho - estimates values of \eqn{\rho}, for a lag \eqn{h}.
#' }
#' @examples
#'
#' set.seed(2019)
#' o <- c(rpois(n = 100, lambda = 15),rpois(n = 50, lambda = 25))
#' mod_2states <- pois.HMM.mle(o = o, m = 2, lambda0 = c(20 ,25),
#'                             A0 = matrix(c(0.9, 0.1, 0.1, 0.9), 2, 2, byrow = TRUE),
#'                             pi0 = c(0.5,0.5), stationary = TRUE)
#' pois.HMM.moments(mod = mod_2states)
#'
#' @export
pois.HMM.moments <- function(mod, lag.max = 10) {
  requireNamespace("expm", quietly = TRUE)
  lambda <- mod$lambda
  A <- mod$A
  pi <- mod$pi
  EO <- pi %*% lambda
  VARO <- EO + pi %*% diag(lambda) %*% lambda - EO^2
  COR <- rep(NA, lag.max)
  for(i in 1:lag.max){
    nume <- (pi %*% diag(lambda) %*% expm::`%^%`(x = A, k = i) %*% lambda - EO^2)
    deno <- (pi %*% diag(lambda) %*% lambda + EO - EO^2)
    COR[i] <- nume / deno
  }
  return(list(mu = EO, sigma = VARO, rho = COR))
}

#' Poisson Hidden Markov Model
#'
#' Fit a Bayesian Poisson hidden Markov model with stan.
#'
#' @param y observed series.
#' @param m the number of states.
#' @param ... Arguments passed to \code{rstan::sampling} (e.g. iter, chains).
#' @return An object of class \code{stanfit} returned by \code{rstan::sampling}.
#' @examples
#'
#' set.seed(2019)
#' y <- c(rpois(n = 100, lambda = 15),rpois(n = 50, lambda = 25),rpois(n = 75, lambda = 32))
#' PHHMM_3states <- bayes.PHMM(y = y, m = 3, chains = 2, iter = 1000)
#' print(PHHMM_3states, digits_summary = 3)
#'
#' @export
bayes.PHMM <- function(y, m = 2, ...){
  standata <- list(N = length(y), m = m, y = y)
  out <- rstan::sampling(stanmodels$PHMM, data = standata, ...)
  return(out)
}

#' ZIP-HMM (Zero inflated State == 1)
#'
#' Fit a Bayesian Zero inflated Poisson hidden Markov model with stan,
#' where zero-inflation only happens in State 1.
#'
#' @param y observed series.
#' @param m the number of states.
#' @param ... Arguments passed to \code{rstan::sampling} (e.g. iter, chains).
#' @return An object of class \code{stanfit} returned by \code{rstan::sampling}.
#' @export
bayes.ZIPHMM1 <- function(y, m = 2, ...){
  standata <- list(N = length(y), m = m, y = y)
  out <- rstan::sampling(stanmodels$ZIPHMM_INF_1, data = standata,...)
  return(out)
}

#' ZIP-HMM (Zero inflated all States)
#'
#' Fit a Bayesian Zero inflated Poisson hidden Markov model with stan.
#'
#' @param y observed series.
#' @param m the number of states.
#' @param ... Arguments passed to \code{rstan::sampling} (e.g. iter, chains).
#' @return An object of class \code{stanfit} returned by \code{rstan::sampling}.
#' @export
bayes.ZIPHMM <- function(y, m = 2, ...){
  standata <- list(N = length(y), m = m, y = y)
  out <- rstan::sampling(stanmodels$ZIPHMM, data = standata,...)
  return(out)
}
