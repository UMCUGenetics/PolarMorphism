#############################################################
#   This is a modified version of the pvonmises function    #
#   originally written by:                                  #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Modified by:                                            #
#   Author: Joanna von Berg                                 #
#   Email: j.vonberg-24@umcutrecht.nl                       #
#############################################################


pvm.mu0.2 <- function(q, kappa, tol) {
  flag <- TRUE
  p <- 1
  sum <- 0
  while (flag) {
    term <- (besselI(x=kappa, nu=p, expon.scaled = FALSE) * sin(p*q))/p
    sum <- sum + term
    p <- p + 1
    if (all(abs(term) < tol))
      flag <- FALSE
  }
  return(q/(2*pi) + sum/(pi * besselI(x=kappa, nu=0, expon.scaled = FALSE)))
}

#' PvonmisesRad
#'
#' PvonmisesRad does things but then faster. (write later)
#' @export
PvonmisesRad.2 <- function(q, kappa, tol) {
  q <- q %% (2 * pi)
  n <- length(q)
  # mu <- mu %% (2 * pi)
  mu <- 0
  result <- rep(NA, n)
  if (mu == 0) {
    result <- pvm.mu0.2(q, kappa, tol)
  } else {
    # result <- rep(NA, n)
    # lt.mu.lo <- pvm.mu0((-mu)%%(2*pi), kappa, tol)
    # eq.mu.lo <- 1 - lt.mu.lo
    # gt.mu.lo <- 1 - lt.mu.lo
    # result[q < mu] <- pvm.mu0((q[q < mu] - mu) %% (2 * pi), kappa, tol) - lt.mu.lo
    # result[q == mu] <- eq.mu.lo
    # result[q > mu] <- pvm.mu0(q[q > mu] - mu, kappa, tol) + gt.mu.lo
  }
  return(result)
}



