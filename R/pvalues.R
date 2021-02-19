#' PvalueForR
#'
#' PvalueForR performs p-value estimation using the chi distribution
#' @param r a vector of distances
#' @param p degrees of freedom (the number of original GWAS)
#' @param exp.r a vector of expected distances (for each SNP, the maximum of the absolute z-scores of the original GWASs)
#' @export
PvalueForR <- function(r, p, exp.r){
  chi.pval <- pchisq(q = r^2, df = p, ncp = exp.r^2, lower.tail = F)
  return(chi.pval)
}
#'
#'
#'
# Vpvonmises <- Vectorize(circular::pvonmises, vectorize.args = c("q", "kappa"))

#' PvalueForTheta
#'
#' PvalueForTheta performs p-value estimation using the von Mises distribution
#' @param theta a vector of angles
#' @param kappa a vector of kappa (either estimated yourself or provided by EstimateKappa)
#' @export
PvalueForTheta <- function(relative.theta, kappa, tol){
  vm.pval <- circular::pvonmises(q = relative.theta, mu = 0, kappa = kappa, tol = tol)/4
  return(vm.pval)
}
