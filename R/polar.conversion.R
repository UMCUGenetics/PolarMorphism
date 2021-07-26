#' CalcLambda
#'
#' Calculates the great circle distance between a SNP's position and the expected position under the null hypothesis of trait-specific effect.
#' The great circle distance is the angle between the vector from the origin of the p-dimensional sphere the the SNP, and the vector from the origin to the expected position.
#' Works per SNP (row-wise).
#' @param inputvec vector of effect sizes (z-scores) of 1 SNP on p traits
#' @returns The angle in radians, not normalized to describe a full circle
CalcLambda <- function(inputvec){
  inputvec <- abs(inputvec)
  mindim <- which.max(inputvec)
  rsquared <- sum(inputvec^2)
  r <- sqrt(rsquared)
  zerovec <- rep(0, length(inputvec))
  zerovec[mindim] <- r
  distsquared <- sum((inputvec - zerovec)^2)
  res <- acos((1 - (distsquared/(2*rsquared))))
  return(res)
}

#' PolarCoords
#'
#' Given two vectors for x- and y-coordinates, performs polar transformation. Returns the angle theta as a number between -pi and +pi.
#' @param x,y vector of coordinates
#' @export
PolarCoords <- function(x, y, z = NULL){
  if(is.null(z)){
    theta <- atan2(y = y, x = x)%%(2*pi)
    return(list(r = sqrt(x^2+y^2), theta = theta))
  }
  #phi <- atan2(y = sqrt(x^2 + y^2), x = z)%%(2*pi)
  r <- sqrt(x^2 + y^2 + z^2)
  phi <- acos(z/r)
  theta <- atan2(y = y, x = x)%%(2*pi)
  return(list(r = r, theta = theta, phi = phi))
}

#' ConvertToPolar
#'
#' ConvertToPolar joins two tibbles with summary statistics by a common SNP id column,
#' with alleles already harmonized, for instance with AlleleFlip.
#' It then converts the original z-scores for trait 1 and 2 to polar coordinates, describing the angle theta with the x-axis (trait 1), and the distance r.
#' @param df1,df2 tibble, with at least columns snpid, a1, a2, beta, se, freq, pval
#' @param snpid Either a character indicating the common column name or,
#' in the case of two differently named columns,
#' a named character vector with the name of the common column in the first tibble as a name.
#' @param covsnps A list of snpid's, indicating which rows (SNPs) should be included
#' for calculation of the covariance matrix. High confidence SNPs are recommended.
#' @export
ConvertToPolar <- function(df1, df2, snpid, whiten = F, covsnps = c(), mahalanobis.threshold = 5){
  df <- dplyr::inner_join(df1, df2, by = snpid, suffix = c(".1", ".2"))
  if(nrow(df) == 0){return()}
  rm(df1)
  rm(df2)
  df %>%
    dplyr::select(starts_with("z")) %>%
    as.matrix() -> zmat
  if(whiten){
    if(length(covsnps) == 0){
      covsnps <- df$snpid
    }
    mahala <- mahalanobis(x = zmat[df$snpid %in% covsnps,], center = F, cov = cov(zmat[df$snpid %in% covsnps,]))
    S <- cov(zmat[df$snpid %in% covsnps,][mahala < mahalanobis.threshold^2,])
    zmat.white <- tcrossprod(zmat, whitening::whiteningMatrix(Sigma = S, method = "ZCA-cor"))
    rm(S)
  }else{zmat.white <- zmat; zmat <- NULL}
  polw <- PolarCoords(x = zmat.white[,1], y = zmat.white[,2])
  res <- dplyr::as_tibble(cbind(polw$r, polw$theta)) #both come from the whitened matrix now
  rm(polw)
  colnames(res) <- c("r", "theta")
  res <- res %>%
    dplyr::mutate(snpid = df$snpid) %>%
    dplyr::mutate(theta = theta)  %>%
    dplyr::mutate(theta.trans = (theta*4)%%(2*pi)) %>%
    dplyr::mutate(pval.1 = df$pval.1) %>%
    dplyr::mutate(pval.2 = df$pval.2) %>%
    dplyr::mutate(beta.1 = df$beta.1) %>%
    dplyr::mutate(beta.2 = df$beta.2) %>%
    dplyr::mutate(se.1 = df$se.1) %>%
    dplyr::mutate(se.2 = df$se.2) %>%
    dplyr::mutate(z.whitened.1 = zmat.white[,1]) %>%
    dplyr::mutate(z.whitened.2 = zmat.white[,2])
  res$theta.trans[res$theta.trans > pi] <- res$theta.trans[res$theta.trans > pi] - (2*pi)
  return(res)
}
