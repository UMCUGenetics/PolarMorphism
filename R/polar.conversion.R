#' PolarCoords
#'
#' Calculates the central angle between a SNP's position and the expected position under the null hypothesis of trait-specific effect.
#' The central angle is the angle between the vector from the origin of the p-dimensional sphere to the SNP, and the vector from the origin to the expected position.
#' Works per SNP (row-wise).
#' @param inputdf vector of effect sizes (z-scores) of 1 SNP on p traits
#' @returns a list with the distance r, and the angle in radians, normalized to describe a full circle
#' @export
PolarCoords <- function(inputdf){
  inputdf <- abs(inputdf)
  rsquared <- apply(X = inputdf, MARGIN = 1, FUN = function(row){return(sum(row^2))})
  r <- sqrt(rsquared)
  zerovec <- as.data.frame(matrix(0, nrow = nrow(inputdf), ncol = ncol(inputdf)))
  zerovec$maxdim <- apply(X = inputdf, MARGIN = 1, FUN = which.max)
  zerovec$r <- r
  zerovec <- as.data.frame(apply(zerovec, 1, function(row){row[row["maxdim"]] <- row["r"]; return(row)}))
  zerovec <- dplyr::select(.data = zerovec, c(-"maxdim", -"r"))
  distsquared <- apply(X = (inputdf - zerovec)^2, MARGIN = 1, FUN = sum)
  angle <- acos((1 - (distsquared/(2*rsquared))))
  p <- ncol(inputdf)
  correction.factor <- (2*pi)/(acos(1 - ((p - sqrt(p))/(p))))
  angle <- angle * correction.factor
  return(list(r = r, angle = angle))
}

# PolarCoords <- function(x = NULL, y = NULL, z = NULL, inputvec = NULL){
#   if(is.null(inputved)){
#     theta <- atan2(y = y, x = x)%%(2*pi)
#     return(list(r = sqrt(x^2+y^2), theta = theta))
#   }
#
# }

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
ConvertToPolar <- function(dfnames, snpid, whiten = F, covsnps = c(), mahalanobis.threshold = 5){
  df <- get(dfnames[1], envir = .GlobalEnv)
  for(i in 2:length(dfnames)){
    df.next <- get(dfnames[i], envir = .GlobalEnv)
    df <- dplyr::inner_join(df, df.next, by = snpid, suffix = c("", paste0(".", i)))
  }
  df <- dplyr::rename(.data = df, "z.1" = "z")

  if(nrow(df) == 0){return()}
  rm(list = dfnames)
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
  polw <- PolarCoords(inputdf = zmat.white)
  res <- dplyr::as_tibble(cbind(polw$r, polw$angle)) #both come from the whitened matrix now
  rm(polw)
  colnames(res) <- c("r", "angle")
  res <- res %>%
    dplyr::mutate(snpid = df$snpid) %>%
    dplyr::mutate(angle = angle)  %>%
    dplyr::mutate(angle.trans = (angle*4)%%(2*pi)) %>%
    dplyr::mutate(pval.1 = df$pval.1) %>%
    dplyr::mutate(pval.2 = df$pval.2) %>%
    dplyr::mutate(beta.1 = df$beta.1) %>%
    dplyr::mutate(beta.2 = df$beta.2) %>%
    dplyr::mutate(se.1 = df$se.1) %>%
    dplyr::mutate(se.2 = df$se.2) %>%
    dplyr::mutate(z.whitened.1 = zmat.white[,1]) %>%
    dplyr::mutate(z.whitened.2 = zmat.white[,2])
  res$angle.trans[res$angle.trans > pi] <- res$angle.trans[res$angle.trans > pi] - (2*pi)
  return(res)
}
