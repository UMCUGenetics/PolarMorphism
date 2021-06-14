
#' PolarCoords
#'
#' Given two vectors for x- and y-coordinates, performs polar transformation. Returns the angle theta as a number between -pi and +pi.
#' @param x @param y vector of coordinates
PolarCoords <- function(x, y, z = NULL){
  if(is.null(z)){
    theta <- atan2(y = y, x = x)%%(2*pi)
    return(list(r = sqrt(x^2+y^2), theta = theta))
  }
  phi <- atan2(y = sqrt(x^2 + y^2), x = z)%%(2*pi)
  theta <- atan2(y = y, x = x)%%(2*pi)
  return(list(r = sqrt(x^2 + y^2 + z^2), theta = theta, phi = phi))
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
#' @param N A list of sample sizes for each trait
#' @param sample.prev In the case of binary traits: the sample prevalence. In the case of continuous traits: NA
#' @param population.prev In the case of binary traits: the population prevalence. In the case of continuous traits: NA
#' @param covsnps A list of logical vectors (with length equal to the number of rows of the dataframes), indicating which rows (SNPs) should be included
#' for calculation of the covariance matrix. High confidence SNPs are recommended.
#' @export
ConvertToPolar <- function(df1, df2, snpid, trait.names, whiten = F, ld.correct = F, ld.path = "~/LDscores/eur_w_ld_chr/", covsnps = c()){
  df <- dplyr::inner_join(df1, df2, by = snpid, suffix = c(".1", ".2"))
  if(nrow(df) == 0){return()}
  rm(df1)
  rm(df2)
  if(ld.correct){
    df <- df[df$snpid %in% unlist(readr::read_table2(paste0(ld.path, "w_hm3.snplist"), col_types = readr::cols_only("SNP" = "c"))),]
    ldscores <- readr::read_table2(paste0(ld.path, "1.l2.ldscore.gz"), col_names = T)[,c(2,6)]
    for(i in 2:22){
      ldscores <- rbind(ldscores, readr::read_table2(paste0(ld.path, i, ".l2.ldscore.gz"), col_names = T)[,c(2,6)])
    }
    df %>%
      dplyr::inner_join(ldscores, by = c("snpid" = "SNP")) -> df
    df$z.1 <- (abs(df$z.1) - (lm(data = df, formula = z.1 ~ L2)$coefficients[2] * df$L2)) * sign(df$z.1)
    df$z.2 <- (abs(df$z.2) - (lm(data = df, formula = z.2 ~ L2)$coefficients[2] * df$L2)) * sign(df$z.2)
  }
  df %>%
    dplyr::select(starts_with("z")) %>%
    as.matrix() -> zmat
  if(whiten){
    mahala <- mahalanobis(x = zmat[df$snpid %in% covsnps,], center = F, cov = cov(zmat[df$snpid %in% covsnps,]))
    S <- cov(zmat[df$snpid %in% covsnps,][mahala < 25,])
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
