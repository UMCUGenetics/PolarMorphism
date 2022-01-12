#' PolarCoords
#'
#' Calculates the central angle between a SNP's position and the expected position under the null hypothesis of trait-specific effect.
#' The central angle is the angle between the vector from the origin of the p-dimensional sphere to the SNP, and the vector from the origin to the expected position.
#' Works per SNP (row-wise).
#' @param inputdf vector of effect sizes (z-scores) of 1 SNP on p traits
#' @returns a list with the distance r, and the angle in radians, normalized to describe a full circle
PolarCoords <- function(inputdf, debug = F){
  p <- ncol(inputdf)
  if(p < 2){
    print(paste("You supplied data for", p, "trait(s). Please supply data for 2 or more traits."))
    return()
  }
  inputdf <- abs(inputdf)
  rsquared <- apply(X = inputdf, MARGIN = 1, FUN = function(row){return(sum(row^2))})
  r <- sqrt(rsquared)
  if(p == 2 & !debug){
    angle <- atan2(y = inputdf[,2], x = inputdf[,1])%%(2*pi)
    angle <- (angle * 4)%%(2*pi)
    angle[angle > pi] <- angle[angle > pi] - (2*pi)
  }else{
    zerovec <- as.data.frame(matrix(0, nrow = nrow(inputdf), ncol = ncol(inputdf)))
    zerovec$maxdim <- apply(X = abs(inputdf), MARGIN = 1, FUN = which.max)
    zerovec$r <- r
    colnms <- colnames(zerovec)
    zerovec2 <- apply(zerovec, 1, function(row){row[row["maxdim"]] <- row["r"]; return(row)})
    zerovec <- as.data.frame(t(zerovec2))
    colnames(zerovec) <- colnms
    zerovec <- dplyr::select(.data = zerovec, c(-"maxdim", -"r"))
    angle <- apply(inputdf*zerovec, MARGIN = 1, sum)
    angle <- acos(angle/rsquared)
    correction.factor <- (2*pi)/acos(sqrt(p)/p)
    angle <- ((angle * correction.factor)%%(2*pi))/2
  }
  return(list(r = r, angle = angle))
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
ConvertToPolar <- function(dfnames, snpid, whiten = F, covsnps = c(), mahalanobis.threshold = 5, whitening.method = "ZCA-cor", LDcorrect = F, ld.path = "~/LDscores/eur_w_ld_chr/"){
  df <- get(dfnames[1], envir = .GlobalEnv)
  for(i in 2:length(dfnames)){
    df.next <- get(dfnames[i], envir = .GlobalEnv)
    df <- dplyr::inner_join(df, df.next, by = snpid, suffix = c(paste0(".", (i-1)), paste0(".", i)))
  }
  #rm(df.next)

  #df <- dplyr::rename(.data = df, "z.1" = "z")
  df <- df[complete.cases(df[,grep("^z", colnames(df))]),]
  if(nrow(df) == 0){return()}
  if(LDcorrect){df <- LDCorrect(df, ld.path = ld.path)}
  df %>%
    dplyr::select(starts_with("z")) %>%
    as.matrix() -> zmat
  if(whiten){
    if(length(covsnps) == 0){
      covsnps <- df$snpid
    }
    mahala <- mahalanobis(x = zmat[df$snpid %in% covsnps,], center = F, cov = cov(zmat[df$snpid %in% covsnps,]))
    print(paste("number of SNPs used for cov:", sum(mahala < mahalanobis.threshold^2)))
    S <- cov(zmat[df$snpid %in% covsnps,][mahala < mahalanobis.threshold^2,])
    zmat.white <- tcrossprod(zmat, whitening::whiteningMatrix(Sigma = S, method = whitening.method))
    rm(S)
  }else{zmat.white <- zmat; zmat <- NULL}
  polw <- PolarCoords(inputdf = zmat.white, debug = F)
  res <- dplyr::as_tibble(cbind(df[,snpid], polw$r, polw$angle, zmat.white)) #both come from the whitened matrix now
  rm(polw)
  colnames(res) <- c(snpid, "r", "angle", paste0("z.whitened.", 1:length(dfnames)))
  df %>%
    select(starts_with(c(snpid, "chr", "bp", "pos", "beta", "se", "pval", "a1", "a2"))) %>%
    dplyr::inner_join(res, by = snpid) -> res
  # res <- res %>%
  #   dplyr::mutate(snpid = df$snpid) %>%
  #   dplyr::mutate(angle = angle)  %>%
  #   #dplyr::mutate(angle.trans = (angle*4)%%(2*pi)) %>%
  #   dplyr::mutate(pval.1 = df$pval.1) %>%
  #   dplyr::mutate(pval.2 = df$pval.2) %>%
  #   dplyr::mutate(beta.1 = df$beta.1) %>%
  #   dplyr::mutate(beta.2 = df$beta.2) %>%
  #   dplyr::mutate(se.1 = df$se.1) %>%
  #   dplyr::mutate(se.2 = df$se.2) %>%
  #   dplyr::mutate(z.whitened.1 = zmat.white[,1]) %>%
  #   dplyr::mutate(z.whitened.2 = zmat.white[,2])
  #res$angle.trans[res$angle.trans > pi] <- res$angle.trans[res$angle.trans > pi] - (2*pi)
  return(res)
}
