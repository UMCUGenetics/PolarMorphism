#' ConvertToPolar
#' 
#' ConvertToPolar joins two tibbles with summary statistics by a common SNP id column,
#' with alleles already harmonized, for instance with AlleleFlip.
#' @param df1 @param df2 tibble, with at least columns snpid, a1, a2, beta, se, freq, pval
#' @param snpid Either a character indicating the common column name or, 
#' in the case of two differently named columns, 
#' a named character vector with the name of the common column in the first tibble as a name.
#' @export
ConvertToPolar <- function(df1, df2, snpid){
  df <- inner_join(df1, df2, by = snpid, suffix = c(".1", ".2"))
  df <- df[complete.cases(df),]
  rm(df1)
  rm(df2)
  if(nrow(df) == 0){return()}
  df %>%
    dplyr::select(starts_with("z")) %>%
    as.matrix() -> zmat
  zmat.white <- whiten(X = zmat, center = F, method = "ZCA")
  pol <- CartToPol(x = zmat[,1], y = zmat[,2])
  polw <- CartToPol(x = zmat.white[,1], y = zmat.white[,2])
  res <- as_tibble(cbind(polw$r, pol$theta))
  rm(pol)
  rm(polw)
  colnames(res) <- c("r", "theta")
  res <- res %>%
    mutate(snpid = df$snpid) %>%
    mutate(theta = theta%%(2*pi))  %>%
    mutate(fold = abs(abs(theta - pi) - (0.5*pi))) %>%
    mutate(traitspecificity = (fold - 0.25*pi)/0.25*pi) %>%
    mutate(snpid = df$snpid) %>%
    mutate(pval.1 = df$pval.1) %>%
    mutate(pval.2 = df$pval.2) %>%
    mutate(beta.1 = df$beta.1) %>%
    mutate(beta.2 = df$beta.2) %>% 
    mutate(se.1 = df$se.1) %>%
    mutate(se.2 = df$se.2)
  return(res)
}