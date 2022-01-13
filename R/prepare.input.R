#' AlleleFlip
#'
#' Alleleflip joins a tibble with summary statistics and a tibble with SNP meta information,
#' by a common SNP id column, and aligns reference and alternative alleles to the meta tibble.
#' @param sumstats A tibble, with at least columns snpid, a1, a2, beta, se, freq, pval.
#' @param snps A tibble, with at least columns snpid, a1, a2.
#' @param snpid Either a character indicating the common column name or,
#' in the case of two differently named columns,
#' a named character vector with the name of the common column in the first tibble as a name.
#' @export
AlleleFlip <- function(sumstats, snps, snpid = "snpid", only.a2 = F){
  sumstats <- snps %>%
    dplyr::inner_join(sumstats, by = snpid, suffix = c("", ".b"))
  if(!only.a2){sumstats$a1.b <- toupper(sumstats$a1.b)}
  sumstats$a2.b <- toupper(sumstats$a2.b)
  if(identical(sumstats$a2, sumstats$a2.b)){
    sumstats <- sumstats %>%
      dplyr::select(-ends_with(".b")) %>%
      dplyr::mutate(z = beta/se)
  }else{
    same <- sumstats$a2 == sumstats$a2.b
    if(only.a2){flip <- sumstats$a2 != sumstats$a2.b
    }else{
      flip <- sumstats$a2 == sumstats$a1.b}
    if(!only.a2){sumstats$freq <- ifelse(flip, 1 - sumstats$freq, sumstats$freq)}
    sumstats$beta <- ifelse(flip, -sumstats$beta, sumstats$beta)
    sumstats$z <- sumstats$beta/sumstats$se
    sumstats <- sumstats[same | flip,]
  }
  return(sumstats)
}

#' LDCorrect
#'
#' LDCorrect performs LD correction given a data-frame or tibble containing at least snpid's and z-scores (column names starting with "z"),
#' and a path to a LD scores file (as formatted by the LDSC software)
#' @param df A tibble, with at least columns snpid and z-scores
#' @param ld.path the path to the LD scores files
#' @export
LDCorrect <- function(df, ld.path){
  #df <- df[df$snpid %in% unlist(readr::read_table2(paste0(ld.path, "w_hm3.snplist"), col_types = readr::cols_only("SNP" = "c"))),]
  ldscores <- readr::read_table2(paste0(ld.path, "1.l2.ldscore.gz"), col_names = T)[,c(2,6)]
  for(i in 2:22){
    ldscores <- rbind(ldscores, readr::read_table2(paste0(ld.path, i, ".l2.ldscore.gz"), col_names = T)[,c(2,6)])
  }
  df %>%
    dplyr::inner_join(ldscores, by = c("snpid" = "SNP")) -> df
  rm(ldscores)
  for(zcolnum in grep(pattern = "^z", colnames(df))){
    print(zcolnum)
    zvec <- unlist(df[,zcolnum])
    sign.zvec <- sign(zvec)
    fit <- lm(formula = abs(zvec) ~ df$L2)
    zvec <- (fit$residuals + fit$coefficients[1])
    df[,zcolnum] <- zvec * sign.zvec
  }
  rm(zvec)
  df <- dplyr::select(.data = df, -"L2")
  return(df)
}









