
#' AlleleFlip
#' 
#' Alleleflip() joins a tibble with summary statistics and a tibble with SNP meta information, by a common SNP id column, and aligns reference and alternative alleles to the meta tibble.
#' @param sumstats a tibble, with at least columns snpid, a1, a2, beta, se, freq
#' @param snps a tibble, with at least columns snpid, a1, a2
#' @param snpid either a character indicating the common column name or, in the case of two differently named columns, a named character vector with the name of the common column in the first tibble as a name
#' @export
AlleleFlip <- function(sumstats, snps, snpid = "snpid"){
  df <- snps %>%
    inner_join(df, by = by, suffix = c("", ".b"))
  df$a1.b <- toupper(df$a1.b); df$a2.b <- toupper(df$a2.b)
  if(identical(df$a1, df$a1.b)){
    df <- df %>%
      select(-ends_with(".b")) %>%
      mutate(z = beta/se)
  }else{
    same <- df$a1 == df$a1.b
    flip <- df$a1 == df$a2.b
    df$freq <- ifelse(flip, 1 - df$freq, df$freq)
    df$beta <- ifelse(flip, -df$beta, df$beta)
    df <- df %>% mutate(z = beta/se)
    df$z[!same & !flip] <- NA
  }
  return(df)
}