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
    if(only.a2){flip <- sumstats$a2 != sumstats$a2.b}else{flip <- sumstats$a2 == sumstats$a1.b}
    if(!only.a2){sumstats$freq <- ifelse(flip, 1 - sumstats$freq, sumstats$freq)}
    sumstats$beta <- ifelse(flip, -sumstats$beta, sumstats$beta)
    sumstats$z <- sumstats$beta/sumstats$se
    sumstats$z[!same & !flip] <- NA
  }
  return(sumstats)
}

#' Import SNPs
#'
#' Imports a tibble with SNP meta information, and subsets to single nucleotide variants.
#' @param file A path for the SNP-file.
#' @param a1,a2 The name of the column for the reference and alternative alleles.
#' @export
ImportSNPs <- function(file = "~/3.Athero-integration/mtag/mafgt0.1.infogt0.95.ngt10000.freqSElt0.05.snps", a1 = "a1", a2 = "a2"){
  snps <- readr::read_table2(file, col_names = T)
  bases <- c("A", "C", "G", "T")
  snps <- snps %>%
    dplyr::filter((dplyr::pull(snps, var = a1) %in% bases) & (dplyr::pull(snps, var = a2) %in% bases))
}

