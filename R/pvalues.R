#' PvalueForR
#'
#' PvalueForR performs p-value estimation using the chi distribution
#' @param r a vector of distances
#' @param p degrees of freedom (the number of original GWAS)
#' @param exp.r a vector of expected distances
#' @export
PvalueForR <- function(r, p = 2, exp.r = NA){
  if(is.na(exp.r)){chi.pval <- pchisq(q = r^2, df = p, lower.tail = F)
  }else{chi.pval <- pchisq(q = r^2, df = p, ncp = exp.r^2, lower.tail = F)} # ifelse(r > exp.r, T, F)}
  return(chi.pval)
}

# TermFuncKappa is the modified bessel function of order i for the given kappa, where q is the angle to be tested, multiplied by sin(i*q)/i
GetTermFunc <- function(kappa){
  result <- function(q, i){
    result2 <- (besselI(x = kappa, nu = i, expon.scaled = FALSE) * sin(i * q))/i
    return(result2)
  }
  return(result)
}

pvalangle <- function(relative.angle, kappa, tol = 1e-20, stepsize = 10){
  # We need this function for the infinite sum that approaches the integral of the probability density function.
  TermFuncKappa <- GetTermFunc(kappa = kappa)
  # We initialize a vector of the same length as relative.angle, to store the sum results in
  sum <- rep(0, length(relative.angle))
  flag <- TRUE
  # We don't add 1 term but *stepsize* number of terms at a time to save compute time
  iters <- 1:stepsize
  while(flag){
    # terms is a stepsize x m matrix, where m is the number of SNPs (length of relative.angle)
    terms <- mapply(FUN = TermFuncKappa, relative.angle, MoreArgs = list(iters))
    # we sum each column of terms to get the sum per SNP
    sum <- sum + apply(X = terms, MARGIN = 2, FUN = sum)
    # we continue summing until the last term that we calculated is so small that it doesn't make much of a difference
    # we now wait until this is true for each SNP, maybe it could save time to only go on calculating for the SNPs that need it
    # but if we need to make new matrices each time that happens, it's probably not faster
    if(all(abs(terms[nrow(terms),]) < tol, na.rm = T)){
      flag <- FALSE
    }else{
      # we update the vector iters so it will calculate the *stepsize* next steps
      iters <- iters + stepsize
    }
  }
  # Once we have the sums, we do the calculation below to get the approximate value of the integral of the probability function for each SNP (i.e. the p-values)
  res <- ((relative.angle/(2*pi)) + sum/(pi * besselI(x = kappa, nu = 0, expon.scaled = FALSE)) + 0.5) * 2
  res[res == -1] <- 0
  return(res)
}

#' PvalueForAngle
#'
#' PvalueForAngle performs p-value estimation using the von Mises distribution
#' @param angle.trans a vector of fourfold transformed angles
#' @param r a vector of distances
#' @export
PvalueForAngle <- function(angle.trans, r, tol = 1e-20, kappa.file = "~/git/PolarMorphism/R/kappas.4foldtransform.Rda"){
  angle.trans <- -abs(angle.trans)
  load(kappa.file)
  kappas.list <- kappas.list2
  angle.pval <- rep(NA, length(angle.trans))
  kappas.list[1,"kappa"] <- 0
  ind <- c(seq(from = 0, to = 5, by = 0.05), seq(from = 5.10, to = 10, by = 0.10), seq(from = 10.50, to = 15, by = 0.50), seq(from = 16, to = 20, by = 1))
  kappas.list <- kappas.list[kappas.list$xmu %in% ind,]

  for(i in 2:nrow(kappas.list)){
    kappas.list[i,"x.lo"] <- mean(c(unlist(kappas.list[i-1,"xmu"]), unlist(kappas.list[i, "xmu"])))
  }
  kappas.list$x.hi <- c(kappas.list$x.lo[2:nrow(kappas.list)], kappas.list$xmu[nrow(kappas.list)])

  intervals <- findInterval(x = r, vec = kappas.list$x.lo)
  ivals <- unique(intervals)
  ivals <- ivals[order(ivals)]
  for(ival in ivals){
    kappa <- kappas.list$kappa[ival]
    ind <- intervals == ival
    angle.pval[ind] <- abs(pvalangle(relative.angle = angle.trans[ind], kappa = kappa, tol = tol))
  }
  return(angle.pval)
}


