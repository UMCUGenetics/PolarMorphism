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

GetTermFunc <- function(kappa){
  result <- function(q, i){
    result2 <- (besselI(x = kappa, nu = i, expon.scaled = FALSE) * sin(i * q))/i
    return(result2)
  }
  return(result)
}

pvalangle <- function(relative.angle, kappa, tol = 1e-20, stepsize = 10){
  relative.angle <- relative.angle %% (2*pi)
  TermFuncKappa <- GetTermFunc(kappa = kappa)
  sum <- rep(0, length(relative.angle))
  flag <- TRUE
  iters <- 1:stepsize
  while(flag){
    terms <- mapply(FUN = TermFuncKappa, relative.angle, MoreArgs = list(iters))
    sum <- sum + apply(X = terms, MARGIN = 2, FUN = sum)
    if(all(abs(terms[nrow(terms),]) < tol)){
      flag <- FALSE
    }else{
      iters <- iters + stepsize
    }
  }
  res <- (relative.angle/(2*pi) + sum/(pi * besselI(x = kappa, nu = 0, expon.scaled = FALSE)) - 0.5) * 2
  res[res == -1] <- 0
  return(res)
}

# pvalangle <- function(relative.angle, kappa, tol){
#   vm.pval1 <- (CircStats::pvm(angle = -abs(relative.angle), mu = 0, kappa = kappa, acc = tol) - 0.5)*2
#   vm.pval2 <- PvonmisesRad.2(q = -abs(relative.angle), kappa = kappa, tol = tol)
#   return(vm.pval)
# }

#' PvalueForAngle
#'
#' PvalueForAngle performs p-value estimation using the von Mises distribution
#' @param angle.trans a vector of fourfold transformed angles
#' @param r a vector of distances, either expected distance (z.max) or actual distance
#' @export
PvalueForAngle <- function(angle.trans, r, tol = 1e-20){
  load("~/git/poly.gwas.integration/R/kappas.4foldtransform.Rda")
  kappas.list <- kappas.list2
  angle.pval <- rep(NA, length(angle.trans))
  #kappas.list <- kappas.list[11:401,]
  #kappas.list[1,"x.lo"] <- 0
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
  #angle.pval[intervals == 1] <- 1
  for(ival in ivals){
    kappa <- kappas.list$kappa[ival]
    ind <- intervals == ival
    angle.pval[ind] <- abs(pvalangle(relative.angle = angle.trans[ind], kappa = kappa, tol = tol))
  }
  return(angle.pval)
}



