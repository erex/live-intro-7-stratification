density.difference.ds <- function(dsobject, dsobject2) {
  # Purpose: compute difference in density estimates, significance test
  #           and confidence interval for the estimated difference
  #          This version based upon non-independence arising from
  #          estimation of a pooled detection function.
  # Input: dsobject - object created by `ds()` with strata and pooled detection function
  #                   assumes only two strata in the analysis
  #    optional: dsobject2 - if this exists, this means two strata do *not* share a common
  #                          detection, instead there are two distinct detection functions
  #                          represented by the two `ds` objects
  # Output: printed table of input values and table of difference,
  #         SE(difference), test statistic, P(test statistic), t-based CI
  # Reference: Section 3.6.5 of Buckland et al. (2001:84-86)
  satter <- function(n, k, p, cvn, cvP, cvEs) {
    # Eqn on p.121, Sect 4.8.1
    #  Inputs: 
    #     n - detections
    #     k - transects
    #     p - parameters in detection function
    #     cvn - encounter rate CV
    #     cvEs - mean group size CV
    #     cvP - cv in detection function parameters
    #  Output: 
    #    Satterthwaite degrees of freedom
    cvD <- sqrt(cvn^2 + cvP^2 + cvEs^2)
    satdf <- cvD^4 / (cvn^4/(k-1) + cvP^4/(n-p) + cvEs^4/(n-1))
    return(satdf)
  }
  #   merge stratum-specific values into vectors for easier calculation  
  secondstrat <- ifelse(nargs()==1, 2, 1)
  L <- dsobject$dht$clusters$summary[1:secondstrat,"Effort"]
  n <- dsobject$dht$clusters$summary[1:secondstrat,"n"]
  ncv <- dsobject$dht$clusters$summary[1:secondstrat, "se.ER"] / dsobject$dht$clusters$summary[1:secondstrat, "ER"]
  k <- dsobject$dht$clusters$summary[1:secondstrat,"k"]
  Es <- dsobject$dht$Expected.S[1:secondstrat, "Expected.S"]
  Escv <- dsobject$dht$Expected.S[1:secondstrat, "se.Expected.S"] / Es
  parm <- length(dsobject$ddf$par)
  f0 <- as.numeric(1/ (dsobject$ddf$fitted[1] * dsobject$ddf$meta.data$width))
  f0cv <- as.numeric(summary(dsobject)$ds$average.p.se / summary(dsobject)$ds$average.p)
  
  if(nargs()==1) {
    M <- n * Es / (2*L)  # Eqn 3.107
    D <- M * f0          # Eqn 3.106
    Dcv <- sqrt(ncv^2 + Escv^2 + f0cv^2)
    Ddiff <- (M[1] - M[2]) * f0  # Eqn 3.108
    Mvar <- M^2 * (ncv^2 + Escv^2)   # Eqn 3.111
    sumMvar <- sum(Mvar)             # Eqn 3.110
    Ddiffvar <- f0^2 * sumMvar + (M[1]-M[2]) * f0^2*f0cv^2  # Eqn 3.109
    satterth <- satter(n = n, k = k, p = parm, cvn = ncv, 
                       cvP = f0cv, cvEs = Escv)   # p. 121
  } else {
    L <- c(L, dsobject2$dht$clusters$summary[1,"Effort"])
    n <- c(n, dsobject2$dht$clusters$summary[1,"n"])
    ncv <- c(ncv, dsobject2$dht$clusters$summary[1, "se.ER"] / dsobject2$dht$clusters$summary[1, "ER"])
    k <- c(k, dsobject2$dht$clusters$summary[1,"k"])
    Es <- c(Es, dsobject2$dht$Expected.S[1, "Expected.S"])
    Escv <- c(Escv, dsobject2$dht$Expected.S[1, "se.Expected.S"] / Es[2])
    parm <- c(parm, length(dsobject2$ddf$par))
    D <- c(dsobject$dht$individuals$D$Estimate, dsobject2$dht$individuals$D$Estimate)
    Dcv <- c(dsobject$dht$individuals$D$cv, dsobject2$dht$individuals$D$cv)
    Ddiff <- D[1] - D[2]
    satterth <- c(dsobject$dht$clusters$D$df, dsobject2$dht$clusters$D$df)
    f0 <- c(f0, as.numeric(1/(dsobject2$dht$clusters$average.p*dsobject2$ddf$meta.data$width)))
    f0cv <- c(f0cv, as.numeric(summary(dsobject2)$ds$average.p.se / summary(dsobject2)$ds$average.p))
  }
  Dvar <- (D^2 * Dcv^2)
  if(nargs()!=1) Ddiffvar <- sum(Dvar)
  diffdf <- sum(Dvar)^2 / sum(Dvar^2/satterth)  # Eqn 3.101
  tstat <- Ddiff / sqrt(Ddiffvar)  # Eqn 3.100
  twotailedp <- pt(-tstat, diffdf) + pt(tstat, diffdf, lower.tail = FALSE)
  cibounds <- Ddiff + c(-1, +1) * qt(0.975, diffdf) * sqrt(Ddiffvar)
  bigtable <- data.frame(n.detect=n, cv.ER=ncv, line.length=L, n.transects=k,
                         group.size=Es, group.size.cv=Escv, nparm.detect=parm,
                         D.hat=D, D.hat.cv=Dcv, f0=f0, f0cv=f0cv)
  outstats <- data.frame(D1.minus.D2=Ddiff, SE.difference=sqrt(Ddiffvar),
                         t.statistic=tstat, df.t.stat=diffdf,
                         P.value=twotailedp, LCB=cibounds[1], UCB=cibounds[2])
  print(round(bigtable, 4))
  print(round(outstats, 4))
}