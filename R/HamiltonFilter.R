HamiltonFilter <- function(bigt, m, p, h, e, sig2, Qhat) {
  # Loop over regimes (1 to h) to calculate univariate Normal
  # or multivariate Normal density given parameter values
  ylik <- matrix(0, bigt - p, h)
  ypwlik <- array(0, c(bigt - p, h, h))
  for (iterh in seq_len(h)) {
    detsig2 <- abs(det(sig2[, , iterh]))
    invsig2 <- solve(sig2[, , iterh])
    for (itert in seq_len(bigt - p)) {
      tmpfit <- e[itert, , iterh]
      matmultwo <- tmpfit %*% invsig2 %*% tmpfit
      ylik[itert, iterh] <- exp(-(m / 2) * log(2 * pi) - (0.5 * log(detsig2)) - (0.5 * matmultwo))
    }
  }

  for (iterh1 in seq_len(h)) {
    for (iterh2 in seq_len(h)) {
      ypwlik[, iterh1, iterh2] <- Qhat[iterh1, iterh2] * ylik[, iterh2]
    }
  }

  ssAmat <- matrix(1, h + 1, h)
  ssEvec <- c(rep(0, h), 1)
  ssAmat[seq_len(h), ] <- diag(h) - t(Qhat)

  Invmatss <- solve(t(ssAmat) %*% ssAmat)

  pSt1_t1 <- Invmatss %*% t(ssAmat) %*% ssEvec

  f <- 0.0
  pytSt1St_t1_itert <- matrix(0, h, h)
  filtprSt1St <- array(0, c(bigt - p, h, h))

  for (its in seq_len(bigt - p)) {
    for (i in seq_len(h)) {
      for (j in seq_len(h)) {
        pytSt1St_t1_itert[i, j] <- pSt1_t1[i] * ypwlik[its, i, j]
      }
    }
    filt_llfval <- sum(pytSt1St_t1_itert)
    f <- f + log(filt_llfval)

    pytSt1St_t1_itert <- pytSt1St_t1_itert / filt_llfval

    filtprSt1St[its, , ] <- pytSt1St_t1_itert
    pSt1_t1 <- colSums(pytSt1St_t1_itert)
  }

  # integrate over St-1
  filtprSt <- rowSums(filtprSt1St, dim = 2)
  return(list(filtprSt = filtprSt, f = f, e = e))
}



ForwardFilter <- function (nvar, bigRk, bigK, bigT, nbeta, sigdraw, xidraw, llh) {
  llht <- matrix(0, bigT, bigK)
  lhma <- matrix(0, bigT, bigK)

  for (itk in seq_len(bigK)) {
    # Setup sigma matrices
    tmpsig <- sigdraw[, , itk]
    tmpsigld <- log(abs(det(tmpsig)))
    tmpsiginv <- solve(tmpsig)
    for (itt in seq_len(bigT)) {
      # tmpe = 1 by nvar, with residuals at time itt
      tmpe <- bigRk[itt, seq_len(nvar), itk]
      matmultwo <- tmpe %*% tmpsiginv %*% tmpe

      llht[itt, itk] <- -(nvar / 2) * log(2 * pi) - (1 / 2) * tmpsigld - (1 / 2) * matmultwo
    }
  }

  # max-adjusted llht for numerical stabilization
  summaxllh <- 0
  for (itt in seq_len(bigT)) {
    # set temporary max to first regime
    tmpmax <- max(llht[itt, ])
    summaxllh <- summaxllh + tmpmax
    for (itk in seq_len(bigK)) {
      lhma[itt, itk] <- exp(llht[itt, itk] - tmpmax)
    }
  }

  # Forward filtering
  sp <- 0
  pfilt <- matrix(0, bigT + 1, bigK)

  if (bigK > 1) {
    pfilt[1, ] <- rep(1 / bigK, bigK)

    transxi <- t(xidraw)

    for (itt in seq_len(bigT)) {
      ptmp <- transxi %*% pfilt[itt, ] * lhma[itt, ]

      # now to normalize, sum up and divide. these probabilities become
      # pfilt values at next iteration, t+1
      st <- sum(ptmp)
      transptmp <- t(ptmp)
      pfilt[itt + 1, ] <- transptmp / st

      sp <- sp + log(st)
    }
  }

  llh <- sp + summaxllh

  return(list(pfilt = pfilt, llh = llh))
}


