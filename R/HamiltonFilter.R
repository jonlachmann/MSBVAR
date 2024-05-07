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
  filt_llfval <- numeric(bigt - p)
  filtprSt1St <- array(0, c(bigt - p, h, h))

  for (its in seq_len(bigt - p)) {
    filt_llfval[its] <- 0
    for (i in seq_len(h)) {
      for (j in seq_len(h)) {
        pytSt1St_t1_itert[i, j] <- pSt1_t1[i] * ypwlik[its, i, j]
        filt_llfval[its] <- filt_llfval[its] + pytSt1St_t1_itert[i, j]
      }
    }
    f <- f + log(filt_llfval[its])

    pytSt1St_t1_itert <- pytSt1St_t1_itert / filt_llfval[its]

    filtprSt1St[its, , ] <- pytSt1St_t1_itert
    for (i in seq_len(h)) {
      pSt1_t1[i] <- 0
      for (j in seq_len(h)) {
        pSt1_t1[i] <- pSt1_t1[i] + pytSt1St_t1_itert[j, i]
      }
    }
  }

  # integrate over St-1
  filtprSt <- rowSums(filtprSt1St, dim = 2)
  return(list(filtprSt = filtprSt, f = f, e = e))
}
