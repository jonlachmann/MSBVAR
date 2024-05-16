#' Calculate weighted residuals over the different regimes based on their probability
#' @method residuals MSVAR
#' @export
residuals.MSVAR <- function (x) {
  fitted <- fitted(x)
  resids <- tail(x$init.model$y, nrow(fitted)) - fitted
  return(resids)
}

#' Calculate fitted values
#' @method fitted MSVAR
#' @export
fitted.MSVAR <- function (x) {
  Z <- tail(embed(x$init.model$y, x$p), nrow(x$fp))
  fitted <- matrix(0, nrow = nrow(Z), ncol = x$m)
  for (i in seq_len(x$h)) {
    # Get the coefficient matrix to use
    A <- x$hreg$Bk[, , i]
    fitted_i <- cbind(Z, 1) %*% A
    for (j in seq_len(ncol(fitted))) {
      fitted[, j] <- fitted[, j] + fitted_i[, j] * x$fp[, i]
    }
  }
  return(fitted)
}
