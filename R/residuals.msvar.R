#' Calculate weighted residuals over the different regimes based on their probability
#' @export
residuals.msvar <- function (x) {
  resids <- array(0, dim = dim(x$hreg$e)[-3])
  for (i in seq_len(dim(resids)[2])) {
    resids[, i] <- rowSums(x$hreg$e[, i, ] * x$fp)
  }
  return(resids)
}

#' Calculate fitted values
#' @export
fitted.msvar <- function (x) {
  residuals <- residuals.msvar(x)
  fitted <- x$init.model$Y + residuals
  return(fitted)
}