#' Create a forecast based on a msvar model
#' @param x The msvar model object
#' @param h The forecast horizon
#' @param samples The number of forecast samples to produce
#' @param Z_init The Z vector to use for prediction, defaults to the last part of Y if left as null
#' @return An array of forecast samples
#' @method predict MSVAR
#' @export
predict.MSVAR <- function (x, h, samples, Z_init = NULL) {
  state_vec <- seq_len(x$h)

  # Decompose sigma matrices for mvnorm draws
  sigmaU <- lapply(state_vec, function (i) chol(x$hreg$Sigmak[, , i]))

  # Get the Z vector if it is not provided
  if (is.null(Z_init)) {
    Z_init <- tail(embed(x$init.model$Y, x$p), 1)
  }

  # Draw states at time T
  states <- sample(state_vec, samples, replace = T, prob = tail(x$fp, 1))
  forecasts <- array(NA, dim = c(h, x$m, samples))
  for (i in seq_len(samples)) {
    # Simulate the path of states
    state1 <- sample(state_vec, 1, replace = T, prob = x$Q[states[i], ])
    path <- c(state1, rep(0, h - 1))
    for (j in seq_len(h)[-1]) {
      path[j] <- sample(state_vec, 1, replace = T, prob = x$Q[path[j - 1], ])
    }
    # Create the forecast
    fcst <- matrix(NA, h, x$m)

    # Get the initial Z vector
    Z <- Z_init
    for (j in seq_len(h)) {
      # Get the coefficient matrix to use
      A <- x$hreg$Bk[, , path[j]]
      fcst[j, ] <- t(A) %*% c(Z, 1)

      # Add error sample
      fcst[j, ] <- fcst[j, ] + t(sigmaU[[path[j]]]) %*% rnorm(x$m)

      # Modify Z
      Z <- c(fcst[j, ], Z[seq_len((x$m - 1) * x$p)])
    }
    # Store forecast
    forecasts[, , i] <- fcst
  }
  return(forecasts)
}

