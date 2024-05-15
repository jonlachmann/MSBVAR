## SS.R:  R-side functions for the SS class objects that describe the
##        state-space for an MSBVAR model
##
## Date:  20071115 : Initial version
##        20081014 : Revised to use the bit package for storage.
##        20100617 : Updated to move sampler to gibbs.msbvar.  This
##                   file only now handles the posterior
##                   post-processing of the SS classed objects.
## Description:
##
## State spaces computed and stored in compressed
## form to reduce the computation and memory overhead involved in
## computing and storing state space draws.  R side storage of these
## state spaces is done using the bit package.  Functions here are
## used to summarize and plot the posterior state spaces.


######################################################################
# Summary functions for MS state-spaces from the BHLK filter or other
# estimators for the state-space
######################################################################

######################################################################
# Sums of the number of draws in each state as a function of the number
# of periods TT
######################################################################

sum.SS <- function(x, ...) {
  h <- x$h
  N <- length(x$ss.sample)
  ss <- x$ss.sample
  TTh <- virtual(x$ss.sample[[1]])$Length
  TT <- TTh / (h - 1)

  # Now find the sums
  sums <- apply(matrix(unlist(lapply(x$ss.sample, as.integer)), nrow = TTh, ncol = N), 1, sum)
  sums <- matrix(sums, TT, h - 1)
  sums <- cbind(sums, rep(N, TT) - rowSums(sums)) # Recover the h'th regime
  return(sums)
}


######################################################################
# Mean regime probability for each state
######################################################################
mean.SS <- function(x, ...) {
  sums <- sum.SS(x)
  N <- length(x$ss.sample)
  return(sums / N)
}

######################################################################
# Plot the mean posterior regime probabilities
######################################################################

plot.SS <- function(x, ylab = "State Probabilities", ...) {
  tmp <- mean.SS(x)
  shift <- x$p / attr(x, "freq")
  plot(
    ts(tmp,
      start = attr(x, "start") + shift,
      end = attr(x, "end"), frequency = attr(x, "freq")
    ),
    ###            deltat=attr(x, "freq")),
    plot.type = "single", col = 1:ncol(tmp),
    ylim = c(0, 1), ylab = ylab, ...
  )
  abline(h = 0.5, lty = 2, ...)
}
