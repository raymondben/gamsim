#' Draw samples from a fitted GAM model
#'
#' Note well: the "enforce" behaviour is thoroughly experimental. The constraint (decreasing/increasing/etc) is evaluated on the sampled output, so for this to make any sort of sense you need to provide \code{newdata} appropriately.
#'
#' @param object : a fitted object inheriting from "gam" or "scam"
#' @param N positive integer: number of simulations to generate
#' @param newdata data.frame: a data frame giving the predictor variables from which to simulate
#' @param enforce string: enforce sample curves to be "decreasing", "increasing", "nondecreasing", "nonincreasing", or "none" for neither
#'
#' @return numeric matrix with N cols and as many rows as newdata
#'
#' @examples
#' x <- data.frame(t = seq(0, pi/2, length.out = 101))
#' x$y <- sin(x$t)*25 + 5*rnorm(101)
#'
#' library(mgcv)
#'
#' ## fit GAM
#' fit <- gam(y ~ s(t), data = x)
#'
#' ## simulations from fitted object
#' draws <- get_draws(fit, N = 5, newdata = x)
#'
#' plot(x$t, draws[, 1], type = "l")
#' for (k in 2:5) lines(x$t, draws[, k])
#'
#' \dontrun{
#'   ## with constraint to be monotone increasing
#'   library(scam)
#'   fit <- scam(y ~ s(t, bs = "mpi"), data = x)
#'
#'   ## simulations from fitted object
#'   draws <- get_draws(fit, N = 5, newdata = x, enforce = "increasing")
#'
#'   plot(x$t, draws[, 1], type = "l")
#'   for (k in 2:5) lines(x$t, draws[, k])
#' }
#'
#' @export
get_draws <- function(object, N, newdata, enforce = "none") {
    assert_that(inherits(object, c("gam", "scam")))
    assert_that(is.string(enforce))
    enforce <- match.arg(tolower(enforce), c("decreasing", "increasing", "nondecreasing", "nonincreasing", "none"))
    do_get_draws(object = object, N = N, newdata = newdata, enforce = enforce, strict_N = TRUE)
}

do_get_draws <- function(object, N, newdata, enforce, strict_N, Cg) {
    if (missing(Cg)) {
        ## evaluate basis function at new data values
        Cg <- predict(object, newdata = newdata, type = "lpmatrix")
    }
    ## random draws of model coefs
    sims <- rmvn(N, mu = if (inherits(object,"scam")) as.numeric(object$coefficients.t) else coef(object),
                    sigma = if (inherits(object,"scam")) object$Vp.t else vcov(object))
    draws <- Cg %*% t(sims) ## N draws from the model posterior
    if (enforce %in% c("nondecreasing", "nonincreasing", "decreasing", "increasing")) {
        ## this seems a bit sledgehammer
        chk <- switch(enforce,
                      decreasing = apply(diff(draws) < 0, 2, all),
                      increasing = apply(diff(draws) > 0, 2, all),
                      nondecreasing = apply(diff(draws) >= 0, 2, all),
                      nonincreasing = apply(diff(draws) <= 0, 2, all),
                      stop("unexpected 'enforce' value"))
        draws <- draws[, chk]
        if (strict_N) {
            ## keep sampling until we get the number of draws that we requested
            while (ncol(draws) < N) {
                ##cat(sprintf("%s: got %d\n",names(object$model)[1],ncol(draws)))
                draws <- cbind(draws, do_get_draws(object = object, N = N, enforce = enforce, strict_N = FALSE, Cg = Cg))
            }
            draws <- draws[, seq_len(N)]
        }
    }
    unname(draws)
}


## testing code only
## enforce the constraint on the model parms, not the resulting curve
do_get_draws2 <- function(object, N, newdata, enforce, strict_N) {
    sims <- do_get_sims(object, N, enforce, strict_N = TRUE)
    Cg <- predict(object, newdata = newdata, type = "lpmatrix")
    draws <- Cg %*% t(sims) ## N draws from the model posterior
    unname(draws)
}

do_get_sims <- function(object, N, enforce, strict_N) {
    ## random draws of model coefs
    sims <- rmvn(N, mu = if (inherits(object,"scam")) as.numeric(object$coefficients.t) else coef(object),
                 sigma = if (inherits(object,"scam")) object$Vp.t else vcov(object))
    if (enforce %in% c("nondecreasing", "nonincreasing", "decreasing", "increasing")) {
        chk <- switch(enforce,
                      decreasing = apply(sims[, -1], 1, is_decreasing), ## the [-1] throws away the intercept, so this only works with one smooth term
                      increasing = apply(sims[, -1], 1, is_increasing),
                      nondecreasing = apply(sims[, -1], 1, is_nondecreasing),
                      nonincreasing = apply(sims[, -1], 1, is_nonincreasing),
                      stop("unexpected 'enforce' value"))
        sims <- sims[chk, ]
        if (strict_N) {
            ## keep sampling until we get the number of sims that we requested
            while (nrow(sims) < N) {
                ##cat(sprintf("%s: got %d\n", names(object$model)[1], nrow(sims)))
                sims <- rbind(sims, do_get_sims(object = object, N = N, enforce = enforce, strict_N = FALSE))
            }
            sims <- sims[seq_len(N), ]
        }
    }
    sims
}


## multivariate normal random deviates
rmvn <- function(n, mu, sigma) {
    L <- mroot(sigma)
    m <- ncol(L)
    t(mu + L %*% matrix(rnorm(m*n), m, n))
}

is_decreasing <- function(z) all(diff(z) < 0)
is_nonincreasing <- function(z) all(diff(z) <= 0)
is_increasing <- function(z) all(diff(z) > 0)
is_nondecreasing <- function(z) all(diff(z) >= 0)
