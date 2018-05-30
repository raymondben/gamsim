#' Draw samples from a fitted GAM model
#'
#' @param object : a fitted object inheriting from "gam" or "scam" 
#' @param N positive integer: number of simulations to generate
#' @param newdata data.frame: a data frame giving the predictor variables from which to simulate
#' @param enforce string: enforce sample curves to be "decreasing", "increasing", or "none" for neither
#'
#' @return numeric matrix with N cols and as many rows as newdata
#'
#' @examples
#' x <- data.frame(t=seq(0, pi/2, length.out=101))
#' x$y <- sin(x$t)*25 + 5*rnorm(101)
#'
#' library(mgcv)
#'
#' ## fit GAM
#' fit <- gam(y~s(t), data=x)
#'
#' ## simulations from fitted object
#' draws <- get_draws(fit, N=5, newdata=x)
#'
#' plot(x$t, draws[,1], type="l")
#' for (k in 2:5) lines(x$t, draws[, k])
#'
#' @export
get_draws <- function(object, N, newdata, enforce="none") {
    assert_that(inherits(object, c("gam", "scam")))
    assert_that(is.string(enforce))
    enforce <- match.arg(tolower(enforce), c("decreasing", "increasing", "none"))
    do_get_draws(object=object, N=N, newdata=newdata, enforce=enforce, strict_N=TRUE)
}

do_get_draws <- function(object, N, newdata, enforce, strict_N) {
    ## evaluate basis function at new data values
    Cg <- predict(object, newdata=newdata, type="lpmatrix")
    ## random draws of model coefs
    sims <- rmvn(N,
                 mu=if (inherits(object,"scam")) as.numeric(object$coefficients.t) else coef(object),
                 sigma=if (inherits(object,"scam")) object$Vp.t else vcov(object))
    draws <- Cg %*% t(sims) ## N draws from the model posterior
    if (enforce %in% c("decreasing", "increasing")) {
        ## this seems a bit sledgehammer
        chk <- if (enforce=="decreasing") {
                   apply(diff(draws)>0,2,any) ## find draws that are not decreasing
               } else {
                   apply(diff(draws)>0,2,any) ## find draws that are not increasing
               }
        draws <- draws[, !chk]
        if (strict_N) {
            ## keep sampling until we get the number of draws that we requested
            while (ncol(draws)<N) {
                ##cat(sprintf("%s: got %d\n",names(object$model)[1],ncol(draws)))
                draws <- cbind(draws, do_get_draws(object=object, N=N, newdata=newdata, enforce=enforce, strict_N=FALSE))
            }
            draws <- draws[, seq_len(N)]
        }
    }
    unname(draws)
}


## multivariate normal random deviates
rmvn <- function(n, mu, sigma) {
    L <- mroot(sigma)
    m <- ncol(L)
    t(mu + L %*% matrix(rnorm(m*n), m, n))
}
