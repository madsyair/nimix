## dist-fsst-core.R -----------------------------------------------------------
## FSST d/p/q/r -- maintainer reference implementation (A. S. Choir),
## adopted verbatim. Batch B neo-normal family (Choir 2020).

#' @include neonorm-utils.R
NULL

#' Fernandez-Steel Skew t Distribution
#'
#' Density, distribution function, quantile function, and random generation
#' for the Fernandez-Steel Skew t (FSST) distribution.
#'
#' @name fsst-distribution
#'
#' @param x,q Numeric vector of quantiles.
#' @param p Numeric vector of probabilities.
#' @param n Integer. Number of observations to generate.
#' @param mu Numeric. Location parameter (default = 0).
#' @param sigma Numeric. Scale parameter, must be positive (default = 1).
#' @param alpha Numeric. Skewness parameter, must be positive (default = 1,
#'   symmetric). Follows the Fernandez-Steel convention shared by all skew
#'   families in nimix: \code{alpha} is the FS skewness \eqn{\gamma}, so
#'   \eqn{P(X > \mu) = \alpha^2/(1 + \alpha^2)} and \code{alpha > 1} skews right.
#'   \eqn{\alpha = 1} gives the standard Student-t distribution with
#'   \eqn{\nu} degrees of freedom.
#' @param nu Numeric. Degrees of freedom, must be positive (default = 5).
#'   As \eqn{\nu \to \infty}, the FSST approaches the FSSN distribution.
#' @param log,log.p Logical. If TRUE, probabilities/densities are given as log.
#'   Default is FALSE.
#' @param lower.tail Logical. If TRUE (default), probabilities are
#'   \eqn{P[X \le x]}, otherwise \eqn{P[X > x]}.
#'
#' @return
#' \code{dfsst} gives the density, \code{pfsst} gives the distribution
#' function, \code{qfsst} gives the quantile function, and
#' \code{rfsst} generates random deviates.
#'
#' @details
#' The Fernandez-Steel Skew t distribution applies the Fernandez-Steel
#' skewing mechanism to a Student-t kernel with \eqn{\nu} degrees of freedom.
#'
#' \deqn{f(y|\mu,\sigma,\alpha,\nu) =
#'   \frac{2}{\sigma(\alpha + 1/\alpha)} t_\nu(z/\alpha) \quad\text{if } y < \mu}
#'
#' \deqn{f(y|\mu,\sigma,\alpha,\nu) =
#'   \frac{2}{\sigma(\alpha + 1/\alpha)} t_\nu(\alpha z) \quad\text{if } y \ge \mu}
#'
#' where \eqn{z = (y - \mu)/\sigma} and \eqn{t_\nu} is the standard Student-t
#' density with \eqn{\nu} degrees of freedom.
#'
#' The mean exists for \eqn{\nu > 1} and the variance exists for \eqn{\nu > 2}.
#'
#' @references
#' Fernandez, C., Osiewalski, J., & Steel, M. F. (1995). Modeling and inference
#' with v-spherical distributions. Journal of the American Statistical
#' Association, 90(432), pp 1331-1340.
#'
#' @examples
#' dfsst(0, mu = 0, sigma = 1, alpha = 1, nu = 5)
#'
#' x <- seq(-5, 5, by = 0.1)
#' plot(x, dfsst(x, alpha = 1, nu = 5), type = "l",
#'      ylab = "Density", main = "FSST Densities")
#' lines(x, dfsst(x, alpha = 2, nu = 3), col = "red")
#' lines(x, dfsst(x, alpha = 0.5, nu = 10), col = "blue")
#'
#' pfsst(c(-1, 0, 1), mu = 0, sigma = 1, alpha = 2, nu = 5)
#' qfsst(c(0.025, 0.5, 0.975), mu = 0, sigma = 1, alpha = 2, nu = 5)
#'
#' set.seed(123)
#' r <- rfsst(1000, mu = 0, sigma = 1, alpha = 2, nu = 5)
#' hist(r, breaks = 30, freq = FALSE, main = "FSST Random Samples")
#' curve(dfsst(x, alpha = 2, nu = 5), add = TRUE, col = "red", lwd = 2)
#'
#' @keywords distribution
#' @importFrom stats dt pt qt
#' @export
dfsst <- function(x, mu = 0, sigma = 1, alpha = 1, nu = 5, log = FALSE) {
    msgs <- .check_params_warn(
        list(mu = mu, sigma = sigma, alpha = alpha, nu = nu),
        c("sigma", "alpha", "nu"))
    if (length(msgs) > 0) return(rep(NaN, length(x)))

    ## Harmonised skewness convention (Ferreira & Steel 1998, 2007): the
    ## exported alpha IS the FS skewness gamma, so alpha > 1 skews right and
    ## P(X > mu) = alpha^2 / (1 + alpha^2). The maintainer's reference kernel
    ## below is written in the alpha = 1/gamma parameterisation, so we invert
    ## once, here, keeping d/p/q/r mutually consistent.
    alpha <- 1 / alpha

    a <- alpha
    nuf <- nu
    z <- (x - mu) / sigma

    log_const <- log(2) - log(sigma) - log(a + 1 / a)

    loglik <- rep(NA_real_, length(x))
    left <- x < mu
    if (any(left, na.rm = TRUE)) {
        loglik[left] <- log_const + dt(z[left] / a, df = nuf, log = TRUE)
    }
    if (any(!left, na.rm = TRUE)) {
        loglik[!left] <- log_const + dt(a * z[!left], df = nuf, log = TRUE)
    }

    if (log) loglik else exp(loglik)
}

.fsst_cutpoint <- function(alpha) {
    a2 <- alpha * alpha
    a2 / (1 + a2)
}

#' @rdname fsst-distribution
#' @export
pfsst <- function(q, mu = 0, sigma = 1, alpha = 1, nu = 5,
                  lower.tail = TRUE, log.p = FALSE) {
    msgs <- .check_params_warn(
        list(mu = mu, sigma = sigma, alpha = alpha, nu = nu),
        c("sigma", "alpha", "nu"))
    if (length(msgs) > 0) return(rep(NaN, length(q)))

    ## Harmonised skewness convention (Ferreira & Steel 1998, 2007): the
    ## exported alpha IS the FS skewness gamma, so alpha > 1 skews right and
    ## P(X > mu) = alpha^2 / (1 + alpha^2). The maintainer's reference kernel
    ## below is written in the alpha = 1/gamma parameterisation, so we invert
    ## once, here, keeping d/p/q/r mutually consistent.
    alpha <- 1 / alpha

    a <- alpha
    a2 <- a * a
    nuf <- nu
    z <- (q - mu) / sigma

    cdf <- rep(NA_real_, length(q))
    left <- q < mu
    if (any(left, na.rm = TRUE)) {
        cdf[left] <- (2 * a2 / (1 + a2)) * stats::pt(z[left] / a, df = nuf)
    }
    if (any(!left, na.rm = TRUE)) {
        cdf[!left] <- (a2 + 2 * stats::pt(a * z[!left], df = nuf) - 1) / (1 + a2)
    }

    if (!lower.tail) cdf <- 1 - cdf
    if (log.p) cdf <- log(cdf)

    cdf
}

#' @rdname fsst-distribution
#' @export
qfsst <- function(p, mu = 0, sigma = 1, alpha = 1, nu = 5,
                  lower.tail = TRUE, log.p = FALSE) {
    msgs <- .check_params_warn(
        list(mu = mu, sigma = sigma, alpha = alpha, nu = nu),
        c("sigma", "alpha", "nu"))
    if (length(msgs) > 0) return(rep(NaN, length(p)))

    ## Harmonised skewness convention (Ferreira & Steel 1998, 2007): the
    ## exported alpha IS the FS skewness gamma, so alpha > 1 skews right and
    ## P(X > mu) = alpha^2 / (1 + alpha^2). The maintainer's reference kernel
    ## below is written in the alpha = 1/gamma parameterisation, so we invert
    ## once, here, keeping d/p/q/r mutually consistent.
    alpha <- 1 / alpha

    pp <- if (log.p) exp(p) else p
    if (!lower.tail) pp <- 1 - pp

    if (any(pp < 0 | pp > 1, na.rm = TRUE)) {
        stop("p must be strictly between 0 and 1")
    }

    np <- length(pp)
    mu <- rep_len(mu, np); sigma <- rep_len(sigma, np)
    alpha <- rep_len(alpha, np); nu <- rep_len(nu, np)
    a <- alpha
    a2 <- a * a
    nuf <- nu
    cut <- a2 / (1 + a2)

    z <- numeric(length(pp))
    z[is.na(pp)] <- NA
    left <- pp <= cut & !is.na(pp)
    right <- !left & !is.na(pp)

    if (any(left, na.rm = TRUE)) {
        prob_left <- (1 + a2[left]) * pp[left] / (2 * a2[left])
        z[left] <- a[left] * stats::qt(prob_left, df = nuf[left])
    }
    if (any(right, na.rm = TRUE)) {
        prob_right <- ((1 + a2[right]) * pp[right] + 1 - a2[right]) / 2
        z[right] <- stats::qt(prob_right, df = nuf[right]) / a[right]
    }

    mu + sigma * z
}

#' @rdname fsst-distribution
#' @export
rfsst <- function(n, mu = 0, sigma = 1, alpha = 1, nu = 5) {
    msgs <- .check_params_warn(
        list(mu = mu, sigma = sigma, alpha = alpha, nu = nu),
        c("sigma", "alpha", "nu"))
    if (length(msgs) > 0) return(rep(NaN, n))


    ## alpha is passed straight through: qfsst() already applies the
    ## harmonised FS convention, so inverting here too would cancel it out.
    n <- ceiling(n)
    u <- stats::runif(n)
    qfsst(u, mu = mu, sigma = sigma, alpha = alpha, nu = nu)
}
