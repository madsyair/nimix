## dist-fssn-core.R -----------------------------------------------------------
## FSSN d/p/q/r -- maintainer reference implementation (A. S. Choir),
## adopted verbatim. Batch B neo-normal family (Choir 2020).

#' @include neonorm-utils.R
NULL

#' Fernandez-Steel Skew Normal Distribution
#'
#' Density, distribution function, quantile function, and random generation
#' for the Fernandez-Steel Skew Normal (FSSN) distribution.
#'
#' @name fssn-distribution
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
#'   \eqn{\alpha = 1} gives the standard normal distribution.
#' @param log,log.p Logical. If TRUE, probabilities/densities are given as log.
#'   Default is FALSE.
#' @param lower.tail Logical. If TRUE (default), probabilities are
#'   \eqn{P[X \le x]}, otherwise \eqn{P[X > x]}.
#'
#' @return
#' \code{dfssn} gives the density, \code{pfssn} gives the distribution
#' function, \code{qfssn} gives the quantile function, and
#' \code{rfssn} generates random deviates.
#'
#' @details
#' The Fernandez-Steel Skew Normal distribution is a special case of the
#' FOSSEP distribution with \eqn{\theta = 2}. It applies the Fernandez-Steel
#' skewing mechanism to a standard normal kernel.
#'
#' \deqn{f(y|\mu,\sigma,\alpha) = \frac{2}{\sigma(\alpha + 1/\alpha)}
#'   \phi(z/\alpha) \quad\text{if } y < \mu}
#'
#' \deqn{f(y|\mu,\sigma,\alpha) = \frac{2}{\sigma(\alpha + 1/\alpha)}
#'   \phi(\alpha z) \quad\text{if } y \ge \mu}
#'
#' where \eqn{z = (y - \mu)/\sigma} and \eqn{\phi} is the standard normal
#' density.
#'
#' @references
#' Fernandez, C., Osiewalski, J., & Steel, M. F. (1995). Modeling and inference
#' with v-spherical distributions. Journal of the American Statistical
#' Association, 90(432), pp 1331-1340.
#'
#' @examples
#' dfssn(0, mu = 0, sigma = 1, alpha = 1)
#'
#' x <- seq(-5, 5, by = 0.1)
#' plot(x, dfssn(x, alpha = 0.5), type = "l",
#'      ylab = "Density", main = "FSSN Densities")
#' lines(x, dfssn(x, alpha = 1), col = "red")
#' lines(x, dfssn(x, alpha = 2), col = "blue")
#'
#' pfssn(c(-1, 0, 1), mu = 0, sigma = 1, alpha = 2)
#' qfssn(c(0.025, 0.5, 0.975), mu = 0, sigma = 1, alpha = 2)
#'
#' set.seed(123)
#' r <- rfssn(1000, mu = 0, sigma = 1, alpha = 2)
#' hist(r, breaks = 30, freq = FALSE, main = "FSSN Random Samples")
#' curve(dfssn(x, alpha = 2), add = TRUE, col = "red", lwd = 2)
#'
#' @keywords distribution
#' @export
dfssn <- function(x, mu = 0, sigma = 1, alpha = 1, log = FALSE) {
    msgs <- .check_params_warn(
        list(mu = mu, sigma = sigma, alpha = alpha),
        c("sigma", "alpha"))
    if (length(msgs) > 0) return(rep(NaN, length(x)))

    ## Harmonised skewness convention (Ferreira & Steel 1998, 2007): the
    ## exported alpha IS the FS skewness gamma, so alpha > 1 skews right and
    ## P(X > mu) = alpha^2 / (1 + alpha^2). The maintainer's reference kernel
    ## below is written in the alpha = 1/gamma parameterisation, so we invert
    ## once, here, keeping d/p/q/r mutually consistent.
    alpha <- 1 / alpha

    a <- alpha
    z <- (x - mu) / sigma

    log_const <- log(2) - log(sigma) - log(a + 1 / a) -
        0.5 * log(2 * pi)

    loglik <- rep(NA_real_, length(x))
    left <- x < mu
    loglik[left] <- log_const - 0.5 * (z[left] / a)^2
    loglik[!left] <- log_const - 0.5 * (a * z[!left])^2

    if (log) loglik else exp(loglik)
}

.fssn_cutpoint <- function(alpha) {
    a2 <- alpha * alpha
    a2 / (1 + a2)
}

#' @rdname fssn-distribution
#' @importFrom stats pnorm
#' @export
pfssn <- function(q, mu = 0, sigma = 1, alpha = 1,
                  lower.tail = TRUE, log.p = FALSE) {
    msgs <- .check_params_warn(
        list(mu = mu, sigma = sigma, alpha = alpha),
        c("sigma", "alpha"))
    if (length(msgs) > 0) return(rep(NaN, length(q)))

    ## Harmonised skewness convention (Ferreira & Steel 1998, 2007): the
    ## exported alpha IS the FS skewness gamma, so alpha > 1 skews right and
    ## P(X > mu) = alpha^2 / (1 + alpha^2). The maintainer's reference kernel
    ## below is written in the alpha = 1/gamma parameterisation, so we invert
    ## once, here, keeping d/p/q/r mutually consistent.
    alpha <- 1 / alpha

    a <- alpha
    a2 <- a * a
    z <- (q - mu) / sigma
    cut <- a2 / (1 + a2)

    cdf <- rep(NA_real_, length(q))
    left <- q < mu
    if (any(left, na.rm = TRUE)) {
        cdf[left] <- (2 * a2 / (1 + a2)) * stats::pnorm(z[left] / a)
    }
    if (any(!left, na.rm = TRUE)) {
        cdf[!left] <- (a2 + 2 * stats::pnorm(a * z[!left]) - 1) / (1 + a2)
    }

    if (!lower.tail) cdf <- 1 - cdf
    if (log.p) cdf <- log(cdf)

    cdf
}

#' @rdname fssn-distribution
#' @importFrom stats qnorm
#' @export
qfssn <- function(p, mu = 0, sigma = 1, alpha = 1,
                  lower.tail = TRUE, log.p = FALSE) {
    msgs <- .check_params_warn(
        list(mu = mu, sigma = sigma, alpha = alpha),
        c("sigma", "alpha"))
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
    mu <- rep_len(mu, np); sigma <- rep_len(sigma, np); alpha <- rep_len(alpha, np)
    a <- alpha
    a2 <- a * a
    cut <- a2 / (1 + a2)

    z <- numeric(length(pp))
    z[is.na(pp)] <- NA
    left <- pp <= cut & !is.na(pp)
    right <- !left & !is.na(pp)

    if (any(left, na.rm = TRUE)) {
        prob_left <- (1 + a2[left]) * pp[left] / (2 * a2[left])
        z[left] <- a[left] * stats::qnorm(prob_left)
    }
    if (any(right, na.rm = TRUE)) {
        prob_right <- ((1 + a2[right]) * pp[right] + 1 - a2[right]) / 2
        z[right] <- stats::qnorm(prob_right) / a[right]
    }

    mu + sigma * z
}

#' @rdname fssn-distribution
#' @export
rfssn <- function(n, mu = 0, sigma = 1, alpha = 1) {
    msgs <- .check_params_warn(
        list(mu = mu, sigma = sigma, alpha = alpha),
        c("sigma", "alpha"))
    if (length(msgs) > 0) return(rep(NaN, n))


    ## alpha is passed straight through: qfssn() already applies the
    ## harmonised FS convention, so inverting here too would cancel it out.
    n <- ceiling(n)
    u <- stats::runif(n)
    qfssn(u, mu = mu, sigma = sigma, alpha = alpha)
}
