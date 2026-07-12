## dist-fossep-core.R -----------------------------------------------------------
## FOSSEP d/p/q/r -- maintainer reference implementation (A. S. Choir),
## adopted verbatim. Batch B neo-normal family (Choir 2020).

#' @include neonorm-utils.R
NULL

#' Fernandez-Osiewalski-Steel Skew Exponential Power Distribution
#'
#' Density, distribution function, quantile function, and random generation
#' for the Fernandez-Osiewalski-Steel Skew Exponential Power (FOSSEP)
#' distribution.
#'
#' @name fossep-distribution
#'
#' @param x,q Numeric vector of quantiles.
#' @param p Numeric vector of probabilities.
#' @param n Integer. Number of observations to generate.
#' @param mu Numeric. Location parameter (default = 0).
#' @param sigma Numeric. Scale parameter, must be positive (default = 1).
#' @param alpha Numeric. Skewness parameter, must be positive (default = 2).
#'   Follows the Fernandez-Steel convention shared by all skew families in
#'   nimix: \code{alpha} is the FS skewness \eqn{\gamma}, so
#'   \eqn{P(X > \mu) = \alpha^2/(1 + \alpha^2)} and \code{alpha > 1} skews right.
#' @param theta Numeric. Kurtosis/exponential-power parameter, must be
#'   positive (default = 2).
#' @param log,log.p Logical. If TRUE, probabilities/densities are given as log.
#'   Default is FALSE.
#' @param lower.tail Logical. If TRUE (default), probabilities are
#'   \eqn{P[X \le x]}, otherwise \eqn{P[X > x]}.
#'
#' @return
#' \code{dfossep} gives the density, \code{pfossep} gives the distribution
#' function, \code{qfossep} gives the quantile function, and
#' \code{rfossep} generates random deviates.
#'
#' @details
#' The Fernandez-Osiewalski-Steel Skew Exponential Power distribution with
#' parameters \eqn{\mu}, \eqn{\sigma}, \eqn{\alpha}, and \eqn{\theta} has density:
#'
#' \deqn{f(y|\mu,\sigma,\alpha,\theta) = \frac{c}{\sigma}
#'   \exp\left(-\frac{1}{2} |\alpha z|^\theta\right) \quad\text{if } y < \mu}
#'
#' \deqn{f(y|\mu,\sigma,\alpha,\theta) = \frac{c}{\sigma}
#'   \exp\left(-\frac{1}{2} |z/\alpha|^\theta\right) \quad\text{if } y \ge \mu}
#'
#' where \eqn{z = (y - \mu)/\sigma},
#' \eqn{c = \frac{\alpha\theta}{(1+\alpha^2)\,2^{1/\theta}\,\Gamma(1/\theta)}}.
#'
#' When \eqn{\theta = 2}, it reduces to the Fernandez-Steel Skew Normal (FSSN).
#'
#' @references
#' Fernandez, C., Osiewalski, J., & Steel, M. F. (1995). Modeling and inference
#' with v-spherical distributions. Journal of the American Statistical
#' Association, 90(432), pp 1331-1340.
#'
#' Rigby, R. A., Stasinopoulos, M. D., Heller, G. Z., & De Bastiani, F. (2019).
#' Distributions for Modeling Location, Scale, and Shape: Using GAMLSS in R.
#' CRC Press.
#'
#' @examples
#' dfossep(0, mu = 0, sigma = 1, alpha = 2, theta = 2)
#'
#' x <- seq(-5, 5, by = 0.1)
#' plot(x, dfossep(x, alpha = 0.5, theta = 1.5), type = "l",
#'      ylab = "Density", main = "FOSSEP Densities")
#' lines(x, dfossep(x, alpha = 2, theta = 3), col = "red")
#'
#' pfossep(c(-1, 0, 1), mu = 0, sigma = 1, alpha = 2, theta = 2)
#' qfossep(c(0.025, 0.5, 0.975), mu = 0, sigma = 1, alpha = 2, theta = 2)
#'
#' set.seed(123)
#' r <- rfossep(1000, mu = 0, sigma = 1, alpha = 2, theta = 2)
#' hist(r, breaks = 30, freq = FALSE, main = "FOSSEP Random Samples")
#' curve(dfossep(x, alpha = 2, theta = 2), add = TRUE, col = "red", lwd = 2)
#'
#' @keywords distribution
#' @importFrom stats pgamma qgamma
#' @export
dfossep <- function(x, mu = 0, sigma = 1, alpha = 2, theta = 2, log = FALSE) {
    msgs <- .check_params_warn(
        list(mu = mu, sigma = sigma, alpha = alpha, theta = theta),
        c("sigma", "alpha", "theta"))
    if (length(msgs) > 0) return(rep(NaN, length(x)))

    a <- alpha
    b <- theta
    z <- (x - mu) / sigma

    loglik <- rep(NA_real_, length(x))
    left <- x < mu
    if (any(left, na.rm = TRUE)) {
        loglik[left] <- -0.5 * (a * abs(z[left]))^b
    }
    right <- x >= mu
    if (any(right, na.rm = TRUE)) {
        loglik[right] <- -0.5 * (abs(z[right]) / a)^b
    }

    loglik <- loglik - log(sigma) + log(a) - log1p(a^2) -
        (1 / b) * log(2) - lgamma(1 + 1 / b)

    if (log) loglik else exp(loglik)
}

#' @rdname fossep-distribution
#' @export
pfossep <- function(q, mu = 0, sigma = 1, alpha = 2, theta = 2,
                    lower.tail = TRUE, log.p = FALSE) {
    msgs <- .check_params_warn(
        list(mu = mu, sigma = sigma, alpha = alpha, theta = theta),
        c("sigma", "alpha", "theta"))
    if (length(msgs) > 0) return(rep(NaN, length(q)))

    a <- alpha
    b <- theta
    k <- a^2
    z <- (q - mu) / sigma
    z1 <- a * z / (2^(1 / b))
    z2 <- z / (a * (2^(1 / b)))

    cdf <- rep(NA_real_, length(q))
    left <- q < mu
    if (any(left, na.rm = TRUE)) {
        s1_left <- abs(z1[left])^b
        cdf[left] <- 1 - stats::pgamma(s1_left, 1 / b, 1)
    }
    right <- q >= mu
    if (any(right, na.rm = TRUE)) {
        s2_right <- abs(z2[right])^b
        cdf[right] <- 1 + k * stats::pgamma(s2_right, 1 / b, 1)
    }
    cdf <- cdf / (1 + k)

    if (!lower.tail) cdf <- 1 - cdf
    if (log.p) cdf <- log(cdf)

    cdf
}

#' @rdname fossep-distribution
#' @export
qfossep <- function(p, mu = 0, sigma = 1, alpha = 2, theta = 2,
                    lower.tail = TRUE, log.p = FALSE) {
    msgs <- .check_params_warn(
        list(mu = mu, sigma = sigma, alpha = alpha, theta = theta),
        c("sigma", "alpha", "theta"))
    if (length(msgs) > 0) return(rep(NaN, length(p)))

    pp <- if (log.p) exp(p) else p
    if (!lower.tail) pp <- 1 - pp

    if (any(pp < 0 | pp > 1, na.rm = TRUE)) {
        stop("p must be strictly between 0 and 1")
    }

    np <- length(pp)
    mu <- rep_len(mu, np); sigma <- rep_len(sigma, np)
    alpha <- rep_len(alpha, np); theta <- rep_len(theta, np)
    a <- alpha
    b <- theta
    k <- a^2

    eps <- .Machine$double.eps
    q <- rep(NA_real_, length(pp))
    left <- pp < (1 / (1 + k))
    if (any(left, na.rm = TRUE)) {
        q[left] <- mu[left] - (sigma[left] * (2^(1 / b[left])) / a[left]) *
            (stats::qgamma(pmax(pmin(1 - pp[left] * (1 + k[left]), 1 - eps), eps),
                           shape = 1 / b[left], scale = 1)^(1 / b[left]))
    }
    right <- pp >= (1 / (1 + k))
    if (any(right, na.rm = TRUE)) {
        q[right] <- mu[right] + sigma[right] * a[right] * (2^(1 / b[right])) *
            (stats::qgamma(pmax(pmin(
                (pp[right] * (1 + k[right]) - 1) / pmax(k[right], eps), 1 - eps), eps),
                shape = 1 / b[right], scale = 1)^(1 / b[right]))
    }

    q[pp == 0] <- -Inf
    q[pp == 1] <- Inf

    q
}

#' @rdname fossep-distribution
#' @export
rfossep <- function(n, mu = 0, sigma = 1, alpha = 2, theta = 2) {
    msgs <- .check_params_warn(
        list(mu = mu, sigma = sigma, alpha = alpha, theta = theta),
        c("sigma", "alpha", "theta"))
    if (length(msgs) > 0) return(rep(NaN, n))

    n <- ceiling(n)
    u <- stats::runif(n)
    qfossep(u, mu = mu, sigma = sigma, alpha = alpha, theta = theta)
}
