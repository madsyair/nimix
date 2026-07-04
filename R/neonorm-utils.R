## neonorm-utils.R -------------------------------------------------------------
## Numerically stable primitives for the neo-normal (MSNBurr) family.
## Reference implementation contributed by the package maintainer
## (A. S. Choir); see Iriawan (2000) and Choir (2020). Kept internal.

#' Numerically Stable log(1 + exp(x))
#'
#' Computes log(1 + exp(x)) avoiding overflow for large positive x and
#' underflow for large negative x.
#'
#' @param x Numeric vector.
#'
#' @return Numeric vector of the same length as `x`.
#'
#' @details
#' Piecewise approximation for numerical stability:
#' \itemize{
#'   \item \eqn{x \ge 33}: returns \eqn{x} (above this, \eqn{e^x} overflows
#'         double precision; \eqn{1 + e^x \approx e^x})
#'   \item \eqn{x \le -37}: returns \eqn{e^x} (below this,
#'         \eqn{\log(1 + \epsilon) \approx \epsilon} and \eqn{1 + e^x}
#'         evaluates to exactly 1 in double precision)
#'   \item otherwise: returns \code{\link[base]{log1p}(\link[base]{exp}(x))}
#' }
#'
#' @keywords internal
log1pexp <- function(x) {
    ifelse(x <= -37, exp(x),
           ifelse(x >= 33, x, log1p(exp(x))))
}

#' Numerically stable log(exp(a) + exp(b))
#'
#' Computes log(exp(a) + exp(b)) avoiding overflow. Equivalent to
#' \eqn{\max(a,b) + \log(1 + \exp(-|a - b|))}.
#'
#' @param a Numeric vector.
#' @param b Numeric vector.
#'
#' @return Numeric vector of the same length.
#'
#' @keywords internal
log_sum_exp <- function(a, b) {
    max_ab <- pmax(a, b)
    max_ab + log1p(exp(-abs(a - b)))
}

#' Numerically Stable log(1 - exp(x))
#'
#' Computes log(1 - exp(x)) for x < 0. Used in upper-tail probability
#' computations.
#'
#' @param x Numeric vector, must be < 0.
#'
#' @return Numeric vector of the same length as `x`.
#'
#' @keywords internal
log1mexp <- function(x) {
    if (any(x >= 0, na.rm = TRUE)) {
        stop("log1mexp requires x < 0")
    }
    result <- x
    mid <- (x >= -0.693) & (x < 0)
    low <- x < -0.693
    result[mid] <- log(-expm1(x[mid]))
    result[low] <- log1p(-exp(x[low]))
    result
}

#' Compute log-Omega for MSNBurr and MSNBurr-IIa
#'
#' Internal function for the log normalizing constant of the MSNBurr
#' and MSNBurr-IIa distributions.
#'
#' @param alpha Positive shape parameter.
#'
#' @return Numeric value of log(omega).
#'
#' @keywords internal
.log_omega_msnburr <- function(alpha) {
    LOG_SQRT_2PI <- 0.5 * log(2 * pi)
    limit <- .Machine$double.xmin * 10
    ifelse(alpha < limit,
           -(alpha + 1) * log(alpha) - LOG_SQRT_2PI,
           (alpha + 1) * log1p(1 / alpha) - LOG_SQRT_2PI)
}

#' Compute Log-Omega for GMSNBurr
#'
#' Internal function for the log normalizing constant of the GMSNBurr
#' distribution.
#'
#' @param alpha Positive shape parameter.
#' @param theta Positive shape parameter.
#'
#' @return Numeric value of log(omega).
#'
#' @keywords internal
.log_omega_gmsnburr <- function(alpha, theta) {
    log_beta <- lbeta(alpha, theta)
    LOG_SQRT_2PI <- 0.5 * log(2 * pi)
    limit <- .Machine$double.xmin * 10
    ifelse(alpha < limit,
           -LOG_SQRT_2PI + log_beta + alpha * (log(theta) - log(alpha)),
           -LOG_SQRT_2PI + log_beta -
               theta * (log(theta) - log(alpha)) +
               (alpha + theta) * log1p(theta / alpha))
}

#' Compute Omega for GMSNBurr
#'
#' Computes the normalizing constant omega for the GMSNBurr distribution.
#'
#' @param alpha Positive shape parameter.
#' @param theta Positive shape parameter.
#'
#' @return Numeric value of omega.
#'
#' @keywords internal
omega_gmsnburr <- function(alpha, theta) {
    if (anyNA(c(alpha, theta))) {
        stop("alpha and theta must not be NA")
    }
    if (any(alpha <= 0)) {
        stop("alpha must be positive")
    }
    if (any(theta <= 0)) {
        stop("theta must be positive")
    }
    exp(.log_omega_gmsnburr(alpha, theta))
}

#' Validate Distribution Parameters
#'
#' Internal function to validate parameters common to all neo-normal
#' distributions.
#'
#' @param mu Location parameter.
#' @param sigma Scale parameter (must be positive).
#' @param alpha Shape parameter (must be positive).
#' @param theta Shape parameter (must be positive).
#'
#' @return Invisible NULL if valid, otherwise stops with error.
#'
#' @keywords internal
.validate_params <- function(mu, sigma, alpha, theta) {
    if (anyNA(c(mu, sigma, alpha, theta))) {
        stop("Parameters must not contain NA values")
    }
    if (any(sigma <= 0)) {
        stop("sigma must be positive")
    }
    if (any(alpha <= 0)) {
        stop("alpha must be positive")
    }
    if (any(theta <= 0)) {
        stop("theta must be positive")
    }
    invisible(NULL)
}

MAX_EXP_ARG <- log(.Machine$double.xmax) - 10

.check_params_warn <- function(param_list, positive, .quiet = FALSE) {
    msgs <- character(0)
    if (anyNA(unlist(param_list))) {
        nm <- names(param_list)
        msgs <- c(msgs, paste(nm, "must not be NA", sep = " ", collapse = ", "))
    }
    for (pnm in positive) {
        v <- param_list[[pnm]]
        if (any(v <= 0, na.rm = TRUE)) {
            msgs <- c(msgs, paste0(pnm, " must be positive"))
        }
    }
    if (length(msgs) > 0 && !.quiet) {
        warning(paste(msgs, collapse = "; "))
    }
    msgs
}
