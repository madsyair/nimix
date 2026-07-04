## cluster-profile.R -------------------------------------------------------------
## Data-side cluster profiling: assign each observation to its MAP cluster and
## describe the observed data within each cluster. Complements summary(), which
## reports fitted component PARAMETERS; clusterProfile() reports what the data
## in each cluster actually looks like.

#' Profile the clusters of a fitted mixture
#'
#' Assigns each observation to its maximum a posteriori (MAP) cluster -- the
#' component it occupies most often across the retained MCMC draws -- and
#' summarises the observed data within each cluster: cluster size, proportion,
#' and per-variable within-cluster mean, standard deviation, and median. This
#' characterises what each cluster \emph{is} in terms of the data, complementing
#' \code{\link{summary}}, which reports the fitted component parameters, and
#' \code{plot(fit, type = "cluster")}, which shows the partition visually.
#'
#' Cluster ids are recoded 1, 2, ... by descending size (largest cluster
#' first). Because mixture labels are not identified, this data-side profile
#' does not require \code{\link{relabel}}: it summarises a partition, and a
#' partition is label-invariant.
#'
#' @param fit A \code{\linkS4class{FitResult}} from \code{\link{nimixClust}}.
#'   Regression fits are also accepted: the response and each covariate are
#'   profiled per regime.
#' @param variables Optional character vector selecting which data columns to
#'   profile (multivariate / regression); default profiles all.
#' @return A data.frame with one row per occupied cluster and columns
#'   \code{cluster}, \code{size}, \code{proportion}, followed by
#'   \code{<var>_mean}, \code{<var>_sd}, \code{<var>_median} for each variable.
#' @examples
#' \dontrun{
#' fit <- nimixClust(y, K = 3, method = "fixedk")
#' clusterProfile(fit)
#' }
#' @seealso \code{\link{summary}}, \code{plot(fit, type = "cluster")}
#' @export
clusterProfile <- function(fit, variables = NULL) {
  if (!methods::is(fit, "FitResult"))
    stop("clusterProfile() expects a FitResult.", call. = FALSE)

  alloc <- fit@clusterAllocation
  n <- ncol(alloc)

  # MAP cluster per observation (most frequent occupied component across draws)
  mapCl <- apply(alloc, 2L, function(col) {
    tb <- table(col)
    as.integer(names(tb)[which.max(tb)])
  })
  # recode to consecutive ids ordered by descending size
  sizes0 <- sort(table(mapCl), decreasing = TRUE)
  recode <- stats::setNames(seq_along(sizes0), names(sizes0))
  cl <- unname(recode[as.character(mapCl)])
  K <- length(sizes0)

  # assemble the data matrix to profile
  isReg <- isRegressionSpec(fit@distSpec)
  if (isReg) {
    X <- fit@prior$X
    resp <- as.numeric(fit@data)
    dat <- cbind(response = resp,
                 if (!is.null(X)) as.matrix(X)[, setdiff(
                   colnames(X), "(Intercept)"), drop = FALSE])
  } else if (.dataDimOf(fit@data) == 1L) {
    dat <- matrix(as.numeric(fit@data), ncol = 1L,
                  dimnames = list(NULL, "y"))
  } else {
    dat <- as.matrix(fit@data)
    if (is.null(colnames(dat)))
      colnames(dat) <- paste0("V", seq_len(ncol(dat)))
  }
  if (!is.null(variables)) {
    keep <- intersect(variables, colnames(dat))
    if (!length(keep))
      stop("None of 'variables' match the data columns: ",
           paste(colnames(dat), collapse = ", "), call. = FALSE)
    dat <- dat[, keep, drop = FALSE]
  }

  out <- data.frame(cluster = seq_len(K),
                    size = as.integer(sizes0),
                    proportion = as.numeric(sizes0) / n)
  for (v in colnames(dat)) {
    xs <- dat[, v]
    out[[paste0(v, "_mean")]]   <- vapply(seq_len(K),
      function(k) mean(xs[cl == k]), numeric(1))
    out[[paste0(v, "_sd")]]     <- vapply(seq_len(K),
      function(k) stats::sd(xs[cl == k]), numeric(1))
    out[[paste0(v, "_median")]] <- vapply(seq_len(K),
      function(k) stats::median(xs[cl == k]), numeric(1))
  }
  class(out) <- c("nimixProfile", "data.frame")
  attr(out, "mapCluster") <- cl
  out
}

#' @export
print.nimixProfile <- function(x, ...) {
  cat("Cluster profile (", nrow(x), " clusters, MAP allocation):\n", sep = "")
  df <- x; class(df) <- "data.frame"; attr(df, "mapCluster") <- NULL
  df$proportion <- round(df$proportion, 3)
  num <- vapply(df, is.numeric, logical(1))
  df[num] <- lapply(df[num], function(v) signif(v, 4))
  print(df, row.names = FALSE)
  invisible(x)
}
