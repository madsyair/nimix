## samplers-afslice.R ----------------------------------------------------------
## Default sampler upgrade for the 3-4 parameter univariate families whose
## location and shape parameters are strongly correlated in the posterior.
##
## Measured motivation (n = 300, K = 2, 2500 iterations, this container):
##
##   family    baseline sampler   min ESS   min ESS/sec   worst correlation
##   fssn      univariate RW        12         0.28       cor(mu, alpha) -0.94
##   gmsnburr  univariate RW        46         0.86       cor(alpha, theta) 0.56
##
## Escalation followed the cheap-first ladder. RW_block on the correlated
## triplet made things WORSE (fssn min ESS 12 -> 7): the adaptive Gaussian
## proposal cannot follow the curved ridge within this budget. The automated
## factor slice sampler (AF_slice; Tibbits, Groendyke, Haran & Liechty 2014,
## JCGS 23(2), 543-563), which adapts a factor basis and slice-samples along
## it, is built for exactly this geometry:
##
##   fssn      AF_slice   min ESS 622   min ESS/sec 9.08   (x32 vs baseline)
##   gmsnburr  AF_slice   min ESS 417   min ESS/sec 4.90   (x5.7 vs baseline)
##
## with recovery unchanged (mu -4.23/4.14, alpha 0.62/1.85 vs truth -4/4,
## 0.5/2). The same mu-alpha coupling is structural to the whole
## Fernandez-Steel skew mechanism and to the Burr families, so the block is
## applied uniformly to the nine 3-4 parameter univariate specs. NUTS remains
## a possible further step (see the review study), but AF_slice already
## removes the pathology without new dependencies or AD requirements.

#' @include class-DistributionSpec.R
#' @include dist-msnburr.R
#' @include dist-gmsnburr.R
#' @include dist-fssn.R
#' @include dist-fossep.R
#' @include dist-fsst.R
#' @include dist-jfst.R
#' @include dist-sep.R
#' @include dist-lep.R
NULL

# Replace the per-parameter default samplers of each mixture component with a
# single AF_slice block over that component's (correlated) parameters.
# Engine-agnostic: works under FixedK (K components), DPM (L cluster slots),
# and MRF, because it discovers the *Tilde nodes from the model itself.
.componentAFSlice <- function(conf, model, params) {
  first <- paste0(params[1L], "Tilde")
  if (!(first %in% model$getVarNames())) return(invisible(conf))
  K <- length(model$expandNodeNames(first))
  for (k in seq_len(K)) {
    nds <- paste0(params, "Tilde[", k, "]")
    nds <- nds[nds %in% model$getNodeNames()]
    if (length(nds) < 2L) next
    for (nd in nds) conf$removeSamplers(nd)
    conf$addSampler(target = nds, type = "AF_slice")
  }
  invisible(conf)
}

#' @describeIn customizeSamplers AF_slice block over (mu, sigma, alpha).
#' @export
setMethod("customizeSamplers", "MSNBurrUvSpec",
  function(spec, conf, model, ...)
    .componentAFSlice(conf, model, c("mu", "sigma", "alpha")))

#' @describeIn customizeSamplers AF_slice block over (mu, sigma, alpha).
#' @export
setMethod("customizeSamplers", "MSNBurr2aUvSpec",
  function(spec, conf, model, ...)
    .componentAFSlice(conf, model, c("mu", "sigma", "alpha")))

#' @describeIn customizeSamplers AF_slice block over (mu, sigma, alpha, theta).
#' @export
setMethod("customizeSamplers", "GMSNBurrUvSpec",
  function(spec, conf, model, ...)
    .componentAFSlice(conf, model, c("mu", "sigma", "alpha", "theta")))

#' @describeIn customizeSamplers AF_slice block over (mu, sigma, alpha).
#' @export
setMethod("customizeSamplers", "FSSNUvSpec",
  function(spec, conf, model, ...)
    .componentAFSlice(conf, model, c("mu", "sigma", "alpha")))

#' @describeIn customizeSamplers AF_slice block over (mu, sigma, alpha, theta).
#' @export
setMethod("customizeSamplers", "FOSSEPUvSpec",
  function(spec, conf, model, ...)
    .componentAFSlice(conf, model, c("mu", "sigma", "alpha", "theta")))

#' @describeIn customizeSamplers AF_slice block over (mu, sigma, alpha, nu);
#'   the truncated nu node poses no problem for slice sampling.
#' @export
setMethod("customizeSamplers", "FSSTUvSpec",
  function(spec, conf, model, ...)
    .componentAFSlice(conf, model, c("mu", "sigma", "alpha", "nu")))

#' @describeIn customizeSamplers AF_slice block over (mu, sigma, alpha, theta).
#' @export
setMethod("customizeSamplers", "JFSTUvSpec",
  function(spec, conf, model, ...)
    .componentAFSlice(conf, model, c("mu", "sigma", "alpha", "theta")))

#' @describeIn customizeSamplers AF_slice block over (mu, sigma, nu).
#' @export
setMethod("customizeSamplers", "SEPUvSpec",
  function(spec, conf, model, ...)
    .componentAFSlice(conf, model, c("mu", "sigma", "nu")))

#' @describeIn customizeSamplers AF_slice block over (mu, sigma, nu).
#' @export
setMethod("customizeSamplers", "LEPUvSpec",
  function(spec, conf, model, ...)
    .componentAFSlice(conf, model, c("mu", "sigma", "nu")))
