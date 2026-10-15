# spict-plots.R - ggplot2/ggplotFL plots for SPiCT fits
# mpb/R/spict-plots.R
#
# The plotSpict* functions are thin wrappers: the fit is converted with
# spict2cmp() and plotted with the generic cmp* functions (cmp-plots.R),
# which also compare SPiCT with other assessment methods.

.spictNeed=function() {
  if (!requireNamespace("spict", quietly=TRUE))
    stop("package 'spict' is required", call.=FALSE)
  invisible(TRUE)}

.spictEnd=function(fit)
  floor(max(c(fit$inp$timeC, unlist(fit$inp$timeI)), na.rm=TRUE))

.spictAnnual=function(fit, nm, label=nm) {
  x  =spict::get.par(nm, fit, exp=TRUE)
  tm =if (nm=="logCpred") fit$inp$timeCpred else fit$inp$time
  sel=abs(tm-round(tm))<1e-8 & tm<=.spictEnd(fit)
  data.frame(qname=label, year=round(tm[sel]), data=x[sel, "est"],
             lower=x[sel, "ll"], upper=x[sel, "ul"], sd=x[sel, "sd"],
             row.names=NULL)}

#' Time series of a SPiCT fit as a data.frame
#'
#' Extracts annual (start of year) estimates and 95\% confidence intervals of
#' biomass, fishing mortality, catch and the MSY ratios from a SPiCT fit.
#'
#' @param fit a fitted SPiCT object (class \code{spictcls})
#' @param quants which quantities, any of \code{"B"}, \code{"F"},
#'   \code{"Catch"}, \code{"BBmsy"}, \code{"FFmsy"}
#'
#' @return a data.frame with columns \code{qname}, \code{year}, \code{data}
#'   (estimate), \code{lower}, \code{upper} and \code{sd} (log scale)
#' @export
spictTS=function(fit, quants=c("B", "F", "Catch", "BBmsy", "FFmsy")) {
  .spictNeed()
  map=c(B="logB", F="logF", Catch="logCpred", BBmsy="logBBmsy", FFmsy="logFFmsy")
  quants=match.arg(quants, names(map), several.ok=TRUE)
  res=do.call(rbind, lapply(quants, function(q) .spictAnnual(fit, map[[q]], q)))
  res$qname=factor(res$qname, levels=quants)
  res}

#' Coerce a SPiCT fit to FLQuants
#'
#' Returns annual estimates of biomass, F, catch, \eqn{B/B_{MSY}} and
#' \eqn{F/F_{MSY}} as \code{FLQuants}. If \code{iter} is given the
#' uncertainty is represented by independent log-normal draws in each year
#' (using the SPiCT standard errors on the log scale), so that
#' \code{ggplotFL::plot()} shows median and ribbons.
#'
#' @param fit a fitted SPiCT object
#' @param quants quantities to return, see \code{\link{spictTS}}
#' @param iter number of iterations; \code{NULL} returns the point estimates
#' @param seed random seed used when \code{iter} is set
#'
#' @return \code{FLQuants}
#' @export
spict2FLQuants=function(fit, quants=c("B", "F", "Catch", "BBmsy", "FFmsy"),
                        iter=NULL, seed=1234) {
  ts=spictTS(fit, quants)
  if (!is.null(iter)) set.seed(seed)
  res=lapply(split(ts, ts$qname, drop=TRUE), function(x) {
    dmns=list(year=as.character(x$year), iter=seq_len(max(1, iter)))
    if (is.null(iter))
      return(FLCore::FLQuant(x$data, dimnames=dmns[1]))
    val=exp(stats::rnorm(nrow(x)*iter, log(x$data), x$sd))
    FLCore::FLQuant(val, dimnames=dmns)})
  FLCore::FLQuants(res)}

#' Plots of a SPiCT fit
#'
#' Wrappers that convert a SPiCT fit with \code{\link{spict2cmp}} and call
#' the generic comparison plots: \code{plotSpictTS} (\code{\link{cmpTS}}),
#' \code{plotSpictStatus} (\code{\link{cmpStatus}}), \code{plotSpictKobe}
#' (\code{\link{cmpKobe}}), \code{plotSpictProduction}
#' (\code{\link{cmpProduction}}), \code{plotSpictFit} (\code{\link{cmpFit}}),
#' \code{plotSpictRuns} (\code{\link{cmpRuns}}, one-step-ahead residuals),
#' \code{plotSpictProcess} (\code{\link{cmpProcess}}), \code{plotSpictACF}
#' (\code{\link{cmpACF}}), \code{plotSpictQQ} (\code{\link{cmpQQ}}),
#' \code{plotSpictRetro} (\code{\link{cmpRetro}}), \code{plotSpictHindcast}
#' (\code{\link{cmpHindcast}}) and \code{plotSpictPriors}
#' (\code{\link{cmpPriors}}). The retrospective and hindcast use
#' \code{spict::retro}, run if \code{fit$retro} is missing. OSA and process
#' residuals are calculated if missing.
#'
#' @param fit a fitted SPiCT object
#' @param scale multiplier for biomass and catch
#' @param type \code{"resid"} (one-step-ahead) or \code{"process"}
#' @param nretroyear number of peels if \code{spict::retro} has to be run
#' @param ... passed to the generic plot
#'
#' @return a ggplot
#' @name plotSpict
NULL

.spictCmp=function(fit, scale=1, osa=FALSE, process=FALSE, retro=0) {
  .spictNeed()
  if (osa && is.null(fit$osar))              fit=spict::calc.osa.resid(fit)
  if (process && is.null(fit$process.resid)) fit=spict::calc.process.resid(fit)
  if (retro>0 && is.null(fit$retro))         fit=spict::retro(fit, nretroyear=retro)
  spict2cmp(fit, scale=scale, peels=retro>0)}

#' @rdname plotSpict
#' @export
plotSpictTS=function(fit, scale=1, ...) cmpTS(.spictCmp(fit, scale), ...)

#' @rdname plotSpict
#' @export
plotSpictStatus=function(fit, ...) cmpStatus(.spictCmp(fit), ...)

#' @rdname plotSpict
#' @export
plotSpictKobe=function(fit, ...) cmpKobe(.spictCmp(fit), ...)

#' @rdname plotSpict
#' @export
plotSpictProduction=function(fit, scale=1, ...) cmpProduction(.spictCmp(fit, scale), ...)

#' @rdname plotSpict
#' @export
plotSpictFit=function(fit, scale=1, ...) cmpFit(.spictCmp(fit, scale), ...)

#' @rdname plotSpict
#' @export
plotSpictRuns=function(fit, ...) cmpRuns(.spictCmp(fit, osa=TRUE), ...)

#' @rdname plotSpict
#' @export
plotSpictProcess=function(fit, ...) cmpProcess(.spictCmp(fit, process=TRUE), ...)

#' @rdname plotSpict
#' @export
plotSpictACF=function(fit, type=c("resid", "process"), ...)
  cmpACF(.spictCmp(fit, osa=TRUE, process=TRUE), type=match.arg(type), ...)

#' @rdname plotSpict
#' @export
plotSpictQQ=function(fit, type=c("resid", "process"), ...)
  cmpQQ(.spictCmp(fit, osa=TRUE, process=TRUE), type=match.arg(type), ...)

#' @rdname plotSpict
#' @export
plotSpictRetro=function(fit, nretroyear=5, ...) cmpRetro(.spictCmp(fit, retro=nretroyear), ...)

#' @rdname plotSpict
#' @export
plotSpictHindcast=function(fit, nretroyear=5, ...)
  cmpHindcast(.spictCmp(fit, retro=nretroyear), ...)

#' @rdname plotSpict
#' @export
plotSpictPriors=function(fit, ...) cmpPriors(.spictCmp(fit), ...)
