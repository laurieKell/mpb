# cmp-convert.R - convert assessment outputs to the common cmpAssess summary
# mpb/R/cmp-convert.R

.flq=function(x, yrs) FLCore::FLQuant(as.numeric(x), dimnames=list(year=as.character(yrs)))

.ptCurve=function(r, k, n, kMax=1.2, nX=201) {
  B=seq(0, kMax*k, length.out=nX)
  p=n-1
  data.frame(stock=B, yield=r/p*B*(1-(B/k)^p))}

# ---- SPiCT -------------------------------------------------------------------

#' Convert a SPiCT fit to a cmpAssess
#'
#' Uses only the exported SPiCT API (\code{spict::get.par}) and the
#' \code{inp}, \code{osar}, \code{process.resid} and \code{retro} elements of
#' the fit. Biomass and F are at the start of each year; \eqn{B/B_{MSY}} and
#' \eqn{F/F_{MSY}} are the SPiCT estimates (with their own CIs); the reference
#' points are the deterministic ones. The index fit is \eqn{qB_t}, extended
#' one year past the last observation for hindcasting. Retrospective peels
#' are taken from \code{fit$retro} (see \code{spict::retro}).
#'
#' @param fit a fitted SPiCT object (class \code{spictcls})
#' @param name method label
#' @param scale multiplier for biomass, catch and index (e.g. 1000 if the fit
#'   is in thousand tonnes)
#' @param process \code{"spict"} to use \code{fit$process.resid} (B and F),
#'   \code{"generic"} for residuals from the production function
#' @param peels convert retrospective peels if present
#'
#' @return a \code{cmpAssess}
#' @export
spict2cmp=function(fit, name="SPiCT", scale=1, process=c("spict", "generic"), peels=TRUE) {
  if (!requireNamespace("spict", quietly=TRUE)) stop("package 'spict' is required")
  process=match.arg(process)
  gp =function(nm) spict::get.par(nm, fit, exp=TRUE)
  ref=function(nm) unname(gp(nm)[1, "est"])
  inp=fit$inp
  end=floor(max(c(inp$timeC, unlist(inp$timeI)), na.rm=TRUE))
  sel=abs(inp$time-round(inp$time))<1e-8 & inp$time<=end
  yrs=round(inp$time[sel])
  ann=function(nm, s=1) {x=gp(nm)[sel, , drop=FALSE]; list(.flq(x[, "est"]*s, yrs),
                                                   .flq(x[, "ll"]*s, yrs), .flq(x[, "ul"]*s, yrs))}
  v  =list(stock=ann("logB", scale), harvest=ann("logF"),
           bbmsy=ann("logBBmsy"), ffmsy=ann("logFFmsy"))
  ts =lapply(v, `[[`, 1); lo=lapply(v, `[[`, 2); up=lapply(v, `[[`, 3)
  ts$catch=.flq(inp$obsC*scale, inp$timeC)

  K=ref("logK")*scale; n=ref("logn"); r=ref("logr")
  rp=c(bmsy=ref("logBmsyd")*scale, fmsy=ref("logFmsyd"), msy=ref("logMSYd")*scale,
       k=K, r=r, n=n)

  B =gp("logB"); q=gp("logq")[, "est"]
  tB=round(inp$time, 6)
  idx=lapply(seq_along(inp$obsI), function(i) {
    ti=c(inp$timeI[[i]], end+1)
    j =match(round(ti, 6), tB)
    data.frame(series=paste("Index", i), year=ti,
               obs=c(inp$obsI[[i]], NA)*scale, hat=q[i]*B[j, "est"]*scale,
               lower=q[i]*B[j, "ll"]*scale, upper=q[i]*B[j, "ul"]*scale)})
  cp=gp("logCpred")
  ci=match(inp$timeC, inp$timeCpred)
  idx=c(list(data.frame(series="Catch", year=inp$timeC, obs=inp$obsC*scale,
                        hat=cp[ci, "est"]*scale, lower=cp[ci, "ll"]*scale, upper=cp[ci, "ul"]*scale)), idx)
  idx=do.call(rbind, idx)
  idx=idx[!is.na(idx$hat) | !is.na(idx$obs), ]

  rsd=NULL
  if (!is.null(fit$osar)) {
    o=fit$osar
    rsd=do.call(rbind, c(list(data.frame(series="Catch", year=o$timeC, residual=o$logCpres)),
      lapply(seq_along(o$logIpres), function(i)
        data.frame(series=paste("Index", i), year=o$timeI[[i]], residual=o$logIpres[[i]]))))}

  prc=NULL
  if (process=="spict" && !is.null(fit$process.resid)) {
    pr=fit$process.resid
    prc=rbind(data.frame(series="Process B", year=pr[[1]], residual=pr[[2]]),
              data.frame(series="Process F", year=pr[[1]], residual=pr[[3]]))}

  prior=.spictPriorDensity(fit, scale)

  pl=NULL
  if (peels && length(fit$retro)) {
    pl=lapply(fit$retro, function(f) {
      if (inherits(f, "try-error") || is.null(f$inp)) return(NULL)
      te=floor(max(c(f$inp$timeC, unlist(f$inp$timeI)), na.rm=TRUE))
      if (te>=end) return(NULL)
      spict2cmp(f, name=name, scale=scale, process="generic", peels=FALSE)})
    pl=pl[!vapply(pl, is.null, logical(1))]
    names(pl)=vapply(pl, function(p) as.character(max(p$index$year[is.finite(p$index$obs)])), character(1))}

  cmpAssess(name, FLCore::FLQuants(ts), FLCore::FLQuants(lo), FLCore::FLQuants(up), rp,
            index=idx, resid=rsd, process=prc, curve=.ptCurve(r, K, n), prior=prior, peels=pl)}

.spictPriorDensity=function(fit, scale=1) {
  pr =fit$inp$priors
  nms=names(pr)[fit$inp$priorsuseflags==1]
  lab=c(logr="r", logK="k", logn="n", logbkfrac="psi", logsdi="sdi",
        logalpha="alpha", logbeta="beta", logm="msy", logq="q")
  res=lapply(nms, function(nm) {
    p=pr[[nm]]
    if (is.list(p)) p=p[[1]]
    po=tryCatch(spict::get.par(nm, fit)[1, c("est", "sd")], error=function(e) NULL)
    if (is.null(po) || !all(is.finite(po))) return(NULL)
    s =if (nm %in% c("logK", "logm")) log(scale) else 0
    mu=c(p[1], po[["est"]])+s
    sd=c(p[2], po[["sd"]])
    rg=range(stats::qnorm(c(0.001, 0.999), mu[1], sd[1]), stats::qnorm(c(0.001, 0.999), mu[2], sd[2]))
    x =seq(rg[1], rg[2], length.out=300)
    l =if (nm %in% names(lab)) lab[[nm]] else sub("^log", "", nm)
    rbind(data.frame(parameter=l, what="Prior",     x=exp(x), density=stats::dnorm(x, mu[1], sd[1])),
          data.frame(parameter=l, what="Posterior", x=exp(x), density=stats::dnorm(x, mu[2], sd[2])))})
  do.call(rbind, res)}

# ---- FLStock + FLBRP ---------------------------------------------------------

#' Convert an age-structured assessment (FLStock and FLBRP) to a cmpAssess
#'
#' Biomass is spawning stock or exploitable biomass, F is either the harvest
#' rate (catch/biomass, comparable with biomass dynamic models) or
#' \eqn{\bar F}, and the reference points and equilibrium yield curve come
#' from the \code{FLBRP}. \eqn{r} and \eqn{n} are the Pella-Tomlinson
#' parameters with the same MSY, \eqn{B_{MSY}} and \eqn{K} (if
#' \code{FLRebuild} is available). Iters in the \code{FLStock} give CIs.
#'
#' @param stk an \code{FLStock}
#' @param brp an \code{FLBRP} for the same stock
#' @param name method label
#' @param biomass \code{"eb"} (exploitable biomass, needs \code{FLRebuild})
#'   or \code{"ssb"}
#' @param harvest \code{"hr"} for catch/biomass or \code{"fbar"}
#' @param index optional \code{FLQuant} of observations; its fit is
#'   \eqn{qB_t} with \eqn{q} the geometric mean of index/biomass
#'
#' @return a \code{cmpAssess}
#' @export
flbrp2cmp=function(stk, brp, name="FLBRP", biomass=c("eb", "ssb"),
                   harvest=c("hr", "fbar"), index=NULL) {
  biomass=match.arg(biomass)
  harvest=match.arg(harvest)
  eb =biomass=="eb"
  if (eb && !requireNamespace("FLRebuild", quietly=TRUE)) stop("biomass='eb' requires FLRebuild")
  B  =if (eb) FLRebuild::ebiomass(stk) else FLCore::ssb(stk)
  C  =FLCore::catch(stk)
  H  =if (harvest=="hr") C/B else FLCore::fbar(stk)
  rf =if (eb) FLRebuild::refptsEB(brp) else FLBRP::refpts(brp)
  qn =if (eb) "eb" else "ssb"
  bmsy=c(rf["msy", qn])[1]; msy=c(rf["msy", "yield"])[1]; k=c(rf["virgin", qn])[1]
  fmsy=if (harvest=="hr") msy/bmsy else c(rf["msy", "harvest"])[1]
  rp =c(bmsy=bmsy, fmsy=fmsy, msy=msy, k=k)
  if (requireNamespace("FLRebuild", quietly=TRUE)) {
    pt=tryCatch(FLRebuild::pellatParams(FLCore::FLPar(msy=msy, bmsy=bmsy, k=k)), error=function(e) NULL)
    if (!is.null(pt)) rp=c(rp, r=c(pt["r"])[1], n=c(pt["p"])[1]+1)}

  eq =brp
  fc =c(FLBRP::refpts(brp)["crash", "harvest"])[1]
  if (!is.finite(fc)) fc=4*c(FLBRP::refpts(brp)["msy", "harvest"])[1]
  FLBRP::fbar(eq)=FLCore::FLQuant(seq(0, 1, length.out=101)*fc)
  crv=data.frame(stock=c(if (eb) FLRebuild::ebiomass(eq) else FLCore::ssb(eq)),
                 yield=c(FLCore::catch(eq)))

  ts =FLCore::FLQuants(stock=B, harvest=H, catch=C)
  idx=NULL
  if (!is.null(index)) {
    yi =intersect(dimnames(index)$year, dimnames(B)$year)
    q  =exp(mean(log(c(index[, yi])/c(B[, yi])), na.rm=TRUE))
    idx=data.frame(series="Index", year=as.numeric(yi), obs=c(index[, yi]), hat=q*c(B[, yi]))}
  cmpAssess(name, ts, refpts=rp, index=idx, curve=crv)}

# ---- biodyn ------------------------------------------------------------------

#' Convert a biodyn to a cmpAssess
#'
#' Biomass, harvest rate and catch from the \code{biodyn}, reference points
#' from \code{mpb::refpts} and the Pella-Tomlinson parameters, index fits
#' and residuals from the \code{diags} slot. The fitted index is
#' \eqn{q(B_t+B_{t+1})/2} for every catch year, so a peel whose catch has
#' been projected one year beyond the last index observation provides a
#' one-step-ahead prediction for hindcasting. JABBA fits can be converted
#' with \code{\link{jabba2cmp}}, or via \code{jabba2biodyn}.
#'
#' @param object a \code{biodyn}
#' @param name method label
#' @param peels optional named list of \code{biodyn} fitted to truncated data
#'
#' @return a \code{cmpAssess}
#' @export
biodyn2cmp=function(object, name="biodyn", peels=NULL) {
  par=object@params
  pv =function(nm) if (nm %in% dimnames(par)$params) stats::median(c(par[nm, ])) else NA
  r=pv("r"); k=pv("k"); p=pv("p")
  rf=mpb::refpts(object)
  rp=c(bmsy=stats::median(c(rf["bmsy", ])), fmsy=stats::median(c(rf["fmsy", ])),
       msy=stats::median(c(rf["msy", ])), k=k, r=r, n=p+1)
  C  =object@catch
  cy =dimnames(C)$year
  S  =object@stock
  d  =object@diags
  idx=rsd=NULL
  last=max(as.numeric(cy))
  if (is.data.frame(d) && nrow(d)>0) {
    last=max(d$year[is.finite(d$obs)])
    sy =as.numeric(dimnames(S)$year)
    s  =apply(S@.Data, 2, stats::median)
    nms=unique(d$name)
    idx=do.call(rbind, lapply(seq_along(nms), function(i) {
      q  =pv(paste0("q", i))
      yr =as.numeric(cy)
      mid=(s[match(yr, sy)]+s[match(yr+1, sy)])/2
      di =d[d$name==nms[i], ]
      data.frame(series=paste("Index", i), year=yr,
                 obs=di$obs[match(yr, di$year)], hat=q*mid)}))
    rsd=data.frame(series=paste("Index", match(d$name, nms)), year=d$year, residual=d$residual)}
  yrs=cy[as.numeric(cy)<=last]
  ts =FLCore::FLQuants(stock=S[, yrs], harvest=(C/S[, cy])[, yrs], catch=C[, yrs])
  pl =NULL
  if (length(peels)) {
    pl=lapply(peels, biodyn2cmp, name=name)
    names(pl)=vapply(pl, function(z) as.character(.cmpTerminal(z)), character(1))}
  cmpAssess(name, ts, refpts=rp, index=idx, resid=rsd,
            curve=if (all(is.finite(c(r, k, p)))) .ptCurve(r, k, p+1), peels=pl)}

# ---- JABBA -------------------------------------------------------------------

#' Convert a JABBA fit to a cmpAssess
#'
#' Uses the \code{timeseries} (posterior median and 95\% CIs of B, F,
#' \eqn{B/B_{MSY}} and \eqn{F/F_{MSY}}), \code{estimates}, \code{diags}
#' (index fits and residuals), \code{pfunc} (production function),
#' \code{timeseries[, , "procB"]} (process deviations) and
#' \code{pars_posterior} elements of a \code{JABBA::fit_jabba} result. Priors
#' for \eqn{r}, \eqn{K} and \eqn{\psi} are taken from the JABBA input when
#' available (log-normal, as \code{c(mean, log sd)}).
#'
#' @param object a JABBA fit, or a list with elements \code{fit} and
#'   \code{input} as returned by \code{FLRebuild::runJABBA}
#' @param name method label
#'
#' @return a \code{cmpAssess}
#' @export
jabba2cmp=function(object, name="JABBA") {
  inp=NULL
  if (!is.null(object$fit)) {inp=object$input; object=object$fit}
  tsr=object$timeseries
  yrs=dimnames(tsr)[[1]]
  get=function(v, s) .flq(tsr[, s, v], yrs)
  qs =c(stock="B", harvest="F", bbmsy="BBmsy", ffmsy="FFmsy")
  ts =lapply(qs, get, s="mu"); lo=lapply(qs, get, s="lci"); up=lapply(qs, get, s="uci")
  ts$catch=.flq(object$catch, object$yr)
  est=object$estimates
  e  =function(nm) if (nm %in% rownames(est)) est[nm, "mu"] else NA
  rp =c(bmsy=e("SBmsy"), fmsy=e("Hmsy"), msy=e("MSY"), k=e("K"), r=e("r"), n=e("m"))
  crv=NULL
  if (!is.null(object$pfunc)) {
    pf =object$pfunc
    crv=data.frame(stock=pf$SB_i, yield=pf$SP)
    crv=crv[!duplicated(crv$stock), ]}
  d  =object$diags
  sr =paste("Index", match(d$name, unique(d$name)))
  idx=data.frame(series=sr, year=d$year, obs=d$obs, hat=d$hat,
                 lower=d$hat.lci, upper=d$hat.uci)
  rsd=data.frame(series=sr, year=d$year, residual=d$residual)
  prc=if ("procB" %in% dimnames(tsr)[[3]])
    data.frame(series="Process B", year=as.numeric(yrs), residual=tsr[, "mu", "procB"])
  prior=NULL
  pp=object$pars_posterior
  if (!is.null(pp)) {
    pn  =intersect(c("K", "r", "psi", "m"), colnames(pp))
    prior=do.call(rbind, lapply(pn, function(p)
      data.frame(parameter=tolower(p), what="Posterior", value=pp[, p])))
    prior$parameter[prior$parameter=="m"]="n"
    jd=inp$jagsdata
    if (!is.null(jd)) for (p in c("r", "K", "psi")) {
      pr=jd[[paste0(p, ".pr")]]
      if (length(pr)==2 && all(is.finite(pr)))
        prior=rbind(prior, data.frame(parameter=tolower(p), what="Prior",
                                      value=stats::qlnorm(stats::ppoints(2000), log(pr[1]), pr[2])))}}
  cmpAssess(name, FLCore::FLQuants(ts), FLCore::FLQuants(lo), FLCore::FLQuants(up), rp,
            index=idx, resid=rsd, process=prc, curve=crv, prior=prior)}
