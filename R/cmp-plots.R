# cmp-plots.R - method-agnostic assessment summaries and comparison plots
# mpb/R/cmp-plots.R
#
# Assessments from different methods (SPiCT, JABBA, biodyn, FLStock+FLBRP, ...)
# are converted to a common "cmpAssess" summary (see cmp-convert.R) and then
# plotted together by the cmp* functions, coloured by method. The prefix
# avoids clashes with the plot/kobe generics in FLCore, ggplotFL, kobe and mpb.

utils::globalVariables(c("method", "qname", "year", "data", "lower", "upper",
  "series", "residual", "xmin", "xmax", "ymin", "ymax", "fill", "outside",
  "lag", "acf", "crit", "stock", "harvest", "yield", "sp", "obs", "hat",
  "pred", "peel", "label", "x", "density", "what", "parameter", "value",
  "bmsy", "msy", "naive", "terminal"))

#' Common assessment summary for comparing methods
#'
#' A \code{cmpAssess} is a list holding the outputs of an assessment in a
#' method-agnostic form, so that the \code{cmp*} plotting functions can
#' compare assessments from different methods. Converters are
#' \code{\link{spict2cmp}}, \code{\link{flbrp2cmp}}, \code{\link{biodyn2cmp}}
#' and \code{\link{jabba2cmp}}.
#'
#' @param name label for the method, used in legends
#' @param ts \code{FLQuants} of annual values with (some of) \code{stock}
#'   (biomass at the start of the year), \code{harvest}, \code{catch},
#'   \code{bbmsy} and \code{ffmsy}; uncertainty can be held as iters
#' @param lower,upper optional \code{FLQuants} of lower and upper confidence
#'   limits, with the same names as \code{ts}
#' @param refpts \code{FLPar} or named numeric with (some of) \code{bmsy},
#'   \code{fmsy}, \code{msy}, \code{k}, \code{r} and \code{n} (Pella-Tomlinson
#'   shape, with \eqn{n=p+1})
#' @param index data.frame of observations and fits with columns
#'   \code{series}, \code{year}, \code{obs}, \code{hat} and optionally
#'   \code{lower}, \code{upper}. \code{hat} may extend beyond the last
#'   observation (one-step-ahead prediction) for hindcasting
#' @param resid data.frame of observation (e.g. one-step-ahead) residuals with
#'   columns \code{series}, \code{year}, \code{residual}; by default
#'   \eqn{\log(obs/hat)} from \code{index}
#' @param process data.frame of process residuals, columns as \code{resid};
#'   by default \eqn{\log(B_{t+1}/(B_t+P(B_t)-C_t))} with \eqn{P} interpolated
#'   from \code{curve}
#' @param curve data.frame with columns \code{stock} and \code{yield}, the
#'   production function (or equilibrium yield curve)
#' @param prior data.frame of priors and posteriors, either densities
#'   (\code{parameter}, \code{what}, \code{x}, \code{density}, with density of
#'   \eqn{\log(x)}) or draws (\code{parameter}, \code{what}, \code{value});
#'   \code{what} is \code{"Prior"} or \code{"Posterior"}
#' @param peels named list of \code{cmpAssess} for retrospective peels, names
#'   are the terminal years
#'
#' @return an object of class \code{cmpAssess}
#' @export
cmpAssess=function(name, ts, lower=NULL, upper=NULL, refpts=NULL, index=NULL,
                   resid=NULL, process=NULL, curve=NULL, prior=NULL, peels=NULL) {
  rp=c(bmsy=NA, fmsy=NA, msy=NA, k=NA, r=NA, n=NA)
  if (!is.null(refpts)) {
    v=if (is(refpts, "FLPar")) stats::setNames(c(refpts[, 1]), dimnames(refpts)[[1]]) else refpts
    v=v[names(v) %in% names(rp)]
    rp[names(v)]=as.numeric(v)}
  norm=function(l) FLCore::FLQuants(lapply(l, .cmpNorm))
  ts   =norm(ts)
  lower=if (is.null(lower)) FLCore::FLQuants() else norm(lower)
  upper=if (is.null(upper)) FLCore::FLQuants() else norm(upper)
  for (r in list(c("bbmsy", "stock", "bmsy"), c("ffmsy", "harvest", "fmsy")))
    if (is.null(ts[[r[1]]]) && !is.null(ts[[r[2]]]) && is.finite(rp[r[3]])) {
      ts[[r[1]]]=ts[[r[2]]]/rp[r[3]]
      if (!is.null(lower[[r[2]]])) lower[[r[1]]]=lower[[r[2]]]/rp[r[3]]
      if (!is.null(upper[[r[2]]])) upper[[r[1]]]=upper[[r[2]]]/rp[r[3]]}
  if (is.null(resid) && !is.null(index)) {
    i=index[is.finite(index$obs) & is.finite(index$hat) & index$obs>0 & index$hat>0, ]
    if (nrow(i)>0) resid=data.frame(series=i$series, year=i$year, residual=log(i$obs/i$hat))}
  res=structure(list(name=name, ts=ts, lower=lower, upper=upper,
                     refpts=FLCore::FLPar(rp), index=index, resid=resid,
                     process=process, curve=curve, prior=prior, peels=peels),
                class="cmpAssess")
  if (is.null(res$process)) res$process=.cmpProcess(res)
  res}

#' @method print cmpAssess
#' @export
print.cmpAssess=function(x, ...) {
  cat("cmpAssess:", x$name, "\n")
  cat("  ts:     ", paste(names(x$ts), collapse=", "), "\n")
  yrs=range(as.numeric(dimnames(x$ts[[1]])$year))
  cat("  years:  ", yrs[1], "-", yrs[2], "\n")
  rp=c(x$refpts)
  names(rp)=dimnames(x$refpts)[[1]]
  cat("  refpts: ", paste(names(rp), signif(rp, 3), sep="=", collapse=", "), "\n")
  for (s in c("index", "resid", "process", "curve", "prior"))
    if (!is.null(x[[s]])) cat("  ", s, ": ", nrow(x[[s]]), " rows\n", sep="")
  if (length(x$peels)) cat("  peels:  ", paste(names(x$peels), collapse=", "), "\n")
  invisible(x)}

# common first dimension ("quant", "all") so FLQuants from different sources combine
.cmpNorm=function(q) {
  d=dimnames(q)
  names(d)[1]="quant"
  d[[1]]="all"
  FLCore::FLQuant(array(q@.Data, dim=dim(q), dimnames=d))}

.cmpRp=function(x, nm) unname(c(x$refpts[nm, 1]))

.cmpList=function(x) {
  if (inherits(x, "cmpAssess")) x=list(x)
  nms=names(x)
  if (is.null(nms)) nms=rep("", length(x))
  nms[nms==""]=vapply(x[nms==""], function(i) i$name, character(1))
  names(x)=nms
  x}

.cmpYrs=function(q) as.numeric(dimnames(q)$year)

# median and 95% interval from an FLQuant, using iters when present
.cmpQ=function(q, fn) {
  a=q@.Data
  if (dim(a)[6]>1) c(apply(a, 2, fn)) else c(a[1, , 1, 1, 1, 1])}

.cmpTSdf=function(x, quants, ci=TRUE) {
  x=.cmpList(x)
  res=lapply(names(x), function(m) {
    a=x[[m]]
    do.call(rbind, lapply(quants, function(q) {
      if (is.null(a$ts[[q]])) return(NULL)
      v=a$ts[[q]]
      d=data.frame(method=m, qname=q, year=.cmpYrs(v),
                   data=.cmpQ(v, stats::median), lower=NA_real_, upper=NA_real_)
      if (ci) {
        if (!is.null(a$lower[[q]])) {
          d$lower=c(a$lower[[q]])[match(d$year, .cmpYrs(a$lower[[q]]))]
          d$upper=c(a$upper[[q]])[match(d$year, .cmpYrs(a$upper[[q]]))]
        } else if (dim(v)[6]>1) {
          d$lower=.cmpQ(v, function(z) stats::quantile(z, 0.025, na.rm=TRUE))
          d$upper=.cmpQ(v, function(z) stats::quantile(z, 0.975, na.rm=TRUE))}}
      d}))})
  res=do.call(rbind, res)
  res$method=factor(res$method, levels=names(x))
  res$qname =factor(res$qname,  levels=quants)
  res}

.cmpLab=c(stock="Biomass", harvest="F", catch="Catch",
          bbmsy="B/B[MSY]", ffmsy="F/F[MSY]")

.cmpTheme=function(base_size=11)
  ggplot2::theme_minimal(base_size=base_size) +
  ggplot2::theme(panel.grid.minor=ggplot2::element_blank(),
                 strip.text=ggplot2::element_text(face="bold"),
                 legend.position="bottom")

# one colour per method, shared by colour and fill so legends merge
.cmpColours=function(p, n, methods=NULL) {
  if (!is.null(methods)) {
    pal=stats::setNames(scales::hue_pal()(length(methods)), methods)
    p=p + ggplot2::scale_colour_manual(values=pal, limits=methods, name=NULL) +
      ggplot2::scale_fill_manual(values=pal, limits=methods, name=NULL)}
  if (n==1) return(p + ggplot2::guides(colour="none", fill="none"))
  p}

#' Time series of one or more assessments
#'
#' Annual estimates with 95\% confidence intervals (ribbons), coloured by
#' method. Dashed lines are the reference points (\eqn{B_{MSY}},
#' \eqn{F_{MSY}}, MSY, or 1 for the ratios) of each method.
#'
#' @param x a \code{cmpAssess} or a named list of them
#' @param quants quantities to plot, from \code{stock}, \code{harvest},
#'   \code{catch}, \code{bbmsy}, \code{ffmsy}
#' @param refs show reference points
#' @param ci show confidence intervals
#'
#' @return a ggplot
#' @export
cmpTS=function(x, quants=c("stock", "harvest", "catch"), refs=TRUE, ci=TRUE) {
  xs=.cmpList(x)
  ts=.cmpTSdf(xs, quants, ci)
  p=ggplot2::ggplot(ts, ggplot2::aes(year, data)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin=lower, ymax=upper, fill=method),
                         data=ts[is.finite(ts$lower), ], alpha=0.2) +
    ggplot2::geom_line(ggplot2::aes(colour=method), linewidth=0.7) +
    ggplot2::facet_grid(qname~., scales="free_y",
                        labeller=ggplot2::as_labeller(.cmpLab, ggplot2::label_parsed)) +
    ggplot2::expand_limits(y=0) +
    ggplot2::labs(x="Year", y=NULL, colour=NULL, fill=NULL) +
    .cmpTheme()
  if (refs) {
    rf=do.call(rbind, lapply(names(xs), function(m)
      data.frame(method=m, qname=c("stock", "harvest", "catch", "bbmsy", "ffmsy"),
                 data=c(.cmpRp(xs[[m]], "bmsy"), .cmpRp(xs[[m]], "fmsy"),
                        .cmpRp(xs[[m]], "msy"), 1, 1))))
    rf=rf[rf$qname %in% quants & is.finite(rf$data), ]
    rf$qname =factor(rf$qname, levels=quants)
    rf$method=factor(rf$method, levels=names(xs))
    p=p + ggplot2::geom_hline(ggplot2::aes(yintercept=data, colour=method),
                              data=rf, linetype=2, linewidth=0.4)}
  .cmpColours(p, length(xs), names(xs))}

#' Stock status of one or more assessments
#'
#' \eqn{B/B_{MSY}} and \eqn{F/F_{MSY}} with confidence intervals.
#'
#' @inheritParams cmpTS
#' @return a ggplot
#' @export
cmpStatus=function(x, ci=TRUE)
  cmpTS(x, c("bbmsy", "ffmsy"), refs=TRUE, ci=ci)

#' Kobe phase plot for one or more assessments
#'
#' Trajectories of \eqn{B/B_{MSY}} against \eqn{F/F_{MSY}} by method on the
#' red/yellow/green quadrants, with the terminal year and its 95\%
#' confidence intervals as crosshairs.
#'
#' @param x a \code{cmpAssess} or a named list of them
#' @param xlim,ylim axis limits; default from the data
#' @param path draw the trajectories
#'
#' @return a ggplot
#' @export
cmpKobe=function(x, xlim=NULL, ylim=NULL, path=TRUE) {
  xs=.cmpList(x)
  ts=.cmpTSdf(xs, c("bbmsy", "ffmsy"))
  bb=ts[ts$qname=="bbmsy", ]
  ff=ts[ts$qname=="ffmsy", ]
  key=paste(bb$method, bb$year)
  fi =match(key, paste(ff$method, ff$year))
  dat=data.frame(method=bb$method, year=bb$year, stock=bb$data, harvest=ff$data[fi],
                 xmin=bb$lower, xmax=bb$upper, ymin=ff$lower[fi], ymax=ff$upper[fi])
  dat=dat[is.finite(dat$stock) & is.finite(dat$harvest), ]
  trm=do.call(rbind, lapply(split(dat, dat$method, drop=TRUE), function(d) d[which.max(d$year), ]))
  if (is.null(xlim)) xlim=c(0, max(2, dat$stock, trm$xmax, na.rm=TRUE)*1.05)
  if (is.null(ylim)) ylim=c(0, max(2, dat$harvest, trm$ymax, na.rm=TRUE)*1.05)
  quad=data.frame(xmin=c(0, 1, 0, 1), xmax=c(1, Inf, 1, Inf),
                  ymin=c(1, 0, 0, 1), ymax=c(Inf, 1, 1, Inf),
                  fill=c("red", "green", "yellow", "yellow"))
  p=ggplot2::ggplot(dat, ggplot2::aes(stock, harvest)) +
    ggplot2::geom_rect(ggplot2::aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
                       data=quad, inherit.aes=FALSE, fill=c("#E8706A", "#7BC67B", "#F3E07A", "#F3E07A"))
  if (path)
    p=p + ggplot2::geom_path(ggplot2::aes(colour=method), linewidth=0.5) +
      ggplot2::geom_point(ggplot2::aes(colour=method), size=0.9, alpha=0.6)
  p=p +
    ggplot2::geom_segment(ggplot2::aes(x=stock, xend=stock, y=ymin, yend=ymax, colour=method),
                          data=trm[is.finite(trm$ymin), ], linewidth=0.8) +
    ggplot2::geom_segment(ggplot2::aes(x=xmin, xend=xmax, y=harvest, yend=harvest, colour=method),
                          data=trm[is.finite(trm$xmin), ], linewidth=0.8) +
    ggplot2::geom_point(ggplot2::aes(fill=method), data=trm, shape=21, size=3.2, colour="black") +
    ggplot2::coord_cartesian(xlim=xlim, ylim=ylim, expand=FALSE) +
    ggplot2::labs(x=expression(B/B[MSY]), y=expression(F/F[MSY]), colour=NULL, fill=NULL,
                  caption=paste("Large points: terminal year",
                                paste(unique(trm$year), collapse="/"), "with 95% CIs.")) +
    .cmpTheme()
  .cmpColours(p, length(xs), names(xs))}

#' Production functions of one or more assessments
#'
#' The production function (or equilibrium yield curve) of each method,
#' with the realised surplus production \eqn{SP_t=B_{t+1}-B_t+C_t} computed
#' from its biomass and catch, and the MSY point.
#'
#' @param x a \code{cmpAssess} or a named list of them
#' @param points show realised surplus production
#'
#' @return a ggplot
#' @export
cmpProduction=function(x, points=TRUE) {
  xs=.cmpList(x)
  crv=do.call(rbind, lapply(names(xs), function(m)
    if (!is.null(xs[[m]]$curve)) data.frame(method=m, xs[[m]]$curve[, c("stock", "yield")])))
  sp =do.call(rbind, lapply(names(xs), function(m) {
    d=.cmpSP(xs[[m]])
    if (!is.null(d)) cbind(method=m, d)}))
  msy=do.call(rbind, lapply(names(xs), function(m)
    data.frame(method=m, stock=.cmpRp(xs[[m]], "bmsy"), yield=.cmpRp(xs[[m]], "msy"))))
  msy=msy[is.finite(msy$stock) & is.finite(msy$yield), ]
  for (v in c("crv", "sp", "msy")) {
    d=get(v)
    if (!is.null(d)) {d$method=factor(d$method, levels=names(xs)); assign(v, d)}}
  p=ggplot2::ggplot(mapping=ggplot2::aes(stock, yield)) +
    ggplot2::geom_hline(yintercept=0, linetype=2, colour="grey50")
  if (points && !is.null(sp))
    p=p + ggplot2::geom_path(ggplot2::aes(colour=method), data=sp, alpha=0.3, linewidth=0.3) +
      ggplot2::geom_point(ggplot2::aes(colour=method), data=sp, alpha=0.5, size=1.2)
  if (!is.null(crv))
    p=p + ggplot2::geom_line(ggplot2::aes(colour=method), data=crv, linewidth=1)
  p=p + ggplot2::geom_point(ggplot2::aes(fill=method), data=msy, shape=21, size=3.5) +
    ggplot2::labs(x="Biomass", y="Surplus production", colour=NULL, fill=NULL,
                  caption="Lines: production functions; small points: realised surplus production; large points: MSY.") +
    .cmpTheme()
  .cmpColours(p, length(xs), names(xs))}

.cmpSP=function(a) {
  if (is.null(a$ts$stock) || is.null(a$ts$catch)) return(NULL)
  B  =.cmpQ(a$ts$stock, stats::median)
  yrs=.cmpYrs(a$ts$stock)
  C  =.cmpQ(a$ts$catch, stats::median)[match(yrs, .cmpYrs(a$ts$catch))]
  n  =length(B)
  ok =diff(yrs)==1
  d  =data.frame(year=yrs[-n], stock=B[-n], yield=B[-1]-B[-n]+C[-n],
                 next.stock=B[-1], catch=C[-n])[ok, ]
  d[is.finite(d$yield), ]}

.cmpProcess=function(a) {
  if (is.null(a$curve)) return(NULL)
  d=.cmpSP(a)
  if (is.null(d) || nrow(d)<3) return(NULL)
  crv=a$curve[order(a$curve$stock), ]
  P  =stats::approx(crv$stock, crv$yield, d$stock, rule=2)$y
  prd=d$stock+P-d$catch
  res=log(d$next.stock/prd)
  res[!is.finite(res)]=NA
  data.frame(series="Process", year=d$year, residual=res)}

#' Observed and fitted values of one or more assessments
#'
#' Observations (points) and fits (lines, with ribbons when available) for
#' each series in the \code{index} element, coloured by method.
#'
#' @param x a \code{cmpAssess} or a named list of them
#' @param series optional subset of series to plot
#'
#' @return a ggplot
#' @export
cmpFit=function(x, series=NULL) {
  xs =.cmpList(x)
  dat=do.call(rbind, lapply(names(xs), function(m) {
    i=xs[[m]]$index
    if (is.null(i)) return(NULL)
    if (is.null(i$lower)) i$lower=NA_real_
    if (is.null(i$upper)) i$upper=NA_real_
    data.frame(method=m, i[, c("series", "year", "obs", "hat", "lower", "upper")])}))
  if (is.null(dat)) stop("no index element in any assessment")
  if (!is.null(series)) dat=dat[dat$series %in% series, ]
  obs=dat[is.finite(dat$obs), ]
  obs=obs[!duplicated(obs[, c("series", "year", "obs")]), ]
  dat$method=factor(dat$method, levels=names(xs))
  p=ggplot2::ggplot(dat, ggplot2::aes(year)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin=lower, ymax=upper, fill=method),
                         data=dat[is.finite(dat$lower), ], alpha=0.2) +
    ggplot2::geom_line(ggplot2::aes(y=hat, colour=method), linewidth=0.7) +
    ggplot2::geom_point(ggplot2::aes(y=obs), data=obs, size=1.1) +
    ggplot2::facet_wrap(~series, scales="free_y", ncol=1) +
    ggplot2::expand_limits(y=0) +
    ggplot2::labs(x="Year", y=NULL, colour=NULL, fill=NULL) +
    .cmpTheme()
  .cmpColours(p, length(xs), names(xs))}

#' Runs test on residuals
#'
#' Runs test for randomness and 3-sigma limits, as in
#' \code{ss3diags::ssruns_sig3}, implemented here so that it does not
#' depend on any assessment package.
#'
#' @param x residuals
#' @param mixing alternative hypothesis, \code{"less"} tests for too few runs
#'   (positive autocorrelation)
#'
#' @return a list with \code{p} (runs test p-value) and \code{lim} (3-sigma
#'   limits)
#' @export
runsTest=function(x, mixing="less") {
  x=x[is.finite(x)]
  n=length(x)
  if (n<3) return(list(p=NA_real_, lim=c(NA, NA)))
  mr =abs(diff(x))
  sig=mean(mr)/1.128
  lim=c(-3, 3)*sig
  s  =sign(x)
  s  =s[s!=0]
  n1 =sum(s>0); n2=sum(s<0)
  if (n1==0 || n2==0) return(list(p=0, lim=lim))
  runs=1+sum(diff(s)!=0)
  mu =2*n1*n2/(n1+n2)+1
  vr =2*n1*n2*(2*n1*n2-n1-n2)/((n1+n2)^2*(n1+n2-1))
  z  =(runs-mu)/sqrt(vr)
  p  =switch(mixing,
             less     =stats::pnorm(z),
             greater  =1-stats::pnorm(z),
             two.sided=2*min(stats::pnorm(z), 1-stats::pnorm(z)))
  list(p=p, lim=lim)}

.cmpResidDF=function(xs, type) {
  res=do.call(rbind, lapply(names(xs), function(m) {
    d=xs[[m]][[type]]
    if (is.null(d) || !nrow(d)) return(NULL)
    data.frame(method=m, series=d$series, year=d$year, residual=d$residual)}))
  if (is.null(res)) stop("no '", type, "' residuals in any assessment")
  res$method=factor(res$method, levels=names(xs))
  res[is.finite(res$residual), ]}

.cmpResidPanel=function(dat, ylab, ncol=NULL) {
  dat$key=paste(dat$method, dat$series, sep=": ")
  if (length(unique(dat$method))==1) dat$key=dat$series
  pts=list(); bnd=list()
  for (s in unique(dat$key)) {
    d =dat[dat$key==s, ]
    rt=runsTest(d$residual)
    lab=paste0(s, " (runs p = ", ifelse(is.na(rt$p), "NA", signif(rt$p, 2)), ")")
    d$label  =lab
    d$outside=d$residual<rt$lim[1] | d$residual>rt$lim[2]
    pts[[s]]=d
    bnd[[s]]=data.frame(label=lab, xmin=min(d$year)-0.5, xmax=max(d$year)+0.5,
                        ymin=rt$lim[1], ymax=rt$lim[2],
                        fill=ifelse(is.na(rt$p), "Not tested", ifelse(rt$p>=0.05, "Passed", "Failed")))}
  pts=do.call(rbind, pts); bnd=do.call(rbind, bnd)
  lv =unique(pts$label)
  pts$label=factor(pts$label, levels=lv); bnd$label=factor(bnd$label, levels=lv)
  ggplot2::ggplot(pts, ggplot2::aes(year, residual)) +
    ggplot2::geom_rect(ggplot2::aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=fill),
                       data=bnd, inherit.aes=FALSE, alpha=0.35) +
    ggplot2::scale_fill_manual(values=c(Passed="#00AA00", Failed="#CC0000",
                                        `Not tested`="grey70"), name="Runs test") +
    ggplot2::geom_hline(yintercept=0, linetype=2, colour="grey40") +
    ggplot2::geom_segment(ggplot2::aes(xend=year, y=0, yend=residual),
                          colour="grey30", linewidth=0.4) +
    ggplot2::geom_point(ggplot2::aes(colour=outside), shape=21, fill="white",
                        size=1.6, stroke=0.5) +
    ggplot2::geom_smooth(method="loess", formula=y~x, se=FALSE,
                         colour="steelblue", linewidth=0.5) +
    ggplot2::scale_colour_manual(values=c(`TRUE`="#CC0000", `FALSE`="grey20"), guide="none") +
    ggplot2::facet_wrap(~label, scales="free_y", ncol=ncol) +
    ggplot2::labs(x="Year", y=ylab,
                  caption="Band: 3-sigma limits; green/red = runs test passed/failed (p >= 0.05).") +
    .cmpTheme()}

#' Observation residuals with runs tests for one or more assessments
#'
#' Residuals (e.g. one-step-ahead for SPiCT, log(obs/hat) otherwise) with
#' the 3-sigma band shaded green if the runs test is passed and red if
#' not; points outside the band are outlined in red. One panel per method
#' and series.
#'
#' @param x a \code{cmpAssess} or a named list of them
#' @param ncol number of panel columns
#'
#' @return a ggplot
#' @export
cmpRuns=function(x, ncol=NULL)
  .cmpResidPanel(.cmpResidDF(.cmpList(x), "resid"), "Residual", ncol)

#' Process residuals with runs tests for one or more assessments
#'
#' @inheritParams cmpRuns
#' @return a ggplot
#' @export
cmpProcess=function(x, ncol=NULL)
  .cmpResidPanel(.cmpResidDF(.cmpList(x), "process"), "Process residual", ncol)

#' Autocorrelation of residuals for one or more assessments
#'
#' @param x a \code{cmpAssess} or a named list of them
#' @param type \code{"resid"} (observation) or \code{"process"} residuals
#' @param lag.max maximum lag
#'
#' @return a ggplot
#' @export
cmpACF=function(x, type=c("resid", "process"), lag.max=15) {
  type=match.arg(type)
  dat=.cmpResidDF(.cmpList(x), type)
  res=do.call(rbind, lapply(split(dat, list(dat$method, dat$series), drop=TRUE), function(d) {
    a=stats::acf(d$residual, lag.max=min(lag.max, nrow(d)-1), plot=FALSE)$acf[, 1, 1]
    data.frame(method=d$method[1], series=d$series[1], lag=seq_along(a)-1, acf=a,
               crit=stats::qnorm(0.975)/sqrt(nrow(d)))}))
  res=res[res$lag>0, ]
  ggplot2::ggplot(res, ggplot2::aes(lag, acf)) +
    ggplot2::geom_hline(yintercept=0, colour="grey60") +
    ggplot2::geom_hline(ggplot2::aes(yintercept=crit), linetype=2, colour="grey40") +
    ggplot2::geom_hline(ggplot2::aes(yintercept=-crit), linetype=2, colour="grey40") +
    ggplot2::geom_col(ggplot2::aes(fill=method), width=0.6, show.legend=FALSE) +
    ggplot2::facet_grid(series~method) +
    ggplot2::labs(x="Lag", y="ACF") +
    .cmpTheme()}

#' Normal QQ plot of residuals for one or more assessments
#'
#' @inheritParams cmpACF
#' @return a ggplot
#' @export
cmpQQ=function(x, type=c("resid", "process")) {
  type=match.arg(type)
  dat=.cmpResidDF(.cmpList(x), type)
  key=paste(dat$method, dat$series)
  sw =sapply(split(dat$residual, key), function(r)
    if (length(r)>=3) stats::shapiro.test(r)$p.value else NA)
  dat$label=paste0("Shapiro p = ", signif(sw[key], 2))
  lab=unique(dat[, c("method", "series", "label")])
  ggplot2::ggplot(dat, ggplot2::aes(sample=residual)) +
    ggplot2::geom_qq_line(colour="grey40", linetype=2) +
    ggplot2::geom_qq(ggplot2::aes(fill=method), shape=21, size=1.5, show.legend=FALSE) +
    ggplot2::geom_text(ggplot2::aes(x=-Inf, y=Inf, label=label), data=lab,
                       inherit.aes=FALSE, hjust=-0.05, vjust=1.3, size=3) +
    ggplot2::facet_grid(series~method, scales="free_y") +
    ggplot2::labs(x="Theoretical quantile", y="Sample quantile") +
    .cmpTheme()}

#' Mohn's rho
#'
#' Mean relative difference between the terminal-year estimate of each
#' retrospective peel and the estimate for the same year in the reference
#' (full) run.
#'
#' @param ref \code{FLQuant} from the full run
#' @param peels list of \code{FLQuant}, one per peel
#'
#' @return numeric
#' @export
mohnsRho=function(ref, peels) {
  rr=.cmpQ(ref, stats::median)
  ry=.cmpYrs(ref)
  d=vapply(peels, function(p) {
    py=.cmpYrs(p)
    ty=max(py[is.finite(.cmpQ(p, stats::median))])
    (.cmpQ(p, stats::median)[py==ty]-rr[ry==ty])/rr[ry==ty]}, numeric(1))
  mean(d, na.rm=TRUE)}

.cmpTerminal=function(a)
  if (!is.null(a$index) && any(is.finite(a$index$obs)))
    max(a$index$year[is.finite(a$index$obs)]) else max(.cmpYrs(a$ts$stock))

#' Retrospective analysis of one or more assessments
#'
#' Peels of each method (\code{peels} element), with the confidence
#' interval of the full run and Mohn's rho annotated per method.
#'
#' @param x a \code{cmpAssess} or a named list of them; methods without peels
#'   are dropped
#' @param quants quantities to plot
#'
#' @return a ggplot
#' @export
cmpRetro=function(x, quants=c("bbmsy", "ffmsy")) {
  xs=.cmpList(x)
  xs=xs[vapply(xs, function(a) length(a$peels)>0, logical(1))]
  if (!length(xs)) stop("no assessment has retrospective peels")
  res=list(); ann=list()
  for (m in names(xs)) {
    a  =xs[[m]]
    end=.cmpTerminal(a)
    base=.cmpTSdf(stats::setNames(list(a), m), quants)
    base$peel=0
    res[[m]]=rbind(base, do.call(rbind, lapply(a$peels, function(p) {
      ty=.cmpTerminal(p)
      d =.cmpTSdf(stats::setNames(list(p), m), quants, ci=FALSE)
      d =d[d$year<=ty, ]
      d$peel=end-ty
      d})))
    rho=vapply(quants, function(q)
      mohnsRho(a$ts[[q]], lapply(a$peels, function(p) {
        ty=.cmpTerminal(p); p$ts[[q]][, as.character(.cmpYrs(p$ts[[q]])[.cmpYrs(p$ts[[q]])<=ty])]})),
      numeric(1))
    ann[[m]]=data.frame(method=m, qname=quants, label=sprintf("Mohn's rho = %.3f", rho))}
  res=do.call(rbind, res); ann=do.call(rbind, ann)
  res$peel=factor(res$peel)
  res$method=factor(res$method, levels=names(xs)); ann$method=factor(ann$method, levels=names(xs))
  res$qname=factor(res$qname, levels=quants);      ann$qname=factor(ann$qname, levels=quants)
  base=res[res$peel=="0" & is.finite(res$lower), ]
  p=ggplot2::ggplot(res, ggplot2::aes(year, data)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin=lower, ymax=upper), data=base, fill="grey70", alpha=0.4) +
    ggplot2::geom_line(ggplot2::aes(colour=peel), linewidth=0.6) +
    ggplot2::geom_text(ggplot2::aes(x=-Inf, y=Inf, label=label), data=ann,
                       inherit.aes=FALSE, hjust=-0.1, vjust=1.5, size=3.2) +
    ggplot2::facet_grid(qname~method, scales="free_y",
                        labeller=ggplot2::labeller(qname=ggplot2::as_labeller(.cmpLab, ggplot2::label_parsed))) +
    ggplot2::scale_colour_viridis_d(end=0.9, name="Peel") +
    ggplot2::expand_limits(y=0) +
    ggplot2::labs(x="Year", y=NULL) +
    .cmpTheme()
  if (all(quants %in% c("bbmsy", "ffmsy")))
    p=p + ggplot2::geom_hline(yintercept=1, linetype=2, colour="grey30")
  p}

#' Hindcast predictions from retrospective peels
#'
#' For each peel with terminal year \eqn{T}, takes the prediction of each
#' series \eqn{h} years ahead (\eqn{T+h}) from the peel's \code{index}
#' element, and the observation for that year and the naive prediction
#' (the observation in year \eqn{T}) from the full run.
#'
#' @param x a \code{cmpAssess} or a named list of them
#' @param h prediction horizon in years
#'
#' @return a data.frame with columns \code{method}, \code{series}, \code{peel}
#'   (terminal year), \code{year}, \code{obs}, \code{pred}, \code{lower},
#'   \code{upper} and \code{naive}
#' @export
cmpHindcastDF=function(x, h=1) {
  xs=.cmpList(x)
  res=lapply(names(xs), function(m) {
    a=xs[[m]]
    if (!length(a$peels) || is.null(a$index)) return(NULL)
    do.call(rbind, lapply(a$peels, function(p) {
      if (is.null(p$index)) return(NULL)
      ty=.cmpTerminal(p)
      do.call(rbind, lapply(unique(p$index$series), function(s) {
        pi=p$index[p$index$series==s & p$index$year==ty+h, ]
        bi=a$index[a$index$series==s, ]
        if (!nrow(pi) || !is.finite(pi$hat[1])) return(NULL)
        data.frame(method=m, series=s, peel=ty, year=ty+h,
                   obs  =bi$obs[match(ty+h, bi$year)], pred=pi$hat[1],
                   lower=if (is.null(pi$lower)) NA else pi$lower[1],
                   upper=if (is.null(pi$upper)) NA else pi$upper[1],
                   naive=bi$obs[match(ty, bi$year)])}))}))})
  do.call(rbind, res)}

#' Mean absolute scaled error
#'
#' MASE of hindcast predictions relative to a naive (random walk)
#' prediction, on the log scale by default. MASE < 1 means the model
#' predicts better than the naive prediction.
#'
#' @param obs observations
#' @param pred model predictions
#' @param naive naive predictions (e.g. the last observation)
#' @param log calculate errors on the log scale
#'
#' @return numeric
#' @export
cmpMASE=function(obs, pred, naive, log=TRUE) {
  if (log) {obs=log(obs); pred=log(pred); naive=log(naive)}
  ok=is.finite(obs) & is.finite(pred) & is.finite(naive)
  mean(abs(obs[ok]-pred[ok]))/mean(abs(obs[ok]-naive[ok]))}

#' Hindcast cross-validation of one or more assessments
#'
#' Observations, the fits of each retrospective peel, and the
#' \eqn{h}-year-ahead hindcast predictions with their confidence intervals,
#' with the MASE (see \code{\link{cmpMASE}}) annotated per method and series.
#'
#' @param x a \code{cmpAssess} or a named list of them
#' @param h prediction horizon in years
#' @param nyears number of years shown before the first hindcast year
#' @param series optional subset of series
#'
#' @return a ggplot
#' @export
cmpHindcast=function(x, h=1, nyears=10, series=NULL) {
  xs=.cmpList(x)
  hc=cmpHindcastDF(xs, h)
  if (is.null(hc)) stop("no hindcasts: assessments need peels with index predictions beyond their last year")
  if (!is.null(series)) hc=hc[hc$series %in% series, ]
  y0=min(hc$peel)-nyears
  lns=do.call(rbind, lapply(unique(hc$method), function(m) {
    a=xs[[m]]
    end=.cmpTerminal(a)
    rbind(data.frame(method=m, a$index[, c("series", "year", "hat")], peel=0)[a$index$year<=end, ],
          do.call(rbind, lapply(a$peels, function(p) {
            ty=.cmpTerminal(p)
            i=p$index[p$index$year<=ty+h, ]
            data.frame(method=m, i[, c("series", "year", "hat")], peel=end-ty)})))}))
  lns=lns[lns$year>=y0 & lns$series %in% unique(hc$series), ]
  obs=do.call(rbind, lapply(unique(hc$method), function(m) {
    i=xs[[m]]$index
    data.frame(method=m, i[is.finite(i$obs), c("series", "year", "obs")])}))
  obs=obs[obs$year>=y0 & obs$series %in% unique(hc$series), ]
  hc$peel=max(sapply(xs[unique(hc$method)], .cmpTerminal))-hc$peel
  ms=do.call(rbind, lapply(split(hc, list(hc$method, hc$series), drop=TRUE), function(d)
    data.frame(method=d$method[1], series=d$series[1],
               label=sprintf("MASE = %.2f", cmpMASE(d$obs, d$pred, d$naive)))))
  lv=sort(unique(c(lns$peel, hc$peel)))
  lns$peel=factor(lns$peel, levels=lv); hc$peel=factor(hc$peel, levels=lv)
  ml=intersect(names(xs), unique(hc$method))
  lns$method=factor(lns$method, levels=ml); hc$method=factor(hc$method, levels=ml)
  obs$method=factor(obs$method, levels=ml); ms$method=factor(ms$method, levels=ml)
  ggplot2::ggplot(lns, ggplot2::aes(year, hat)) +
    ggplot2::geom_line(ggplot2::aes(colour=peel), linewidth=0.6) +
    ggplot2::geom_point(ggplot2::aes(y=obs), data=obs, size=1.5) +
    ggplot2::geom_errorbar(ggplot2::aes(y=pred, ymin=lower, ymax=upper, colour=peel),
                           data=hc[is.finite(hc$lower), ], width=0.2) +
    ggplot2::geom_point(ggplot2::aes(y=pred, colour=peel), data=hc, shape=21,
                        fill="white", size=2.4, stroke=1) +
    ggplot2::geom_text(ggplot2::aes(x=-Inf, y=Inf, label=label), data=ms,
                       inherit.aes=FALSE, hjust=-0.1, vjust=1.5, size=3.2) +
    ggplot2::facet_grid(series~method, scales="free_y") +
    ggplot2::scale_colour_viridis_d(end=0.9, name="Peel") +
    ggplot2::labs(x="Year", y=NULL,
                  caption=paste0("Lines: fits to each peel; open circles: ", h,
                                 "-year-ahead hindcasts.")) +
    .cmpTheme()}

#' Prior and posterior densities for one or more assessments
#'
#' Plots the \code{prior} element of each assessment. Densities are of
#' \eqn{\log(\theta)}, plotted against \eqn{\theta} on a log axis; draws are
#' converted to densities of \eqn{\log(\theta)}. Priors are dashed and
#' posteriors solid, coloured by method.
#'
#' @param x a \code{cmpAssess} or a named list of them
#' @param parameters optional subset of parameters
#'
#' @return a ggplot
#' @export
cmpPriors=function(x, parameters=NULL) {
  xs=.cmpList(x)
  res=do.call(rbind, lapply(names(xs), function(m) {
    d=xs[[m]]$prior
    if (is.null(d)) return(NULL)
    if (!"density" %in% names(d)) d=.cmpDraws2Density(d)
    data.frame(method=m, d[, c("parameter", "what", "x", "density")])}))
  if (is.null(res)) stop("no prior element in any assessment")
  if (!is.null(parameters)) res=res[res$parameter %in% parameters, ]
  res$method=factor(res$method, levels=names(xs))
  res$what  =factor(res$what, levels=c("Prior", "Posterior"))
  res$grp   =paste(res$method, res$what)
  p=ggplot2::ggplot(res, ggplot2::aes(x, density, group=grp)) +
    ggplot2::geom_area(ggplot2::aes(fill=method), data=res[res$what=="Posterior", ],
                       position="identity", alpha=0.2) +
    ggplot2::geom_line(ggplot2::aes(colour=method, linetype=what), linewidth=0.6) +
    ggplot2::scale_linetype_manual(values=c(Prior=2, Posterior=1), name=NULL) +
    ggplot2::scale_x_log10(breaks=scales::breaks_log(4),
                           labels=function(z) format(z, digits=2, drop0trailing=TRUE, trim=TRUE)) +
    ggplot2::facet_wrap(~parameter, scales="free") +
    ggplot2::labs(x="Value (log scale)", y="Density of log(parameter)", colour=NULL, fill=NULL) +
    .cmpTheme()
  .cmpColours(p, length(xs), names(xs))}

.cmpDraws2Density=function(d) {
  do.call(rbind, lapply(split(d, list(d$parameter, d$what), drop=TRUE), function(z) {
    v=log(z$value[is.finite(z$value) & z$value>0])
    if (length(v)<5) return(NULL)
    k=stats::density(v)
    data.frame(parameter=z$parameter[1], what=z$what[1], x=exp(k$x), density=k$y)}))}
