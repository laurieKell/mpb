#' fitPella
#' 
#' @title fitPella 
#' 
#' @description 
#' 
#' @name 
#' @param 
#' 
#' @aliases
#' 
#' @docType method
#' 
#' @rdname 
#' @seealso 
#' 
#' @examples
#' \dontrun{
#' }


  fitPella
Code: function(object, dir = tempdir(), cmdOps = paste("-maxfn 500
                                                       -iprint 0"), lav = FALSE, maxF = 2.5, silent = !TRUE)
  Docs: function(object, index = index, exeNm = "pella", package =
                   "mpb", dir = tempdir(), cmdOps = paste("-maxfn 500
                                                          -iprint 0"), lav = FALSE, maxF = 2.5)
    Argument names in code not in docs:
  silent
Argument names in docs not in code:
  index exeNm package
Mismatches in argument names (first 3):
  Position: 2 Code: dir Docs: index
Position: 3 Code: cmdOps Docs: exeNm
Position: 4 Code: lav Docs: package

Codoc mismatches from documentation object 'oem,FLStock-method':
  \S4method{oem}{FLStock}
Code: function(object, sel = FLQuant(FLQuant(1, dimnames =
                                               dimnames(harvest(object)))), timing = 0.5,
               fish.dependent = TRUE, effort = c("f", "h"), mass =
                 TRUE)
  Docs: function(object, cv = rlnorm(dim(stock(object))[6], FLQuant(0,
                                                                    dimnames = dimnames(stock(object))[-6]), 0.3), timing
                 = 0.5, mult = TRUE, fishDepend = FALSE, effort =
                   c("f", "h"), mass = TRUE, q = FLQuant(cumprod(1 +
                                                                   rep(0, dim(fbar(object))[2])), dimnames =
                                                           dimnames(fbar(object))), sel = FLQuant(FLQuant(1,
                                                                                                          dimnames = dimnames(harvest(object)))), seed = NULL)
    Argument names in code not in docs:
  fish.dependent
Argument names in docs not in code:
  cv mult fishDepend q seed
Mismatches in argument names (first 3):
  Position: 2 Code: sel Docs: cv
Position: 4 Code: fish.dependent Docs: mult
Position: 5 Code: effort Docs: fishDepend

#' plot 
#' 
#' @title  
#' 
#' @description 
#' 
#' @name 
#' @param 
#' 
#' @aliases
#' 
#' @docType method
#' 
#' @rdname 
#' @seealso 
#' 
#' @examples
#' \dontrun{
#' }

Codoc mismatches from documentation object 'plot,biodyn,missing-method':
  \S4method{plot}{biodyn,missing}
Code: function(x, y, probs = c(0.9, 0.75, 0.5, 0.25, 0.1), na.rm =
                 FALSE, type = 7, worm = NULL, fn = list(Stock =
                                                           function(x) stock(x), Harvest = function(x)
                                                             harvest(x), Yield = function(x) catch(x)), facet =
                 facet_wrap(~qname, scales = "free", ncol = 1), ...)
  Docs: function(x, y, probs = c(0.95, 0.75, 0.5, 0.25, 0.05), na.rm =
                   FALSE, type = 7, worm = NULL, fn = list(Stock =
                                                             function(x) stock(x), Harvest = function(x)
                                                               harvest(x), Yield = function(x) catch(x)), facet =
                   facet_wrap(~qname, scales = "free", ncol = 1), ...)
    Mismatches in argument default values:
  Name: 'probs' Code: c(0.9, 0.75, 0.5, 0.25, 0.1) Docs: c(0.95, 0.75, 0.5, 0.25, 0.05)

#' plotIndex
#' 
#' @title  
#' 
#' @description 
#' 
#' @name 
#' @param 
#' 
#' @aliases
#' 
#' @docType method
#' 
#' @rdname 
#' @seealso 
#' 
#' @examples
#' \dontrun{
#' }


#' 
#' 
#' @title  
#' 
#' @description 
#' 
#' @name 
#' @param 
#' 
#' @aliases
#' 
#' @docType method
#' 
#' @rdname 
#' @seealso 
#' 
#' @examples
#' \dontrun{
#' }

#' 
#' 
#' @title  
#' 
#' @description 
#' 
#' @name 
#' @param 
#' 
#' @aliases
#' 
#' @docType method
#' 
#' @rdname 
#' @seealso 
#' 
#' @examples
#' \dontrun{
#' } 
 

#' 
#' 
#' @title  
#' 
#' @description 
#' 
#' @name 
#' @param 
#' 
#' @aliases
#' 
#' @docType method
#' 
#' @rdname 
#' @seealso 
#' 
#' @examples
#' \dontrun{
#' }

Undocumented arguments in documentation object 'diagsFn'
???res???

#' 
#' 
#' @title  
#' 
#' @description 
#' 
#' @name 
#' @param 
#' 
#' @aliases
#' 
#' @docType method
#' 
#' @rdname 
#' @seealso 
#' 
#' @examples
#' \dontrun{
#' }

Documented arguments not in \usage in documentation object 'diagsFn':
  ???object??? ???index???

#' 
#' 
#' @title  
#' 
#' @description 
#' 
#' @name 
#' @param 
#' 
#' @aliases
#' 
#' @docType method
#' 
#' @rdname 
#' @seealso 
#' 
#' @examples
#' \dontrun{
#' }

Undocumented arguments in documentation object 'control<-,biodyn,FLPar-method'
???object??? ???value???

#' 
#' 
#' @title  
#' 
#' @description 
#' 
#' @name 
#' @param 
#' 
#' @aliases
#' 
#' @docType method
#' 
#' @rdname 
#' @seealso 
#' 
#' @examples
#' \dontrun{
#' }

Undocumented arguments in documentation object 'biodyn'
???params???

#' 
#' 
#' @title  
#' 
#' @description 
#' 
#' @name 
#' @param 
#' 
#' @aliases
#' 
#' @docType method
#' 
#' @rdname 
#' @seealso 
#' 
#' @examples
#' \dontrun{
#' }


Undocumented arguments in documentation object 'biodyns,biodyn-method'
???object???

#' 
#' 
#' @title  
#' 
#' @description 
#' 
#' @name 
#' @param 
#' 
#' @aliases
#' 
#' @docType method
#' 
#' @rdname 
#' @seealso 
#' 
#' @examples
#' \dontrun{
#' }


Undocumented arguments in documentation object 'boot,biodyn-method'
???object??? ???run???

#' 
#' 
#' @title  
#' 
#' @description 
#' 
#' @name 
#' @param 
#' 
#' @aliases
#' 
#' @docType method
#' 
#' @rdname 
#' @seealso 
#' 
#' @examples
#' \dontrun{
#' }
#'
 
Documented arguments not in \usage in documentation object 'boot,biodyn-method':
  ???object;???

#' 
#' 
#' @title  
#' 
#' @description 
#' 
#' @name 
#' @param 
#' 
#' @aliases
#' 
#' @docType method
#' 
#' @rdname 
#' @seealso 
#' 
#' @examples
#' \dontrun{
#' }

Duplicated \argument entries in documentation object 'diags':
  ???object??? ???value??? ???...???

#' 
#' 
#' @title  
#' 
#' @description 
#' 
#' @name 
#' @param 
#' 
#' @aliases
#' 
#' @docType method
#' 
#' @rdname 
#' @seealso 
#' 
#' @examples
#' \dontrun{
#' }

Documented arguments not in \usage in documentation object 'diags':
  ???value???

#' 
#' 
#' @title  
#' 
#' @description 
#' 
#' @name 
#' @param 
#' 
#' @aliases
#' 
#' @docType method
#' 
#' @rdname 
#' @seealso 
#' 
#' @examples
#' \dontrun{
#' }

Undocumented arguments in documentation object 'controlFn'
???object???

#' 
#' 
#' @title  
#' 
#' @description 
#' 
#' @name 
#' @param 
#' 
#' @aliases
#' 
#' @docType method
#' 
#' @rdname 
#' @seealso 
#' 
#' @examples
#' \dontrun{
#' }

Documented arguments not in \usage in documentation object 'controlFn':
  ???r??? ???k??? ???p??? ???b0??? ???phaseR??? ???phaseK??? ???phaseP??? ???phaseB0??? ???min??? ???max???

#' 
#' 
#' @title  
#' 
#' @description 
#' 
#' @name 
#' @param 
#' 
#' @aliases
#' 
#' @docType method
#' 
#' @rdname 
#' @seealso 
#' 
#' @examples
#' \dontrun{
#' }

Undocumented arguments in documentation object 'fitPella'
???object??? ???index??? ???exeNm??? ???package??? ???dir??? ???cmdOps??? ???lav??? ???maxF???

#' 
#' 
#' @title  
#' 
#' @description 
#' 
#' @name 
#' @param 
#' 
#' @aliases
#' 
#' @docType method
#' 
#' @rdname 
#' @seealso 
#' 
#' @examples
#' \dontrun{
#' }

Undocumented arguments in documentation object 'fnProfile'
???x???

#' 
#' 
#' @title  
#' 
#' @description 
#' 
#' @name 
#' @param 
#' 
#' @aliases
#' 
#' @docType method
#' 
#' @rdname 
#' @seealso 
#' 
#' @examples
#' \dontrun{
#' }

Documented arguments not in \usage in documentation object 'fnProfile':
  ???x:???

#' 
#' 
#' @title  
#' 
#' @description 
#' 
#' @name 
#' @param 
#' 
#' @aliases
#' 
#' @docType method
#' 
#' @rdname 
#' @seealso 
#' 
#' @examples
#' \dontrun{
#' }

Undocumented arguments in documentation object 'fwd,biodyn,missing-method'
???catch??? ???harvest??? ???stock??? ???hcr??? ???pe??? ???peMult??? ???minF??? ???maxF??? ???bounds???
???lag??? ???end??? ???starvationRations???

#' 
#' 
#' @title  
#' 
#' @description 
#' 
#' @name 
#' @param 
#' 
#' @aliases
#' 
#' @docType method
#' 
#' @rdname 
#' @seealso 
#' 
#' @examples
#' \dontrun{
#' }

Documented arguments not in \usage in documentation object 'fwd,biodyn,missing-method':
  ???ctrl???

#' 
#' 
#' @title  
#' 
#' @description 
#' 
#' @name 
#' @param 
#' 
#' @aliases
#' 
#' @docType method
#' 
#' @rdname 
#' @seealso 
#' 
#' @examples
#' \dontrun{
#' }

Undocumented arguments in documentation object 'hcr,biodyn,missing-method'
???params??? ???yr??? ???byr??? ???hyr??? ???tac??? ???tacMn??? ???bndF??? ???bndTac??? ???iaF??? ???iaTac???
???maxF???

#' 
#' 
#' @title  
#' 
#' @description 
#' 
#' @name 
#' @param 
#' 
#' @aliases
#' 
#' @docType method
#' 
#' @rdname 
#' @seealso 
#' 
#' @examples
#' \dontrun{
#' }

Undocumented arguments in documentation object 'kobe'
???what??? ???probs??? ???year??? ???nwrms??? ???sim??? ???drop???

#' 
#' 
#' @title  
#' 
#' @description 
#' 
#' @name 
#' @param 
#' 
#' @aliases
#' 
#' @docType method
#' 
#' @rdname 
#' @seealso 
#' 
#' @examples
#' \dontrun{
#' }

Documented arguments not in \usage in documentation object 'kobe':
  ???method???

#' 
#' 
#' @title  
#' 
#' @description 
#' 
#' @name 
#' @param 
#' 
#' @aliases
#' 
#' @docType method
#' 
#' @rdname 
#' @seealso 
#' 
#' @examples
#' \dontrun{
#' }

Undocumented arguments in documentation object 'mixnorm'
???n??? ???bin??? ???left???

#' 
#' 
#' @title  
#' 
#' @description 
#' 
#' @name 
#' @param 
#' 
#' @aliases
#' 
#' @docType method
#' 
#' @rdname 
#' @seealso 
#' 
#' @examples
#' \dontrun{
#' }

Undocumented arguments in documentation object 'resample,FLQuant-method'
???x???

#' 
#' 
#' @title  
#' 
#' @description 
#' 
#' @name 
#' @param 
#' 
#' @aliases
#' 
#' @docType method
#' 
#' @rdname 
#' @seealso 
#' 
#' @examples
#' \dontrun{
#' }

Duplicated \argument entries in documentation object 'resample,FLQuant-method':
  ???object??? ???dim??? ???size??? ???replace??? ???...???

#' 
#' 
#' @title  
#' 
#' @description 
#' 
#' @name 
#' @param 
#' 
#' @aliases
#' 
#' @docType method
#' 
#' @rdname 
#' @seealso 
#' 
#' @examples
#' \dontrun{
#' }

Documented arguments not in \usage in documentation object 'resample,FLQuant-method':
  ???object???

#' 
#' 
#' @title  
#' 
#' @description 
#' 
#' @name 
#' @param 
#' 
#' @aliases
#' 
#' @docType method
#' 
#' @rdname 
#' @seealso 
#' 
#' @examples
#' \dontrun{
#' }

Undocumented arguments in documentation object 'fmsy'
???params???

#' 
#' 
#' @title  
#' 
#' @description 
#' 
#' @name 
#' @param 
#' 
#' @aliases
#' 
#' @docType method
#' 
#' @rdname 
#' @seealso 
#' 
#' @examples
#' \dontrun{
#' }

Documented arguments not in \usage in documentation object 'fmsy':
  ???\code{object}???

#' 
#' 
#' @title  
#' 
#' @description 
#' 
#' @name 
#' @param 
#' 
#' @aliases
#' 
#' @docType method
#' 
#' @rdname 
#' @seealso 
#' 
#' @examples
#' \dontrun{
#' }

Undocumented arguments in documentation object 'nll,FLQuant,FLQuant,FLPar-method'
???index???

#' 
#' 
#' @title  
#' 
#' @description 
#' 
#' @name 
#' @param 
#' 
#' @aliases
#' 
#' @docType method
#' 
#' @rdname 
#' @seealso 
#' 
#' @examples
#' \dontrun{
#' }

Undocumented arguments in documentation object 'oem,FLStock-method'
???timing??? ???seed???

#' 
#' 
#' @title  
#' 
#' @description 
#' 
#' @name 
#' @param 
#' 
#' @aliases
#' 
#' @docType method
#' 
#' @rdname 
#' @seealso 
#' 
#' @examples
#' \dontrun{
#' }

Undocumented arguments in documentation object 'plot,biodyn,missing-method'
???...??? ???data???
Duplicated \argument entries in documentation object 'plot,biodyn,missing-method':
  ???x???

#' 
#' 
#' @title  
#' 
#' @description 
#' 
#' @name 
#' @param 
#' 
#' @aliases
#' 
#' @docType method
#' 
#' @rdname 
#' @seealso 
#' 
#' @examples
#' \dontrun{
#' }

Undocumented arguments in documentation object 'plotEql,biodyn,missing-method'
???data???

#' 
#' 
#' @title  
#' 
#' @description 
#' 
#' @name 
#' @param 
#' 
#' @aliases
#' 
#' @docType method
#' 
#' @rdname 
#' @seealso 
#' 
#' @examples
#' \dontrun{
#' }

Documented arguments not in \usage in documentation object 'plotEql,biodyn,missing-method':
  ???x???

#' 
#' 
#' @title  
#' 
#' @description 
#' 
#' @name 
#' @param 
#' 
#' @aliases
#' 
#' @docType method
#' 
#' @rdname 
#' @seealso 
#' 
#' @examples
#' \dontrun{
#' }

Undocumented arguments in documentation object 'plotIndex,FLQuants-method'
???data??? ???facet??? ???...???

#' 
#' 
#' @title  
#' 
#' @description 
#' 
#' @name 
#' @param 
#' 
#' @aliases
#' 
#' @docType method
#' 
#' @rdname 
#' @seealso 
#' 
#' @examples
#' \dontrun{
#' }

Documented arguments not in \usage in documentation object 'plotIndex,FLQuants-method':
  ???x???

#' 
#' 
#' @title  
#' 
#' @description 
#' 
#' @name 
#' @param 
#' 
#' @aliases
#' 
#' @docType method
#' 
#' @rdname 
#' @seealso 
#' 
#' @examples
#' \dontrun{
#' }


Undocumented arguments in documentation object 'plotMSEfn'
???mp??? ???om??? ???brp???
Documented arguments not in \usage in documentation object 'plotMSEfn':
  ???x??? ???y??? ???z???

#' 
#' 
#' @title  
#' 
#' @description 
#' 
#' @name 
#' @param 
#' 
#' @aliases
#' 
#' @docType method
#' 
#' @rdname 
#' @seealso 
#' 
#' @examples
#' \dontrun{
#' }

Undocumented arguments in documentation object 'plotProduction,biodyn,missing-method'
???data???
Documented arguments not in \usage in documentation object 'plotProduction,biodyn,missing-method':
  ???object???

#' 
#' 
#' @title  
#' 
#' @description 
#' 
#' @name 
#' @param 
#' 
#' @aliases
#' 
#' @docType method
#' 
#' @rdname 
#' @seealso 
#' 
#' @examples
#' \dontrun{
#' }

Undocumented arguments in documentation object 'profileFn'
???r??? ???k??? ???p??? ???b0??? ???min??? ???max???
Documented arguments not in \usage in documentation object 'profileFn':
  ???fitted:??? ???which:??? ???range;??? ???fn:??? ???run:???

#' 
#' 
#' @title  
#' 
#' @description 
#' 
#' @name 
#' @param 
#' 
#' @aliases
#' 
#' @docType method
#' 
#' @rdname 
#' @seealso 
#' 
#' @examples
#' \dontrun{
#' }


Undocumented arguments in documentation object 'mseBiodyn'
???eql??? ???start??? ???end??? ???interval??? ???oem??? ???hcrPar??? ???bndF??? ???bndTac??? ???maxF???
???omega??? ???refB??? ???qTrend???
Documented arguments not in \usage in documentation object 'mseBiodyn':
  ???what??? ???mult???

#' 
#' 
#' @title  
#' 
#' @description 
#' 
#' @name 
#' @param 
#' 
#' @aliases
#' 
#' @docType method
#' 
#' @rdname 
#' @seealso 
#' 
#' @examples
#' \dontrun{
#' }

Undocumented arguments in documentation object 'setControl<-,biodyn,FLPar-method'
???min??? ???max???

#' 
#' 
#' @title  
#' 
#' @description 
#' 
#' @name 
#' @param 
#' 
#' @aliases
#' 
#' @docType method
#' 
#' @rdname 
#' @seealso 
#' 
#' @examples
#' \dontrun{
#' }

Undocumented arguments in documentation object 'sim,missing,missing-method'
???p??? ???b0???
Documented arguments not in \usage in documentation object 'sim,missing,missing-method':
  ???model???

#' 
#' 
#' @title  
#' 
#' @description 
#' 
#' @name 
#' @param 
#' 
#' @aliases
#' 
#' @docType method
#' 
#' @rdname 
#' @seealso 
#' 
#' @examples
#' \dontrun{
#' }


Undocumented arguments in documentation object 'timeSeries,FLStock,FLPar-method'
???params??? ???df???

#' timeSeries 
#' 
#' @title timeSeries 
#' 
#' @description 
#' 
#' @name 
#' @param 
#' 
#' @aliases timeSeries timeSeries-method 
#'          timeSeries, aspic,missing-method
#'          timeSeries, biodyn,missing-method
#' 
#' @docType method
#' 
#' @rdname 
#' @seealso 
#' 
#' @examples
#' \dontrun{
#' }


#' control<-
#' 
#' @title control<- 
#' 
#' @description 
#' 
#' @name 
#' @param 
#' 
#' @aliases
#' 
#' @docType method
#' 
#' @rdname 
#' @seealso 
#' 
#' @examples
#' \dontrun{
#' }

#' @title eql
#' 
#' @description 
#' 

#' @title feasible 
#' 
#' @description 
#' 

#' @title fwd
#' 
#' @description 
#' 

#' @title hcr
#' 
#' @description 
#' 

#' @title nll 
#' 
#' @description 
#' 

#' @title oem 
#' 
#' @description 
#' 

#' @title plotEql
#' 
#' @description 
#' 

#' @title plotHcr 
#' 
#' @description 
#' 

#' @title plotProduction 
#' 
#' @description 
#' 

#' @title refptSD
#' 
#' @description 
#' 

#' @title resample 
#' 
#' @description 
#' 

#' @title setControl<- 
#' 
#' @description 
#' 

#' @title sim
#' 
#' @description 
#' 

#' @title xval
#' 
#' @description 
#' 
