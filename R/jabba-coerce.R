utils::globalVariables(c("bdModel","biomass"))

setGeneric('as.biodyn', function(object, ...) standardGeneric('as.biodyn'))

tryIt<-function(x){
  rtn=try(x)
  if ("try-error"%in%is(rtn)) return(NULL)
  return(rtn)}

# jabba2biodyn, as.biodyn(list) and the list coercions are in jabba-coerce2.R
