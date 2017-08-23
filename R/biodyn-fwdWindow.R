setGeneric("fwdWindow", function(x,y,...){
  standardGeneric("fwdWindow")})
setGeneric("fbar", function(object, ...)
  standardGeneric("fbar"))

setMethod('fwdWindow', signature(x='biodyn',y="ANY"),
          function(x,y,start=dims(x)$minyear, end=dims(x)$maxyear, extend=TRUE, frequency=1,...){
            
            object=window(x,end=end)
            
            return(object)})

setMethod("fbar", signature(object="biodyn"),
          function(object, ...) harvest(object))

