



#' @title Show S4 Objects
#' 
#' @description Show S4 object.
#' 
#' @param object S4 objects, e.g., \linkS4class{simon_oc}, \linkS4class{simon_pr}
#' 
#' @returns 
#' The \link[methods]{show} method for S4 objects 
#' does not have a returned value.
#' 
#' @keywords internal
#' @name show_simon
#' @aliases show,simon_oc-method
#' @export
setMethod(f = show, signature = 'simon_oc', definition = function(object) {
  object |> autoplot.simon_oc() |> print()
})



#' @rdname show_simon
#' @aliases show,simon_pr-method
#' @export
setMethod(f = show, signature = 'simon_pr', definition = function(object) {
  object |> 
    as_flextable.simon_pr() |> 
    print()
  # object |> autoplot.simon_pr() |> print() # skip this
})



#' @rdname show_simon
#' @aliases show,ph2simon4-method
#' @export
setMethod(f = show, signature = 'ph2simon4', definition = function(object) {
  
  out <- object |> summary.ph2simon4()
  out$EN[] <- sprintf(fmt = '%.1f', out$EN)
  out$p[] <- sprintf(fmt = '%.1f%%', 1e2*out$p) 
  
  cat('\n Simon\'s 2-Stage Phase II Design\n\n')
  sprintf(fmt = 'Unacceptable Response Rate: %.1f%%\n', 1e2*object@pu) |> cat()
  sprintf(fmt = 'Desirable Response Rate: %.1f%%\n', 1e2*object@pa) |> cat()
  sprintf(fmt = 'Controlled Error Rates: \u03b1 \u2264 %.f%%, \u03b2 \u2264 %.f%%\n', 1e2*object@alpha, 1e2*object@beta) |> cat()
  sprintf(fmt = 'Maximum Sample Size Allowed: %d\n\n', object@nmax) |> cat()
  
  cbind(
    out$design, 
    out$EN, 
    out$p
  ) |> print.default(right = TRUE, quote = FALSE)
  
  cat('\n')
  
})
