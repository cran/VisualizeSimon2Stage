
#' @title Summarize a Simon's Two-Stage Design
#' 
#' @description
#' Summarize a Simon's two-stage design
#' 
#' @param object a \link[clinfun]{ph2simon} or \linkS4class{ph2simon4} object
#' 
#' @param ... parameters of function [ph2simon4()], most importantly `type`
#' 
#' @returns
#' Function [summary.ph2simon()] returns a \link[base]{list} with three (3) elements
#' \describe{
#' \item{`'design'`}{\link[base]{integer} \link[base]{matrix}}
#' \item{`'EN'`}{\link[base]{double} \link[base]{matrix}}
#' \item{`'p'`}{\link[base]{double} \link[base]{matrix}}
#' }
#' 
#' @keywords internal
#' @name summary_ph2simon
#' @export summary.ph2simon4
#' @export
summary.ph2simon4 <- function(object, ...) {
  
  atr <- attributes(object)[c('r1', 'n1', 'r', 'n')]
  
  sm <- .mapply(FUN = simon_pr.ph2simon4, 
                dots = atr, 
                MoreArgs = list(prob = c(object@pu, object@pa)))
  
  ret_EN <- sm |> 
    lapply(FUN = slot, name = 'eN') |> 
    do.call(what = rbind)
  colnames(ret_EN) <- paste0('EN(', c('pu', 'pa'), ')')
  
  ret_p <- sm |> 
    #lapply(FUN = \(i) {
    lapply(FUN = \(i) { # does CRAN accept this as of Spring 2025?
      c(i@frail, i@reject[1L], 1-i@reject[2L])
    }) |> 
    do.call(what = rbind)
  colnames(ret_p) <- c('PET(pu)', 'PET(pa)', '\u03b1', '\u03b2')
  
  ret_design <- do.call(what = cbind, args = atr)
  rownames(ret_design) <- rownames(ret_EN) <- rownames(ret_p) <- object@type
  return(list(
    design = ret_design,
    EN = ret_EN,
    p = ret_p
  ))
}




#' @rdname summary_ph2simon
#' @export summary.ph2simon
#' @export
summary.ph2simon <- function(object, ...) {
  object |> 
    ph2simon4(...) |>
    summary.ph2simon4()
}


