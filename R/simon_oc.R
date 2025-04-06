







#' @rdname simon_oc
#' 
#' @param prob *named* \link[base]{double} \link[base]{vector}, 
#' true response rate(s) \eqn{p} of (multiple) drug(s).
#' The `names(prob)` should be the name(s) of the drug(s).
#' 
#' @param object \link[clinfun]{ph2simon} or \linkS4class{ph2simon4} object
#' 
#' @param ... parameters of function [ph2simon4()], most importantly `type`
#' 
#' @param r1,n1,r,n (optional) \link[base]{integer} scalars, see \linkS4class{ph2simon4}.
#' 
#' @param R \link[base]{integer} scalar, number of simulations.  Default `1e4L`.
#' 
#' @returns 
#' Function [simon_oc()] returns \linkS4class{simon_oc} object.
#' 
#' @keywords internal
#' @export
simon_oc <- function(prob, R, object, ...) UseMethod(generic = 'simon_oc', object = object)

#' @rdname simon_oc
#' @export simon_oc.ph2simon
#' @export
simon_oc.ph2simon <- function(prob, R = 1e4L, object, ...) object |> ph2simon4(...) |> simon_oc.ph2simon4(prob = prob, R = R)
    
    
#' @rdname simon_oc
#' @export simon_oc.ph2simon4
#' @export
simon_oc.ph2simon4 <- function(
    prob, 
    R = 1e4L,
    object, ...,
    r1 = object@r1, n1 = object@n1, r = object@r, n = object@n
) {
  
  if (!(pn <- length(prob)) || !is.numeric(prob) || anyNA(prob) || any(prob < 0, prob > 1)) stop('`prob` must be probabilities')
  if (!length(nm <- names(prob)) || any(!nzchar(nm))) stop('`prob` must be named probabilities')
  
  M <- prob |>
    lapply(FUN = r_simon.ph2simon4, r1 = r1, n1 = n1, r = r, n = n, R = R) |> # [r_simon.ph2simon4] checks if `object` is **one** design
    do.call(what = cbind) 
  # `R x pn` 'matrix' # number of positive responses of each regimen
  
  idx <- max.col(M, ties.method = 'first') # indices of regimen being chosen (i.e., having maximum positive response) at each of the simulated trials
  idx_success <- idx[.rowSums(M > r, m = R, n = pn, na.rm = FALSE) > 0L] # faster than `Rfast::rowAny(M > r)`
  # indices of regimen {having maximum positive response} AND {succeeding in Simon trial}
  
  # eN <- .colMeans(n1 + (M > r1) * (n - n1), m = R, n = pn, na.rm = FALSE) # expected total sample size
  # no need to do this. `eN` can be obtained theoretically
  
  new(Class = 'simon_oc', 
      simon_pr.ph2simon4(prob = prob, n1 = n1, r1 = r1, n = n, r = r),
      maxResp = idx |> tabulate(nbins = pn),
      simon_maxResp = idx_success |> tabulate(nbins = pn))
}














#' @importFrom ggplot2 autolayer aes geom_rect coord_polar labs
#' @export
autolayer.simon_oc <- function(object, ...) {
  pn <- length(prob <- object@prob)
  R <- sum(object@maxResp)
  maxResp <- object@maxResp
  simon_maxResp <- object@simon_maxResp
  ymax <- cumsum(maxResp)
  ymin <- c(0, ymax[-pn])
  nm <- sprintf(
    fmt = '%s; p = %.f%%\nHaving Max # of Responses: %.1f%%\nHaving Max # of Responses & Simon\'s Success: %.1f%%\nExpected Sample Size: %.1f', 
    names(prob), 1e2*prob, 1e2*maxResp/R, 1e2*simon_maxResp/R, object@eN)
  
  list(
    geom_rect(mapping = aes(ymax = ymax, ymin = ymin, xmax = 1, xmin = 0, fill = nm), alpha = .3, stat = 'identity', colour = 'white'),
    geom_rect(mapping = aes(ymax = ymin + simon_maxResp, ymin = ymin, xmax = 1, xmin = 0, fill = nm), stat = 'identity', colour = 'white'),
    coord_polar(theta = 'y'),
    #coord_radial(theta = 'y'), # dont' understand what this is!!!
    labs(fill = sprintf('Regimen\n(%d Simulated Trials)', R))
  )
}




#' @importFrom ggplot2 autoplot ggplot theme theme_void
#' @importFrom grid unit
#' @export
autoplot.simon_oc <- function(object, ...) {
  ggplot() + autolayer.simon_oc(object, ...) + 
    theme_void() +
    theme(
      legend.key.spacing.y = unit(.02, units = 'npc')
    )
}




