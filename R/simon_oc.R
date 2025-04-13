







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
    sort(decreasing = TRUE) |> 
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
  maxResp <- object@maxResp / R
  simon_maxResp <- object@simon_maxResp / R
  max. <- maxResp |> cumsum()
  min. <- c(0, max.[-pn])
  nm. <- sprintf(
    fmt = '%s\np = %.f%%\nE(N)=%.1f', 
    names(prob), 1e2*prob, object@eN)
  
  list(
    geom_textpath(mapping = aes(x = c(0, 0, .5), y = c(1.4, .85, 1.4), label = c(
      'max{+} & Simon\'s Success',
      'max{+}',
      R |> sprintf(fmt = '%d Simulated Trials')
    )), color = 'grey40'),
    
    geom_rect(mapping = aes(xmin = min., xmax = max., ymin = .9, ymax = 1.3, fill = nm.), alpha = .15, color = 'white'),
    geom_rect(mapping = aes(xmin = min., xmax = min. + simon_maxResp, ymin = .9, ymax = 1.3, fill = nm.), alpha = .5, color = 'white'),
    
    geom_textpath(mapping = aes(x = min. + simon_maxResp/2, y = 1.2, label = sprintf(fmt = '%.1f%%', 1e2*simon_maxResp), color = nm.), fontface = 2),
    geom_textpath(mapping = aes(x = (min. + max.)/2, y = 1.02, label = sprintf(fmt = '%.1f%%', 1e2*maxResp), color = nm.), fontface = 2),
    
    ylim(c(0, 1.5))
  
  )
  
}




#' @importFrom ggplot2 autoplot ggplot theme theme_void
#' @importFrom grid unit
#' @export
autoplot.simon_oc <- function(object, ...) {
  ggplot() + autolayer.simon_oc(object, ...) + 
    coord_polar(theta = 'x') +
    theme_void() +
    labs(fill = NULL, color = NULL) +
    theme(legend.position = 'inside') +
    theme(
      legend.key.spacing.y = unit(.02, units = 'npc')
    )
}




