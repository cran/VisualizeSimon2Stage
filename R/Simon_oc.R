

#' @title \linkS4class{Simon_oc}: Operating Characteristics of Simon's Two-Stage Design
#' 
#' @description
#' Operating characteristics of Simon's two-stage design.
#' 
#' @slot .Data \linkS4class{Simon_pr} object
#' 
#' @slot maxResp \link[base]{integer} \link[base]{vector} of length \eqn{p}, 
#' the frequencies of each regime having maximum response.  
#' The summation of `maxResp` is the number of simulation copies.
#' 
#' @slot Simon_maxResp \link[base]{integer} \link[base]{vector} of length \eqn{p}, 
#' the frequencies of each regime having maximum response and success in Simon's two-stage trial.
#' 
#' @include Simon_pr.R
#' @name Simon_oc
#' @importFrom methods setClass
#' @aliases Simon_oc-class
#' @export
setClass(Class = 'Simon_oc', contains = 'Simon_pr', slots = c(
  maxResp = 'integer',
  Simon_maxResp = 'integer'
))








#' @rdname Simon_oc
#' 
#' @param prob *named* \link[base]{numeric} \link[base]{vector}, 
#' true response rate(s) of (multiple) drug(s).
#' The `names(prob)` should be the respective keyword(s) for the drug(s).
#' 
#' @param simon \link[clinfun]{ph2simon} object
#' 
#' @param type \link[base]{character} scalar, type of Simon's two-stage design.
#' Currently supports
#' `'minimax'` (default) for minimum total sample size,
#' `'optimal'` for minimum expected total sample size *under \eqn{p_0}*,
#' `'n1'` for minimum Stage-1 sample size \eqn{n_1},
#' `'maximax'` to use up the user-provided maximum total sample size 
#' (parameter `nmax` of \link[clinfun]{ph2simon})
#' 
#' @param n1,n (optional) \link[base]{integer} scalars, Stage-1 sample size \eqn{n_1} 
#' and total sample size \eqn{n}.  Will be overridden if `simon` is given
#' 
#' @param r1,r (optional) \link[base]{integer} scalars, number of response
#' in Stage-1 \eqn{r_1} and overall \eqn{r} required *exclusively*.
#' In other words, passing Stage-1 means observing \eqn{>r_1} response.
#' Will be overridden if `simon` is given
#' 
#' @param R \link[base]{integer} scalar, number of simulations.  Default `1e4L`.
#' 
#' @param ... potential parameters, currently not in use
#' 
#' @details ..
#' 
#' @returns 
#' 
#' Function [Simon_oc()] returns \linkS4class{Simon_oc} object
#' 
#' @references 
#' \doi{10.1016/0197-2456(89)90015-9}
#' 
#' @examples 
#' library(clinfun)
#' (x = ph2simon(pu = .2, pa = .4, ep1 = .05, ep2 = .1)) 
#' Simon_oc(prob = c(A = .3, B = .2, C = .15), simon = x, type = 'minimax', R = 1e3L)
#' Simon_oc(prob = c(A = .3, B = .2, C = .15), simon = x, type = 'optimal', R = 1e3L)
#' 
#' @importFrom methods new
#' @export
Simon_oc <- function(
    prob, 
    simon, type = c('minimax', 'optimal', 'n1', 'maximax'),
    R = 1e4L,
    n1 = stop('must provide `n1`'), n = stop('must provide `n`'), 
    r1 = stop('must provide `r1`'), r = stop('must provide `r`'), 
    ...
) {
  
  if (!(pn <- length(prob)) || !is.numeric(prob) || anyNA(prob) || any(prob < 0, prob > 1)) stop('`prob` must be probabilities')
  if (!length(nm <- names(prob)) || any(!nzchar(nm))) stop('`prob` must be named probabilities')
  
  if (!missing(simon) && inherits(simon, what = 'ph2simon')) {
    # overwrite user provided `n1`, `n`, `r1`, `r`
    design <- summary.ph2simon(simon)$design[match.arg(type), , drop = TRUE]
    n1 <- design['n1']
    n <- design['n']
    r1 <- design['r1']
    r <- design['r']
  }
  
  M <- do.call(cbind, args = lapply(prob, FUN = r_simon, n1 = n1, n = n, r1 = r1, R = R)) 
  # `R x pn` 'matrix' # number of positive responses of each regimen
  
  idx <- max.col(M, ties.method = 'first') # indices of regimen being chosen (i.e., having maximum positive response) at each of the simulated trials
  idx_success <- idx[.rowSums(M > r, m = R, n = pn, na.rm = FALSE) > 0L] # faster than `Rfast::rowAny(M > r)`
  # indices of regimen {having maximum positive response} AND {succeeding in Simon trial}
  
  # eN <- .colMeans(n1 + (M > r1) * (n - n1), m = R, n = pn, na.rm = FALSE) # expected total sample size
  # no need to do this. `eN` can be obtained theoretically
  
  new(Class = 'Simon_oc', 
      Simon_pr(prob = prob, n1 = n1, r1 = r1, n = n, r = r),
      maxResp = tabulate(idx, nbins = pn),
      Simon_maxResp = tabulate(idx_success, nbins = pn))
}








#' @title Show \linkS4class{Simon_oc} Object
#' 
#' @description Show \linkS4class{Simon_oc} object
#' 
#' @param object \linkS4class{Simon_oc} object
#' 
#' @returns 
#' The \link[methods]{show} method for \linkS4class{Simon_oc} object 
#' does not have a returned value.
#' 
#' @importFrom methods setMethod show signature
#' @export
setMethod(f = show, signature = signature(object = 'Simon_oc'), definition = function(object) {
  cat(Sprintf.Simon_oc(object))
  print(autoplot.Simon_oc(object))
})








#' @importFrom ggplot2 autolayer aes geom_rect coord_polar labs theme guides guide_legend
#' @importFrom grid unit
#' @export
autolayer.Simon_oc <- function(object, ...) {
  pn <- length(prob <- object@prob)
  R <- sum(object@maxResp)
  maxResp <- object@maxResp
  Simon_maxResp <- object@Simon_maxResp
  ymax <- cumsum(maxResp)
  ymin <- c(0, ymax[-pn])
  nm <- sprintf(
    fmt = '%s; p = %.f%%\nHaving Max # of Responses: %.1f%%\nHaving Max # of Responses & Simon\'s Success: %.1f%%\nExpected Sample Size: %.1f', 
    names(prob), 1e2*prob, 1e2*maxResp/R, 1e2*Simon_maxResp/R, object@eN)
  
  list(
    geom_rect(mapping = aes(ymax = ymax, ymin = ymin, xmax = 1, xmin = 0, fill = nm), alpha = .3, stat = 'identity', colour = 'white'),
    geom_rect(mapping = aes(ymax = ymin + Simon_maxResp, ymin = ymin, xmax = 1, xmin = 0, fill = nm), stat = 'identity', colour = 'white'),
    theme(
      legend.key.spacing.y = unit(100, units = 'npc') #?? why not seeing the space??
    ),
    #guides(
    #  fill = guide_legend(byrow = TRUE)
    #),
    coord_polar(theta = 'y'),
    #coord_radial(theta = 'y'), # dont' understand what this is!!!
    labs(fill = sprintf('Regimen\n(%d Simulated Trials)', R))
  )
}




#' @importFrom ggplot2 autoplot ggplot theme_void
#' @export
autoplot.Simon_oc <- function(object, ...) {
  ggplot() + autolayer.Simon_oc(object, ...) + theme_void()
}






#' @title Short Paragraph to Describe a \linkS4class{Simon_oc} Object
#' 
#' @description
#' To create a short paragraph to describe a \linkS4class{Simon_oc} object.
#' 
#' @param model \linkS4class{Simon_oc} object
#' 
#' @param ... additional parameters, currently not in use
#' 
#' @returns 
#' [Sprintf.Simon_oc] returns a \link[base]{noquote} \link[base]{character} scalar.
#' 
#' @export Sprintf.Simon_oc
#' @export
Sprintf.Simon_oc <- function(model, ...) {
  #type <- model[['type']]
  nm <- names(prob <- model@prob)
  N <- sum(model@maxResp)
  maxResp <- model@maxResp / N
  Simon_maxResp <- model@Simon_maxResp / N
  noquote(sprintf(
    fmt = 'We simulated %d trials of each of the %d drugs %s with estimated response rates of %s, respectively, using this design. The percentage of trials with each of these drugs having the highest number of responses are %s. The percentage of trials with each of these drugs both having the highest number of responses and being accepted by the Simon\'s two-stage design are %s.', 
    N, length(prob),
    paste(sQuote(nm), collapse = ', '),
    paste(sprintf(fmt = '%.0f%%', 1e2*prob), collapse = ', '),
    paste(sprintf(fmt = '%.1f%%', 1e2*maxResp), 'for', nm, collapse = ', '),
    paste(sprintf(fmt = '%.1f%%', 1e2*Simon_maxResp), 'for', nm, collapse = ', ')))
}
