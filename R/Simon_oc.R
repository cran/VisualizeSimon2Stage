

#' @title \linkS4class{Simon_oc}: Operating Characteristics of Simon's Two-Stage Design
#' 
#' @description
#' Operating characteristics of Simon's two-stage design.
#' 
# @slot .Data \linkS4class{Simon_pr} object
#' 
#' @slot maxResp \link[base]{integer} \link[base]{vector} of same length as \eqn{p}, 
#' the frequencies of each regime having maximum response.  
#' The summation of `maxResp` is the number of simulation copies.
#' 
#' @slot Simon_maxResp \link[base]{integer} \link[base]{vector} of same length as \eqn{p}, 
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
#' @param prob *named* \link[base]{double} \link[base]{vector}, 
#' true response rate(s) \eqn{p} of (multiple) drug(s).
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
#' (parameter `nmax` of function \link[clinfun]{ph2simon})
#' 
#' @param n1,n (optional) \link[base]{integer} scalars, Stage-1 sample size \eqn{n_1} 
#' and total sample size \eqn{n}.  Overridden if `simon` is given
#' 
#' @param r1,r (optional) \link[base]{integer} scalars, number of response
#' in Stage-1 \eqn{r_1} and overall \eqn{r} required *exclusively*,
#' i.e., passing Stage-1 means observing \eqn{>r_1} response.
#' Overridden if `simon` is given
#' 
#' @param R \link[base]{integer} scalar, number of simulations.  Default `1e4L`.
#' 
#' @param ... potential parameters, currently not in use
#' 
#' @details ..
#' 
#' @returns 
#' 
#' Function [Simon_oc] returns \linkS4class{Simon_oc} object.
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
#' @keywords internal
#' @export
setMethod(f = show, signature = signature(object = 'Simon_oc'), definition = function(object) {
  cat(Sprintf.Simon_oc(object))
  print(autoplot.Simon_oc(object))
})








#' @importFrom ggplot2 autolayer aes geom_rect coord_polar scale_fill_discrete ylim labs
#' @importFrom geomtextpath geom_textpath
#' @export
autolayer.Simon_oc <- function(object, ...) {
  pn <- length(prob <- object@prob)
  R <- sum(object@maxResp)
  maxResp <- object@maxResp
  Simon_maxResp <- object@Simon_maxResp
  max_ <- cumsum(maxResp)
  min_ <- c(0, max_[-pn])
  
  nm <- names(prob)
  legd <- sprintf(fmt = '%s; E(N)=%.1f', names(prob), object@eN)
  nm1 <- sprintf(fmt = 'Pr(%s)=%.f%%', names(prob), 1e2*prob)
  pr_simon <- sprintf(fmt = '%.1f%%', 1e2*Simon_maxResp/R)
  pr <- sprintf(fmt = '%.1f%%', 1e2*maxResp/R)
  
  list(
    geom_rect(mapping = aes(xmax = max_, xmin = min_, ymax = 1, ymin = .7, fill = nm), alpha = .1, stat = 'identity', colour = 'white', show.legend = FALSE),
    geom_rect(mapping = aes(xmax = min_ + Simon_maxResp, xmin = min_, ymax = 1, ymin = .7, fill = nm), alpha = .3, stat = 'identity', colour = 'white'),
    geom_textpath(mapping = aes(x = (min_+max_)/2, y = .65, label = pr, color = nm), size = 3, show.legend = FALSE),
    geom_textpath(mapping = aes(x = (min_+(min_ + Simon_maxResp))/2, y = .95, label = pr_simon, color = nm), size = 3, show.legend = FALSE),
    geom_textpath(mapping = aes(x = (min_+max_)/2, y = 1.1, label = nm1, color = nm), show.legend = FALSE),
    scale_fill_discrete(breaks = nm, labels = legd, name = NULL),
    ylim(c(0,1.2)),
    coord_polar(theta = 'x'),
    labs(caption = sprintf('Based on %d simulations', R))
  )
}




#' @importFrom ggplot2 autoplot ggplot theme theme_void
#' @export
autoplot.Simon_oc <- function(object, ...) {
  ggplot() + autolayer.Simon_oc(object, ...) + 
    theme_void() +
    theme(
      legend.position = 'inside'
    )
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
#' @keywords internal
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
