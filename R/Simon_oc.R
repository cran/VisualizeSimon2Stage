


#' @title S4 class \linkS4class{Simon_oc}
#' 
#' @slot prob \strong{named} \link[base]{numeric} vector of length \eqn{p}, true response rate(s)
#' 
#' @slot maxResp \link[base]{integer} vector of length \eqn{p}, 
#' the frequencies of each regime having maximum response.  
#' The summation of \code{maxResp} is the number of simulation copies.
#' 
#' @slot Simon_maxResp \link[base]{integer} vector of length \eqn{p}, 
#' the frequencies of each regime having maximum response and success in Simon's two-stage trial.
#' 
#' @slot eN \link[base]{numeric} vector of length \eqn{p}, 
#' expected sample sizes by simulation
#' 
#' @export
setClass(Class = 'Simon_oc', slots = c(
  prob = 'numeric',
  maxResp = 'integer',
  Simon_maxResp = 'integer',
  #n1 = 'integer', n = 'integer', r1 = 'integer', r = 'integer',
  eN = 'numeric'
))



#' @title Operating Characteristics of Simon's Two-Stage Design
#' 
#' @description ..
#' 
#' @param prob \strong{named} \link[base]{numeric} vector, true response rate(s)
#' 
#' @param simon \link[clinfun]{ph2simon} object
#' 
#' @param type \link[base]{character} scalar, either
#' \code{'minimax'} for Simon's two-stage design with minimum total sample size (default), 
#' \code{'optimal'} for minimum expected total sample size \strong{under \eqn{p_0}}, 
#' \code{'n1'} for minimum stage-1 sample size,
#' or \code{'maximax'} for maximum total sample size (as provided by user).
#' 
#' @param n1,n (optional) \link[base]{integer} scalars, stage 1 sample size \eqn{n_1} 
#' and total sample size \eqn{n}.  Will be overridden if \code{simon} is given
#' 
#' @param r1,r (optional) \link[base]{integer} scalars, number of positive response
#' in Stage 1 \eqn{r_1} and overall \eqn{r} required \strong{exclusive}.
#' In other words, passing Stage 1 means observing \eqn{>r_1} positive response.
#' Will be overridden if \code{simon} is given
#' 
#' @param N \link[base]{integer} scalar, number of simulations
#' 
#' @param ... potential parameters, currently not in use
#' 
#' @details ..
#' 
#' @return 
#' 
#' \link{Simon_oc} returns \link{Simon_oc} object
#' 
#' @examples 
#' library(clinfun)
#' (x = ph2simon(pu = .2, pa = .4, ep1 = .05, ep2 = .1)) 
#' Simon_oc(prob = c(A = .3, B = .2, C = .15), simon = x, N = 1e3L)
#' 
#' @references 
#' \doi{10.1016/0197-2456(89)90015-9}
#' 
#' @export
Simon_oc <- function(
    prob, 
    simon, type = c('minimax', 'optimal', 'n1', 'maximax'),
    N,
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

  
  sim_simon <- function(prob, n1, n, r1, N, ...) {
    ret <- rbinom(n = N, size = n1, prob = prob) # number of positive responses in stage1  
    id2 <- (ret > r1) # indexes of trials going to stage2
    ret[id2] <- ret[id2] + rbinom(n = sum(id2), size = n - n1, prob = prob) # number of positive responses in stage2
    return(ret)
    # terminate: ret <= r1
    # fail: ret <= r (`r` does not need to be a parameter in this function)
    # success: ret > r
  }
  
  M <- do.call(cbind, args = lapply(prob, FUN = sim_simon, n1 = n1, n = n, r1 = r1, N = N, ...)) 
  # `N x pn` 'matrix' # number of positive responses of each regimen
  idx <- max.col(M, ties.method = 'first') # indexes of regimen being chosen (i.e., having maximum positive response) at each of the simulated trials
  idx_success <- idx[.rowSums(M > r, m = N, n = pn, na.rm = FALSE) > 0L] # faster than `rowAny(M > r)`
  # indexes of regimen {having maximum positive response} AND {succeeding in Simon trial}
  eN <- .colMeans(n1 + (M > r1) * (n - n1), m = N, n = pn, na.rm = FALSE) # expected total sample size
  
  new(Class = 'Simon_oc', 
    maxResp = tabulate(idx, nbins = pn),  #  / N
    Simon_maxResp = tabulate(idx_success, nbins = pn), # / N
    prob = prob,
    #n1 = n1, n = n, r1 = r1, r = r, 
    eN = eN
  )
}



#' @title Show \linkS4class{Simon_oc} Object
#' 
#' @description Show \linkS4class{Simon_oc} object
#' 
#' @param object \linkS4class{Simon_oc} object
#' 
#' @return 
#' The \link[methods]{show} method for \linkS4class{Simon_oc} object 
#' does not have a returned value.
#' 
#' @export
setMethod(f = show, signature = signature(object = 'Simon_oc'), definition = function(object) {
  print(autoplot.Simon_oc(object))
  #cat('\n')
  #cat(.green(as.character.modelText(modelText.Simon_oc(x))), sep ='\n\n')
  #cat('\n')
})









#' @export
autolayer.Simon_oc <- function(object, ...) {
  pn <- length(prob <- object@prob)
  N <- sum(object@maxResp)
  maxResp <- object@maxResp / N
  Simon_maxResp <- object@Simon_maxResp / N
  ymax <- cumsum(maxResp)
  ymin <- c(0, ymax[-pn])
  nm <- names(prob)
  
  list(
    geom_rect(mapping = aes(ymax = ymax, ymin = ymin, xmax = 1, xmin = 0, fill = nm), colour = 'white', alpha = .3, stat = 'identity'),
    geom_rect(mapping = aes(ymax = ymin + Simon_maxResp, ymin = ymin, xmax = 1, xmin = 0, fill = nm), stat = 'identity', colour = 'white'),
    coord_polar(theta = 'y'),
    scale_fill_discrete(
      name = sprintf('Regimen\n(%d Simulated Trials)', N),
      labels = sprintf(
        fmt = '%s; p = %.f%%\nHaving Max # of Responses: %.1f%%\nHaving Max # of Responses & Simon\'s Success: %.1f%%\nExpected Sample Size: %.1f\n', 
        nm, 1e2*prob, 1e2*maxResp, 1e2*Simon_maxResp, object@eN))
  )
}

#' @export
autoplot.Simon_oc <- function(object, ...) {
  ggplot() + autolayer.Simon_oc(object, ...) + theme_void()
}








