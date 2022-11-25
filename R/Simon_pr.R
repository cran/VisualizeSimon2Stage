

#' @title S4 Class \linkS4class{simon_design}
#' 
#' @slot n1,n \link[base]{integer} scalars, stage-1 sample size \eqn{n_1} and total sample size \eqn{n}
#' 
#' @slot r1,r \link[base]{integer} scalars, number of positive response
#' in stage-1 \eqn{r_1} and overall \eqn{r} required \strong{exclusive}.
#' In other words, passing Stage 1 indicates observing \eqn{>r_1} positive responses,
#' and rejecting \eqn{H_0} indicates observing \eqn{>r} positive responses.
#' 
#' @export
setClass(Class = 'simon_design', slots = c( # using name 'simon' will have `initMatrix` error
  r1 = 'integer', n1 = 'integer', r = 'integer', n = 'integer'
), validity = function(object) {
  r1 <- object@r1
  n1 <- object@n1
  r <- object@r
  n <- object@n
  x <- c(r1, n1, r, n)
  if (anyNA(x) || any(x < 0L)) stop('illegal r1, n1, r, n')
  if ((r1>n1) || (r>n) || (r1>r) || (n1>n)) stop('illegal r1, n1, r, n')
})



#' @title S4 Class \linkS4class{Simon_pr}
#' 
#' @description 
#' Probabilities of early termination, failure and success, of Simon's Two-Stage Design.
#' 
#' @slot .Data \eqn{l\times 3} \link[base]{numeric} \link[base]{matrix}, probability of frail 
#' (i.e., early termination), fail (to reject the null) and success (to reject the null), at each
#' response rate \eqn{p} given in \code{@@prob}
#' 
#' @slot eN \link[base]{numeric} \link[base]{vector} of length \eqn{l}, expected sample size(s) \eqn{\textrm{E}(N)}
#' 
#' @slot prob \link[base]{numeric} \link[base]{vector} of length \eqn{l}, response rate(s) \eqn{p}
#' 
#' @slot simon \linkS4class{simon_design} object
#' 
#' @export
setClass(Class = 'Simon_pr', contains = 'matrix', slots = c(
  eN = 'numeric', 
  prob = 'numeric',
  simon = 'simon_design'
))




#' @title S4 Class \linkS4class{Simon_oc}
#' 
#' @slot .Data \linkS4class{Simon_pr} object
#' 
#' @slot maxResp \link[base]{integer} \link[base]{vector} of length \eqn{p}, 
#' the frequencies of each regime having maximum response.  
#' The summation of \code{maxResp} is the number of simulation copies.
#' 
#' @slot Simon_maxResp \link[base]{integer} \link[base]{vector} of length \eqn{p}, 
#' the frequencies of each regime having maximum response and success in Simon's two-stage trial.
#' 
#' @export
setClass(Class = 'Simon_oc', contains = 'Simon_pr', slots = c(
  maxResp = 'integer',
  Simon_maxResp = 'integer'
))











#' @title Probabilities of Simon's Two-Stage Design
#' 
#' @description 
#' Probability of frail 
#' (i.e., early termination), fail (to reject the null) and success (to reject the null) 
#' of a Simon's two-stage design, at given true response rate(s).
#' 
#' @param n1,n \link[base]{integer} scalars, stage-1 sample size \eqn{n_1} and total sample size \eqn{n}
#' 
#' @param r1,r \link[base]{integer} scalars, number of positive response
#' in stage-1 \eqn{r_1} and overall \eqn{r} required \strong{exclusive}.
#' In other words, passing Stage 1 indicates observing \eqn{>r_1} positive responses,
#' and rejecting \eqn{H_0} indicates observing \eqn{>r} positive responses.
#' 
#' @param prob \link[base]{numeric} \link[base]{vector}, true response rate(s) \eqn{p}
#' 
#' @details
#' 
#' Given the Simon's two-stage design with \eqn{n_1}, \eqn{r_1}, \eqn{n} and \eqn{r}, for a response rate
#' \eqn{p}, we have the number of Stage-1 positive responses \eqn{X_1 \sim \textrm{Binom}(n_1, p)} 
#' and the number of Stage-2 positive responses \eqn{X_2 \sim \textrm{Binom}(n-n_1, p)}.  
#' Obviously \eqn{X_1} and \eqn{X_2} are independent.
#' 
#' The probability of early termination is \eqn{\textrm{Pr}(X_1 \leq r_1)}.
#' 
#' The probability of failure to reject \eqn{H_0} is 
#' \deqn{\sum_{s_1 = r_1+1}^{n_1} \textrm{Pr}(X_1=s_1)\cdot\textrm{Pr}(X_2 \leq (r-s_1))}
#' 
#' The probability of rejection of \eqn{H_0} is 
#' \deqn{\sum_{s_1 = r_1+1}^{n_1} \textrm{Pr}(X_1=s_1)\cdot\textrm{Pr}(X_2 > (r-s_1))}
#' 
#' Parameters nomenclature of \code{n1}, \code{n}, \code{r1} and \code{r} follows that of 
#' PASS and \link[clinfun]{ph2simon}.
#' 
#' @return 
#' 
#' \link{Simon_pr} returns \linkS4class{Simon_pr} object.
#' 
#' @references 
#' \doi{10.1016/0197-2456(89)90015-9}
#' 
#' @importFrom stats dbinom pbinom
#' 
#' @examples 
#' Simon_pr(prob = c(.2, .4), n1 = 15L, r1 = 3L, n = 24L, r = 7L)
#' 
#' @export
Simon_pr <- function(prob, n1, n, r1, r) {
  if (length(n1) != 1L || !is.integer(n1) || n1 <= 0L) stop('`n1` must be positive integer scalar')
  if (length(n) != 1L || !is.integer(n) || n <= 0L) stop('`n` must be positive integer scalar')
  if (length(r1) != 1L || !is.integer(r1) || r1 < 0L) stop('`r1` must be positive integer scalar')
  if (length(r) != 1L || !is.integer(r) || r < 0L) stop('`r` must be positive integer scalar')
  
  if (!length(prob) || !is.double(prob) || anyNA(prob) || any(prob < 0, prob > 1)) stop('`prob` must be (0,1) vector')
  
  s1 <- (r1 + 1L):n1 # responses needed in Stage 1 to continue
  
  ret <- do.call(rbind, args = lapply(prob, FUN = function(p) {
    
    # `p`: response rate
    
    p_frail <- pbinom(q = r1, size = n1, prob = p, lower.tail = TRUE) 
    # early termination; Prob(X1 <= r1), where X1 ~ Binom(n1, p)
    
    d1s <- dbinom(x = s1, size = n1, prob = p)
    p_success <- sum(d1s * pbinom(q = r - s1, size = n - n1, prob = p, lower.tail = FALSE))
    # success in rejection H0; \sum_{s1 = r1+1}^{n1} Prob(X1 = s1) * Prob(X2 > (r-s1))
    # .. where X1 ~ Binom(n1, p), X2 ~ Binom(n-n1, p); X1 and X2 are independent
    # ?stats::pbinom ok with negative `q`
    
    p_fail <- 1 - p_frail - p_success
    # stopifnot(isTRUE(all.equal.numeric(p_fail, sum(d1s * pbinom(q = r - s1, size = n - n1, prob = p, lower.tail = TRUE)))))
    # failure to reject H0; \sum_{s1 = r1+1}^{n1} Prob(X1 = s1) * Prob(X2 <= (r-s1))
    # .. where X1 ~ Binom(n1, p), X2 ~ Binom(n-n1, p); X1 and X2 are independent
    
    return(unname(c(p_frail, p_fail, p_success)))
    
  }))
  
  colnames(ret) <- c('P(ET)', 'P(Fail)', 'P(Success)')
  
  new('Simon_pr', ret, 
      eN = ret[,1L] * n1 + (1 - ret[,1L]) * n, 
      prob = prob,
      simon = new(Class = 'simon_design', n1 = n1, n = n, r1 = r1, r = r))
}


#' @title Show \linkS4class{Simon_pr} Object
#' 
#' @description Show \linkS4class{Simon_pr} object
#' 
#' @param object \linkS4class{Simon_pr} object
#' 
#' @return 
#' The \link[methods]{show} method for \linkS4class{Simon_pr} object 
#' does not have a returned value.
#' 
#' @export
setMethod(f = show, signature = signature(object = 'Simon_pr'), definition = function(object) {
  .data <- object@.Data
  .data[] <- sprintf('%.1f%%', 1e2*.data)
  dimnames(.data) <- list('Response Rate & E(N)' = sprintf('%.0f%%; %.1f', 1e2*object@prob, object@eN), 
                          'Probabilities' = c('Early Termination', 'Fail', 'Success'))
  print(.data, right = TRUE, quote = FALSE)
  print(autoplot.Simon_pr(object))
})


#' @export
`[.Simon_pr` <- function(x, i) {
  x@.Data <- `[`(unclass(x), i, j = TRUE, drop = FALSE)
  x@eN <- `[`(x@eN, i = i)
  x@prob <- `[`(x@prob, i = i)
  new(Class = 'Simon_pr', x) 
}



#' @export
autolayer.Simon_pr <- function(object, ...) {
  object <- object[1L] # `[.Simon_pr`
  obj <- unclass(object)
  nm <- paste(c('Early Termination', 'Fail', 'Success'), sprintf(fmt = '%.1f%%', 1e2*obj))
  ymax <- cumsum(obj)
  ymin <- c(0, ymax[-3L])
  list(
    geom_rect(mapping = aes(ymax = ymax, ymin = ymin, xmax = 1, xmin = 0, fill = nm), colour = 'white'),
    coord_polar(theta = 'y'),
    scale_fill_discrete(name = sprintf(fmt = 'Simon\'s 2-Stage\nResponse Rate %.0f%%\nExpected Total # = %.1f', 1e2*object@prob, object@eN))
  )
}

#' @export
autoplot.Simon_pr <- function(object, ...) {
  ggplot() + autolayer.Simon_pr(object, ...) + theme_void()
}













#' @title Operating Characteristics of Simon's Two-Stage Design
#' 
#' @description ..
#' 
#' @param prob \strong{named} \link[base]{numeric} \link[base]{vector}, true response rate(s)
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
#' @param R \link[base]{integer} scalar, number of simulations
#' 
#' @param ... potential parameters, currently not in use
#' 
#' @details ..
#' 
#' @return 
#' 
#' \link{Simon_oc} returns \link{Simon_oc} object
#' 
#' @references 
#' \doi{10.1016/0197-2456(89)90015-9}
#' 
#' @importFrom stats rbinom
#' 
#' @examples 
#' library(clinfun)
#' (x = ph2simon(pu = .2, pa = .4, ep1 = .05, ep2 = .1)) 
#' Simon_oc(prob = c(A = .3, B = .2, C = .15), simon = x, type = 'minimax', R = 1e3L)
#' Simon_oc(prob = c(A = .3, B = .2, C = .15), simon = x, type = 'optimal', R = 1e3L)
#' 
#' @export
Simon_oc <- function(
    prob, 
    simon, type = c('minimax', 'optimal', 'n1', 'maximax'),
    R,
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
  
  
  sim_simon <- function(prob, n1, n, r1, R, ...) {
    ret <- rbinom(n = R, size = n1, prob = prob) # number of positive responses in stage1  
    id2 <- (ret > r1) # indexes of trials going to stage2
    ret[id2] <- ret[id2] + rbinom(n = sum(id2), size = n - n1, prob = prob) # number of positive responses in stage2
    return(ret)
    # terminate: ret <= r1
    # fail: ret <= r (`r` does not need to be a parameter in this function)
    # success: ret > r
  }
  
  M <- do.call(cbind, args = lapply(prob, FUN = sim_simon, n1 = n1, n = n, r1 = r1, R = R, ...)) 
  # `R x pn` 'matrix' # number of positive responses of each regimen
  idx <- max.col(M, ties.method = 'first') # indexes of regimen being chosen (i.e., having maximum positive response) at each of the simulated trials
  idx_success <- idx[.rowSums(M > r, m = R, n = pn, na.rm = FALSE) > 0L] # faster than `Rfast::rowAny(M > r)`
  # indexes of regimen {having maximum positive response} AND {succeeding in Simon trial}
  
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
  R <- sum(object@maxResp)
  maxResp <- object@maxResp
  Simon_maxResp <- object@Simon_maxResp
  ymax <- cumsum(maxResp)
  ymin <- c(0, ymax[-pn])
  nm <- names(prob)
  
  list(
    geom_rect(mapping = aes(ymax = ymax, ymin = ymin, xmax = 1, xmin = 0, fill = nm), alpha = .3, stat = 'identity', colour = 'white'),
    geom_rect(mapping = aes(ymax = ymin + Simon_maxResp, ymin = ymin, xmax = 1, xmin = 0, fill = nm), stat = 'identity', colour = 'white'),
    coord_polar(theta = 'y'),
    scale_fill_discrete(
      name = sprintf('Regimen\n(%d Simulated Trials)', R),
      labels = sprintf(
        fmt = '%s; p = %.f%%\nHaving Max # of Responses: %.1f%%\nHaving Max # of Responses & Simon\'s Success: %.1f%%\nExpected Sample Size: %.1f\n', 
        nm, 1e2*prob, 1e2*maxResp/R, 1e2*Simon_maxResp/R, object@eN))
  )
}

#' @export
autoplot.Simon_oc <- function(object, ...) {
  ggplot() + autolayer.Simon_oc(object, ...) + theme_void()
}










