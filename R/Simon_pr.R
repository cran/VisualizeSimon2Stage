


#' @title S4 class \linkS4class{Simon_pr}
#' 
#' @slot .Data \eqn{p\times 3} \link[base]{numeric} \link[base]{matrix}, probability of frail 
#' (i.e., early termination), fail (to reject the null) and success (to reject the null), at each
#' true response rate given in \code{@@prob}
#' 
#' @slot eN \link[base]{numeric} vector of length \eqn{p}, expected sample size(s)
#' 
#' @slot prob \link[base]{numeric} vector of length \eqn{p}, true response rate(s)
#' 
# @slot n1,n \link[base]{integer} scalars, Stage 1 sample size \eqn{n_1} and total sample size \eqn{n}
# 
# @slot r1,r \link[base]{integer} scalars, number of positive response
# in Stage 1 \eqn{r_1} and overall \eqn{r} required \strong{exclusive}.
# In other words, passing Stage 1 means observing \eqn{>r_1} positive response.
#' 
#' @export
setClass(Class = 'Simon_pr', contains = 'matrix', slots = c( # was 'Simon_pr'
  eN = 'numeric', 
  prob = 'numeric'#,
  #n1 = 'integer',
  #n = 'integer',
  #r1 = 'integer',
  #r = 'integer'
))




#' @title Probabilities of Simon's Two-Stage Design
#' 
#' @description 
#' Probability of frail 
#' (i.e., early termination), fail (to reject the null) and success (to reject the null) 
#' of a Simon's Two-Stage Design, at given true response rate(s).
#' 
#' @param n1,n \link[base]{integer} scalars, Stage 1 sample size \eqn{n_1} and total sample size \eqn{n}
#' 
#' @param r1,r \link[base]{integer} scalars, number of positive response
#' in Stage 1 \eqn{r_1} and overall \eqn{r} required \strong{exclusive}.
#' In other words, passing Stage 1 means observing \eqn{>r_1} positive response.
#' 
#' @param prob \link[base]{numeric} vector, true response rate(s)
#' 
#' @details
#' Parameters nomenclature of \code{n1}, \code{n}, \code{r1} and \code{r} follows that of 
#' PASS and \link[clinfun]{ph2simon}.
#' 
#' @return 
#' 
#' \link{Simon_pr} returns \linkS4class{Simon_pr} object.
#' 
#' @examples 
#' Simon_pr(n1 = 15L, r1 = 3L, n = 24L, r = 7L, prob = c(.2, .3))
#' 
#' @references 
#' \doi{10.1016/0197-2456(89)90015-9}
#' 
#' @export
Simon_pr <- function(prob, n1, n, r1, r) {
  if (length(n1) != 1L || !is.integer(n1) || n1 <= 0L) stop('`n1` must be positive integer scalar')
  if (length(n) != 1L || !is.integer(n) || n <= 0L) stop('`n` must be positive integer scalar')
  if (length(r1) != 1L || !is.integer(r1) || r1 <= 0L) stop('`r1` must be positive integer scalar')
  if (length(r) != 1L || !is.integer(r) || r <= 0L) stop('`r` must be positive integer scalar')
  
  if (!length(prob) || !is.double(prob) || anyNA(prob) || any(prob < 0, prob > 1)) stop('`prob` must be (0,1) vector')
  
  n2 <- n - n1 # Stage 2 sample size
  s1 <- (r1 + 1L):n1 # responses needed in Stage 1 to continue
  s2 <- r - s1 # `s2 + 1L` is the responses needed in Stage 2 to succeed # ?stats::pbinom ok with negative `q` 
  
  tmp <- lapply(prob, FUN = function(pp) {
    d1s <- dbinom(x = s1, size = n1, prob = pp)
    p_success <- sum(d1s * pbinom(q = s2, size = n2, prob = pp, lower.tail = FALSE))
    p_frail <- pbinom(q = r1, size = n1, prob = pp, lower.tail = TRUE) # early termination
    p_fail <- 1 - p_frail - p_success
    # stopifnot(isTRUE(all.equal.numeric(p_fail, sum(d1s * pbinom(q = s2, size = n2, prob = pp, lower.tail = TRUE)))))
    unname(c(p_frail, p_fail, p_success))
  })
  ret <- do.call(rbind, args = tmp)
  new('Simon_pr', ret, 
      eN = ret[,1L] * n1 + (1 - ret[,1L]) * n, 
      prob = prob#,
      #n1 = n1, n = n, r1 = r1, r = r # R error..
      )
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

