



#' @title \linkS4class{Simon_pr}: Probabilities of a Simon's Two-Stage Design
#' 
#' @description 
#' Probability of frail 
#' (i.e., early termination), fail (to reject the null) and success (to reject the null) 
#' of a Simon's two-stage design, at given true response rate(s).
#' 
#' @slot .Data \eqn{l\times 3} \link[base]{numeric} \link[base]{matrix}, probability of frail 
#' (i.e., early termination), fail (to reject the null) and success (to reject the null), at each
#' response rate \eqn{p} given in `@@prob`
#' 
#' @slot eN \link[base]{numeric} \link[base]{vector} of length \eqn{l}, expected sample size(s) \eqn{\textrm{E}(N)}
#' 
#' @slot prob \link[base]{numeric} \link[base]{vector} of length \eqn{l}, response rate(s) \eqn{p}
#' 
#' @name Simon_pr
#' @importFrom methods setClass
#' @aliases Simon_pr-class
#' @export
setClass(Class = 'Simon_pr', contains = 'matrix', slots = c(
  eN = 'numeric', 
  prob = 'numeric'
))







#' @rdname Simon_pr
#' 
#' @param n1,n positive \link[base]{integer} scalars, Stage-1 sample size \eqn{n_1} and total sample size \eqn{n}
#' 
#' @param r1,r non-negative \link[base]{integer} scalars, number of response
#' in Stage-1 \eqn{r_1} and overall \eqn{r} required *exclusively*.
#' In other words, passing Stage-1 indicates observing \eqn{>r_1} responses,
#' and rejecting \eqn{H_0} indicates observing \eqn{>r} responses.
#' 
#' @param prob \link[base]{numeric} \link[base]{vector}, true response rate(s) \eqn{p}
#' 
#' @details
#' 
#' Given the Simon's two-stage design \eqn{(n_1, r_1, n, r)}, for a response rate
#' \eqn{p}, we have the number of Stage-1 positive responses \eqn{X_1 \sim \textrm{Binom}(n_1, p)} 
#' and the number of Stage-2 positive responses \eqn{X_2 \sim \textrm{Binom}(n-n_1, p)}.  
#' Obviously \eqn{X_1} and \eqn{X_2} are independent.
#' 
#' The probability of early termination is \eqn{\textrm{Pr}(X_1 \leq r_1)}.
#' 
#' The probability of failure to reject \eqn{H_0} is 
#' \deqn{\sum_{s_1 = r_1+1}^{n_1} \textrm{Pr}(X_1=s_1)\cdot\textrm{Pr}(X_2 \leq (r-s_1))}
#' 
#' The probability of rejecting \eqn{H_0} is 
#' \deqn{\sum_{s_1 = r_1+1}^{n_1} \textrm{Pr}(X_1=s_1)\cdot\textrm{Pr}(X_2 > (r-s_1))}
#' 
#' Parameters nomenclature of `n1`, `n`, `r1` and `r` follows that of 
#' PASS and \link[clinfun]{ph2simon}.
#' 
#' @returns 
#' 
#' [Simon_pr] returns \linkS4class{Simon_pr} object.
#' 
#' @references 
#' \doi{10.1016/0197-2456(89)90015-9}
#' 
#' \url{https://www.ncss.com/software/pass/}
#' 
#' @examples 
#' Simon_pr(prob = c(.2, .4), n1 = 15L, r1 = 3L, n = 24L, r = 7L)
#' 
#' @importFrom methods new
#' @importFrom stats dbinom pbinom
#' @export
Simon_pr <- function(prob, n1, n, r1, r) {
  if (length(n1) != 1L || !is.integer(n1) || is.na(n1) || n1 <= 0L) stop('`n1` must be positive integer scalar')
  if (length(n) != 1L || !is.integer(n) || is.na(n) || n <= 0L) stop('`n` must be positive integer scalar')
  if (length(r1) != 1L || !is.integer(r1) || is.na(r1) || r1 < 0L) stop('`r1` must be non-negative integer scalar')
  if (length(r) != 1L || !is.integer(r) || is.na(r) || r < 0L) stop('`r` must be non-negative integer scalar')
  
  if (!length(prob) || !is.double(prob) || anyNA(prob) || any(prob < 0, prob > 1)) stop('`prob` must be (0,1) vector')
  
  s1 <- (r1 + 1L):n1 # responses needed in Stage-1 to continue
  
  ret <- do.call(rbind, args = lapply(prob, FUN = function(p) {
    
    # `p`: response rate
    
    p_frail <- pbinom(q = r1, size = n1, prob = p, lower.tail = TRUE) 
    # early termination; Prob(X1 <= r1), where X1 ~ Binom(n1, p)
    
    d1s <- dbinom(x = s1, size = n1, prob = p)
    p_success <- sum(d1s * pbinom(q = r - s1, size = n - n1, prob = p, lower.tail = FALSE))
    # success in rejection H0; \sum_{s1 = r1+1}^{n1} Prob(X1 = s1) * Prob(X2 > (r-s1))
    # .. where X1 ~ Binom(n1, p), X2 ~ Binom(n-n1, p); X1 and X2 are independent
    # ?stats::pbinom is okay with negative `q`
    
    p_fail <- 1 - p_frail - p_success
    # stopifnot(isTRUE(all.equal.numeric(p_fail, sum(d1s * pbinom(q = r - s1, size = n - n1, prob = p, lower.tail = TRUE)))))
    # failure to reject H0; \sum_{s1 = r1+1}^{n1} Prob(X1 = s1) * Prob(X2 <= (r-s1))
    # .. where X1 ~ Binom(n1, p), X2 ~ Binom(n-n1, p); X1 and X2 are independent
    
    return(unname(c(p_frail, p_fail, p_success)))
    
  }))
  
  colnames(ret) <- c('P(ET)', 'P(Fail)', 'P(Success)')
  
  new(Class = 'Simon_pr', ret, 
      eN = ret[,1L] * n1 + (1 - ret[,1L]) * n, 
      prob = prob)
}







#' @title Show \linkS4class{Simon_pr} Object
#' 
#' @description Show \linkS4class{Simon_pr} object
#' 
#' @param object \linkS4class{Simon_pr} object
#' 
#' @returns 
#' The \link[methods]{show} method for \linkS4class{Simon_pr} object 
#' does not have a returned value.
#' 
#' @importFrom methods setMethod show signature
#' @export
setMethod(f = show, signature = signature(object = 'Simon_pr'), definition = function(object) {
  .data <- object@.Data
  .data[] <- sprintf('%.1f%%', 1e2*.data)
  dimnames(.data) <- list('Response Rate & E(N)' = sprintf('%.0f%%; %.1f', 1e2*object@prob, object@eN), 
                          'Probabilities' = c('Early Termination', 'Fail', 'Success'))
  print(.data, right = TRUE, quote = FALSE)
  print(autoplot.Simon_pr(object))
})








#' @importFrom methods new
#' @export
`[.Simon_pr` <- function(x, i) {
  x@.Data <- `[`(unclass(x), i, j = TRUE, drop = FALSE)
  x@eN <- `[`(x@eN, i = i)
  x@prob <- `[`(x@prob, i = i)
  new(Class = 'Simon_pr', x) 
}





#' @importFrom ggplot2 autolayer aes geom_rect coord_polar labs
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
    labs(fill = sprintf(fmt = 'Simon\'s 2-Stage\nResponse Rate %.0f%%\nExpected Total # = %.1f', 1e2*object@prob, object@eN))
  )
}




#' @importFrom ggplot2 autoplot ggplot theme_void
#' @export
autoplot.Simon_pr <- function(object, ...) {
  ggplot() + autolayer.Simon_pr(object, ...) + theme_void()
}







#' @title Random Generator based on Simon's Two-Stage Design
#' 
#' @description
#' Random generator based on Simon's two-stage design.
#' 
#' @param R positive \link[base]{integer} scalar, number of trials \eqn{R}
#' 
#' @param n1,n positive \link[base]{integer} scalars, Stage-1 sample size \eqn{n_1} and total sample size \eqn{n}
#' 
#' @param r1 non-negative \link[base]{integer} scalar, number of response
#' in Stage-1 \eqn{r_1} required *exclusively*.
#' In other words, passing Stage-1 indicates observing \eqn{>r_1} responses
#' 
#' @param prob \link[base]{numeric} scalar, true response rate \eqn{p}
#' 
#' @details
#' Function [r_simon()] generates \eqn{R} copies of the number of responses \eqn{y} in the Simon's two-stage design.
#' The conclusion of the trials are, 
#' \describe{
#' \item{\eqn{y \leq r_1}}{indicates early termination}
#' \item{\eqn{r_1 < y \leq r}}{indicates failure to reject \eqn{H_0}}
#' \item{\eqn{y > r}}{indicates success to reject \eqn{H_0}}
#' }
#' 
#' Here \eqn{r} is not needed to *generate* the random number of responses \eqn{y}.
#' Instead, \eqn{r} is needed to *determine* if the trial is a failure or a success. 
#' Therefore, \eqn{r} is not a parameter in [r_simon].
#' 
#' @returns
#' Function [r_simon()] returns an \link[base]{integer} \link[base]{vector} of length \eqn{R},
#' which are the \eqn{R} copies of the number of responses in the Simon's two-stage design.
#' 
#' @examples
#' library(clinfun)
#' ph2simon(pu = .2, pa = .4, ep1 = .05, ep2 = .1) # using 'Optimal'
#' # set.seed if needed 
#' (ys = r_simon(R = 10L, n1 = 19L, n = 54L, r1 = 4L, prob = .3))
#' cut.default(ys, breaks = c(0, 4L, 15L, 54L), right = TRUE,
#'   labels = c('early-termination', 'fail', 'success'))
#' 
#' @importFrom stats rbinom
#' @export
r_simon <- function(R, n1, n, r1, prob) {
  ret <- rbinom(n = R, size = n1, prob = prob) # number of positive responses in stage1  
  id2 <- (ret > r1) # indices of trials going to stage2
  ret[id2] <- ret[id2] + rbinom(n = sum(id2), size = n - n1, prob = prob) # number of positive responses in stage2
  return(ret)
}





