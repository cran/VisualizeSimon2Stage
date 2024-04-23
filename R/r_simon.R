

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
#' in Stage-1 \eqn{r_1} required *exclusively*,
#' i.e., passing Stage-1 indicates observing \eqn{>r_1} responses
#' 
#' @param prob \link[base]{double} scalar, true response rate \eqn{p}
#' 
#' @details
#' Function [r_simon] generates \eqn{R} copies of the number of responses \eqn{y} in the Simon's two-stage design.
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
#' Function [r_simon] returns an \link[base]{integer} \link[base]{vector} of length \eqn{R},
#' which are the \eqn{R} copies of the number of responses in the Simon's two-stage design.
#' 
#' @examples
#' library(clinfun)
#' ph2simon(pu = .2, pa = .4, ep1 = .05, ep2 = .1) # using 'Optimal'
#' # set.seed if needed 
#' (ys = r_simon(R = 10L, n1 = 19L, n = 54L, r1 = 4L, prob = .3))
#' table(cut.default(ys, breaks = c(0, 4L, 15L, 54L), right = TRUE,
#'   labels = c('early-termination', 'fail', 'success')))
#' 
#' @importFrom stats rbinom
#' @export
r_simon <- function(R, n1, n, r1, prob) {
  ret <- rbinom(n = R, size = n1, prob = prob) # number of positive responses in stage1  
  id2 <- (ret > r1) # indices of trials going to stage2
  ret[id2] <- ret[id2] + rbinom(n = sum(id2), size = n - n1, prob = prob) # number of positive responses in stage2
  return(ret)
}


