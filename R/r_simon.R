

#' @title Random Generator based on Simon's Two-Stage Design
#' 
#' @description
#' Random generator based on Simon's two-stage design.
#' 
#' @param R positive \link[base]{integer} scalar, number of trials \eqn{R}
#' 
#' @param prob \link[base]{double} scalar, true response rate \eqn{p}
#' 
#' @param object a \link[clinfun]{ph2simon} or \linkS4class{ph2simon4} object
#' 
#' @param ... parameters of function [ph2simon4()], most importantly `type`
#' 
#' @param r1,n1,r,n (optional) \link[base]{integer} scalars, see \linkS4class{ph2simon4}.
#' 
#' @details
#' Function [r_simon()] generates \eqn{R} copies of the number of responses \eqn{y} in **one** Simon's two-stage design.
#' The conclusion of the trials are, 
#' \describe{
#' \item{\eqn{y \leq r_1}}{indicates early termination}
#' \item{\eqn{r_1 < y \leq r}}{indicates failure to reject \eqn{H_0}}
#' \item{\eqn{y > r}}{indicates success to reject \eqn{H_0}}
#' }
#' 
#' Here \eqn{r} is not needed to *generate* the random number of responses \eqn{y}.
#' Instead, \eqn{r} is needed to *determine* if the trial is a failure or a success. 
#' Therefore, \eqn{r} is not a parameter of function [r_simon()].
#' 
#' @returns
#' Function [r_simon()] returns an \link[base]{integer} \link[base]{vector} of length \eqn{R},
#' which are the \eqn{R} copies of the number of responses in the Simon's two-stage design.
#' 
#' @keywords internal
#' @name r_simon
#' @export
r_simon <- function(R, prob, object, ...) UseMethod(generic = 'r_simon', object = object)

#' @rdname r_simon
#' @export r_simon.ph2simon
#' @export
r_simon.ph2simon <- function(R, prob, object, ...) object |> ph2simon4(...) |> r_simon.ph2simon4(R = R, prob = prob)
  
#' @rdname r_simon
#' @importFrom stats rbinom
#' @export r_simon.ph2simon4
#' @export
r_simon.ph2simon4 <- function(
    R, prob, 
    object, ...,
    r1 = object@r1, n1 = object@n1, r = object@r, n = object@n
) {
  
  if (length(r1) != 1L || !is.integer(r1) || is.na(r1) || r1 < 0L) stop('`r1` must be non-negative integer scalar')
  if (length(n1) != 1L || !is.integer(n1) || is.na(n1) || n1 <= 0L) stop('`n1` must be positive integer scalar')
  if (length(r) != 1L || !is.integer(r) || is.na(r) || r < 0L) stop('`r` must be non-negative integer scalar')
  if (length(n) != 1L || !is.integer(n) || is.na(n) || n <= 0L) stop('`n` must be positive integer scalar')
  
  ret <- rbinom(n = R, size = n1, prob = prob) # number of responses in stage1  
  id2 <- (ret > r1) # indices of trials going to stage2
  ret[id2] <- ret[id2] + rbinom(n = sum(id2), size = n - n1, prob = prob) # number of positive responses in stage2
  
  dx <- cut.default(ret, breaks = c(0, r1, r, n), right = TRUE, include.lowest = TRUE)
  levels(dx) <- paste(levels(dx), c('Early Termination', 'Fail', 'Success'), sep = '; ')
  attr(ret, which = 'dx') <- dx
  
  return(ret)
  
}


