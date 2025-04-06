

#' @title S4 Object \linkS4class{ph2simon4}
#' 
#' @description
#' One or more Simon's two-stage designs.
#' 
#' @slot r1,r non-negative \link[base]{integer} scalars or \link[base]{vector}s
#' 
#' @slot n1,n positive \link[base]{integer} scalars or \link[base]{vector}s
#' 
#' @slot pu,pa \link[base]{double} scalars
#' 
#' @slot alpha,beta \link[base]{double} scalars
#' 
#' @slot nmax \link[base]{integer}
#' 
#' @slot type \link[base]{character} scalars or \link[base]{vector}s,
#' type of Simon's two-stage design, 
#' 
#' @keywords internal
#' @name ph2simon4
#' @aliases ph2simon4-class
#' @export
setClass(Class = 'ph2simon4', slots = c(
  r1 = 'integer',
  n1 = 'integer',
  r = 'integer',
  n = 'integer',
  pu = 'numeric', pa = 'numeric',
  alpha = 'numeric', beta = 'numeric',
  nmax = 'integer',
  type = 'character'
))




#' @title Probabilities of one Simon's Two-Stage Design
#' 
#' @description 
#' Probabilities of frail 
#' (i.e., early termination) and success (to reject \eqn{H_0}) 
#' of **one** Simon's two-stage design, at given true response rate(s).
#' 
#' @slot frail \link[base]{numeric} scalar or \link[base]{vector},
#' probabilities of frail (i.e., early termination) 
#' at given true response rate(s) \eqn{p}.
#' 
#' @slot reject \link[base]{numeric} scalar or \link[base]{vector},
#' probabilities of success (to reject \eqn{H_0}) 
#' at given true response rate(s) \eqn{p}.
#' 
#' @slot eN \link[base]{numeric} scalar or \link[base]{vector}, 
#' expected sample size(s) \eqn{\textrm{E}(n)} 
#' at given true response rate(s) \eqn{p}.
#' 
#' @slot prob \link[base]{double} scalar or \link[base]{vector}, true response rate(s) \eqn{p}
#' 
#' @keywords internal
#' @name simon_pr
#' @aliases simon_pr-class
#' @export
setClass(Class = 'simon_pr', slots = c(
  frail = 'numeric',
  reject = 'numeric',
  eN = 'numeric', 
  prob = 'numeric'
))





#' @title Operating Characteristics of Simon's Two-Stage Design
#' 
#' @description
#' Operating characteristics of **one** Simon's two-stage design.
#' 
#' @slot maxResp \link[base]{integer} \link[base]{vector} of same length as \eqn{p}, 
#' the frequencies of each regime having maximum response.  
#' The summation of `maxResp` is the number of simulation copies.
#' 
#' @slot simon_maxResp \link[base]{integer} \link[base]{vector} of same length as \eqn{p}, 
#' the frequencies of each regime having maximum response and success in Simon's two-stage trial.
#' 
#' @keywords internal
#' @name simon_oc
#' @aliases simon_oc-class 
#' @export
setClass(Class = 'simon_oc', contains = 'simon_pr', slots = c(
  maxResp = 'integer',
  simon_maxResp = 'integer'
))

