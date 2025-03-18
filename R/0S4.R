

#' @title S4 Object \linkS4class{ph2simon4}
#' 
#' @description
#' One or more Simon's two-stage designs.
#' 
#' @slot r1,r non-negative \link[base]{integer} scalars or \link[base]{vector}s, 
#' number of responses
#' in Stage-1 \eqn{r_1} and overall \eqn{r} required *exclusively*, in other words
#' \itemize{
#' \item {pass Stage-1 if observed \eqn{>r_1} response;}
#' \item {reject \eqn{H_0} if observed \eqn{>r} responses.}
#' }
#' 
#' @slot n1,n positive \link[base]{integer} scalars or \link[base]{vector}s, 
#' Stage-1 sample size \eqn{n_1} 
#' and total sample size \eqn{n}.
#' 
#' @slot pu,pa \link[base]{double} scalars
#' 
#' @slot alpha,beta \link[base]{double} scalars
#' 
#' @slot nmax \link[base]{integer}
#' 
#' @slot type \link[base]{character} scalars or \link[base]{vector}s,
#' type of Simon's two-stage design, one or more values among
#' \describe{
#' \item{`'minimax'`}{(default) minimum total sample size}
#' \item{`'optimal'`}{minimum expected total sample size *under \eqn{p_0}*}
#' \item{`'n1'`}{minimum Stage-1 sample size}
#' \item{`'maximax'`}{to use up the user-provided maximum total sample size, 
#' i.e., parameter `nmax` of function \link[clinfun]{ph2simon}}
#' }
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
#' @details
#' 
#' Given one Simon's two-stage design \eqn{(r_1,n_1,r,n)} 
#' and a true response rate \eqn{p}, 
#' we have the number of Stage-1 positive responses \eqn{X_1 \sim \textrm{Binom}(n_1, p)} 
#' and the number of Stage-2 positive responses \eqn{X_2 \sim \textrm{Binom}(n-n_1, p)}.  
#' Obviously \eqn{X_1} and \eqn{X_2} are independent.
#' 
#' The probability of early termination is 
#' \deqn{p_{\textrm{frail}} = \textrm{Pr}(X_1 \leq r_1)}
#' 
#' The probability of failure to reject \eqn{H_0} is 
#' \deqn{\sum_{s_1 = r_1+1}^{n_1} \textrm{Pr}(X_1=s_1)\cdot\textrm{Pr}(X_2 \leq (r-s_1))}
#' 
#' The probability of successfully rejecting \eqn{H_0} is 
#' \deqn{\sum_{s_1 = r_1+1}^{n_1} \textrm{Pr}(X_1=s_1)\cdot\textrm{Pr}(X_2 > (r-s_1))}
#' 
#' The expected sample size is
#' \deqn{\textrm{E}(n) = p_{\textrm{frail}} \cdot n_1 + (1 - p_{\textrm{frail}}) \cdot n}
#' 
#' Parameters nomenclature of `r1`, `n1`, `r` and `n` follows that of 
#' PASS and function \link[clinfun]{ph2simon}.
#' 
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
#' @name simon_oc
#' @aliases simon_oc-class 
#' @export
setClass(Class = 'simon_oc', contains = 'simon_pr', slots = c(
  maxResp = 'integer',
  simon_maxResp = 'integer'
))

