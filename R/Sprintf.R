#' @title Short Paragraph to Describe a \link[clinfun]{ph2simon} Object
#' 
#' @description
#' To create a short paragraph to describe a \link[clinfun]{ph2simon} object.
#' 
#' @param model a \link[clinfun]{ph2simon} or \linkS4class{ph2simon4} object
#' 
#' @param ... parameters of function [ph2simon4()], most importantly `type`
#' 
#' @returns 
#' Function [Sprintf.ph2simon()] and [Sprintf.ph2simon4()] return a \link[base]{character} scalar.
#' 
#' @keywords internal
#' @name Sprintf_ph2simon
#' @export
Sprintf.ph2simon4 <- function(model, ...) {
  
  r1 <- model@r1
  if (length(r1) != 1L || !is.integer(r1) || is.na(r1) || r1 < 0L) stop('`r1` must be non-negative integer scalar')
  n1 <- model@n1 
  r <- model@r
  n <- model@n
  
  msum <- summary.ph2simon4(model, ...)
  
  sprintf(
    fmt = 'Simon\'s %s two-stage design for testing the null hypothesis p\u207A\u2264%.f%% versus the alternative hypothesis p\u207A>%.f%% with type-I-error rate %.1f%%, as described below, will achieve %.1f%% power at true p\u207A=%.0f%%. The drug will be tested on %d patients in the first stage. The trial will be terminate early if %d or fewer patients respond (early termination probability %.1f%% under the null p\u207A=%.f%%). Otherwise another %d patients will be enrolled in the second stage and the drug will be rejected if the total number of patients responding is %d or fewer. This design requires a maximum sample size of %d patients, with an expected sample size of %.1f patients under the null p\u207A=%.f%%. This design is provided by <u>**`R`**</u> package <u>**`clinfun`**</u>.',
    model@type,  
    1e2*model@pu, 1e2*model@pu, 1e2*model@alpha, 1e2*(1-model@beta), 1e2*model@pa,
    n1, r1, 1e2*msum$p[,'PET(pu)'], 1e2*model@pu, 
    n - n1, r,
    n, msum$EN[,'EN(pa)'], 1e2*model@pu)
}


#' @rdname Sprintf_ph2simon
#' @export
Sprintf.ph2simon <- function(model, ...) model |> ph2simon4(...) |> Sprintf.ph2simon4()







#' @title Short Paragraph to Describe a \linkS4class{simon_oc} Object
#' 
#' @description
#' To create a short paragraph to describe a \linkS4class{simon_oc} object.
#' 
#' @param model \linkS4class{simon_oc} object
#' 
#' @param ... additional parameters, currently not in use
#' 
#' @returns 
#' Function [Sprintf.simon_oc()] returns a \link[base]{character} scalar.
#' 
#' @keywords internal
#' @export
Sprintf.simon_oc <- function(model, ...) {
  nm <- names(prob <- model@prob)
  N <- sum(model@maxResp)
  maxResp <- model@maxResp / N
  simon_maxResp <- model@simon_maxResp / N
  sprintf(
    fmt = 'We simulated %d trials of each of the %d drugs %s with estimated response rates of %s, respectively, using this design. The percentage of trials with each of these drugs having the highest number of responses are %s. The percentage of trials with each of these drugs both having the highest number of responses and being accepted by the Simon\'s two-stage design are %s.', 
    N, length(prob),
    paste(sQuote(nm), collapse = ', '),
    paste(sprintf(fmt = '%.0f%%', 1e2*prob), collapse = ', '),
    paste(sprintf(fmt = '%.1f%%', 1e2*maxResp), 'for', nm, collapse = ', '),
    paste(sprintf(fmt = '%.1f%%', 1e2*simon_maxResp), 'for', nm, collapse = ', '))
}

