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
  
  z <- model |> 
    summary.ph2simon4(...)
  
  # wording copied from
  # https://cancer.unc.edu/biostatistics/program/ivanova/SimonsTwoStageDesign.aspx
  
  sprintf(
    fmt = 'Simon\'s %s two-stage design for testing the null hypothesis p\u207A\u2264%.f%% with type-I-error rate %.1f%%, as described below, achieves %.1f%% power at true p\u207A=%.0f%%. In the first stage, %d patients will be accrued. The trial will be terminated if %d or fewer patients respond. Otherwise, another %d patients will be accrued for a total of %d patients, and the `drug` will be rejected if the total number of responses is %d or fewer. This design is provided by <u>**`R`**</u> package <u>**`clinfun`**</u>.',
    model@type,  
    1e2*model@pu,
    1e2*model@alpha, 
    1e2*(1-z$p[,4L]), 1e2*model@pa,
    n1, r1,
    n - n1, n, r)
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
#' @importFrom scales label_percent
#' @export
Sprintf.simon_oc <- function(model, ...) {
  nm <- names(prob <- model@prob)
  N <- sum(model@maxResp)
  sprintf(
    fmt = 'We simulated %d trials of each of the %d drugs %s with estimated response rates of %s, respectively, using this design. The percentage of trials with each of these drugs having the highest number of responses are %s. The percentage of trials with each of these drugs both having the highest number of responses and being accepted by the Simon\'s two-stage design are %s.', 
    N, length(prob),
    paste0('`', nm, '`', collapse = ', '),
    prob |> label_percent()() |> paste(collapse = ', '),
    (model@maxResp / N) |> label_percent(accuracy = .1, suffix = sprintf('%% for %s', nm))() |> paste(collapse = ', '),
    (model@simon_maxResp / N) |> label_percent(accuracy = .1, suffix = sprintf('%% for %s', nm))() |> paste(collapse = ', ')
  )
}

