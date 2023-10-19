

#' @title Alternate Print Method for a Simon's Two-Stage Design
#' 
#' @description
#' An alternate \link[base]{print} method for \link[clinfun]{ph2simon} object.
#' 
#' @param x a \link[clinfun]{ph2simon} object
#' 
#' @param ... additional parameters, currently not in use
#' 
#' @returns
#' Function [print_ph2simon()] does not have a returned value.
#' 
#' @note
#' We do not overwrite `clinfun:::print.ph2simon`.
#' 
#' @examples
#' library(clinfun)
#' (x = ph2simon(pu = .2, pa = .4, ep1 = .05, ep2 = .1)) 
#' print_ph2simon(x)
#' 
#' @export
print_ph2simon <- function(x, ...) {
  out <- summary.ph2simon(x)
  out$EN[] <- sprintf(fmt = '%.1f', out$EN)
  out$p[] <- sprintf(fmt = '%.1f%%', 1e2*out$p) 
  cat('\n Simon\'s 2-Stage Design\n\n')
  cat(sprintf('Unacceptable Response Rate: %.1f%%\n', 1e2*x[['pu']]))
  cat(sprintf('Desirable Response Rate: %.1f%%\n', 1e2*x[['pa']]))
  cat(sprintf('Controlled Error Rates: \u03b1 \u2264 %.f%%, \u03b2 \u2264 %.f%%\n', 1e2*x[['alpha']], 1e2*x[['beta']]))
  cat(sprintf('Maximum Sample Size Allowed: %d\n\n', x$nmax))
  print.default(cbind(
    out$design, 
    out$EN, 
    out$p
  ), right = TRUE, quote = FALSE, ...)
  cat('\n')
}





#' @title Summarize a Simon's Two-Stage Design
#' 
#' @description
#' Summarize a Simon's two-stage design
#' 
#' @param object \link[clinfun]{ph2simon} object
#' 
#' @param ... potential parameters, currently not in use
#' 
#' @returns
#' [summary.ph2simon] returns a \link[base]{list} with three (3) elements
#' \describe{
#' \item{`'design'`}{\link[base]{integer} \link[base]{matrix}}
#' \item{`'EN'`}{\link[base]{double} \link[base]{matrix}}
#' \item{`'p'`}{\link[base]{double} \link[base]{matrix}}
#' }
#' 
#' @examples
#' library(clinfun)
#' (x = ph2simon(pu = .2, pa = .4, ep1 = .05, ep2 = .1)) 
#' summary(x)
#' 
#' @include Simon_pr.R
#' @export summary.ph2simon
#' @export
summary.ph2simon <- function(object, ...) {
  rid <- c(
    minimax = 1L, # 'minimum total sample size'
    optimal = which.min(object$out[,'EN(p0)']), # 'minimum expected total sample size'
    n1 = which.min(object$out[,'n1']), # 'minimum Stage 1 sample size'
    maximax = which.max(object$out[,'n']) # 'maximum total sample size'
  )
  ret_design <- object$out[rid, c('r1', 'n1', 'r', 'n')]
  storage.mode(ret_design) <- 'integer'
  r1 <- ret_design[,'r1']
  n1 <- ret_design[,'n1']
  r <- ret_design[,'r']
  n <- ret_design[,'n']
  
  sm <- .mapply(FUN = function(n1, n, r1, r) {
    Simon_pr(prob = c(object$pu, object$pa), n1 = n1, n = n, r1 = r1, r = r)
  }, dots = list(n1 = n1, n = n, r1 = r1, r = r), MoreArgs = NULL)
  
  ret_EN <- do.call(rbind, args = lapply(sm, FUN = slot, name = 'eN'))
  colnames(ret_EN) <- paste0('EN(', c('pu', 'pa'), ')')
  
  ret_p <- do.call(rbind, args = lapply(sm, FUN = function(i) {
    i <- unclass(i)
    c(i[1L,1L], i[2L,1L], i[1L,3L], 1-i[2L,3L])
  }))
  colnames(ret_p) <- c('PET(pu)', 'PET(pa)', '\u03b1', '\u03b2')
  
  # In the returned value of \link[clinfun]{ph2simon}, 
  # `PET` indicates probability of early termination (i.e. frail),
  # `EN` indicates expected sample size.
  rownames(ret_design) <- rownames(ret_EN) <- rownames(ret_p) <- names(rid)
  return(list(
    design = ret_design,
    EN = ret_EN,
    p = ret_p
  ))
}



#' @importFrom ggplot2 autolayer aes geom_bar coord_polar xlim labs
#' @export
autolayer.ph2simon <- function(
    object, # return from ?clinfun::ph2simon
    type = c('minimax', 'optimal', 'n1', 'maximax'), 
    n1 = stop('must provide `n1`'), n = stop('must provide `n`'), 
    r1 = stop('must provide `r1`'), r = stop('must provide `r`'), 
    pu = stop('must provide `pu`'), pa = stop('must provide `pa`'),
    ...
) {
  
  if (!missing(object)) { # overwrite `n1`, `n`, `r1`, `r`, `pu`, `pa`
    pu <- object[['pu']]
    pa <- object[['pa']]
    type <- match.arg(type)
    obj <- summary.ph2simon(object)[['design']][type, , drop = TRUE]
    n1 <- obj['n1']
    n <- obj['n']
    r1 <- obj['r1']
    r <- obj['r']
  } else type <- '(Customized)'
  
  sm <- Simon_pr(prob = c(pu, pa), n1 = n1, n = n, r1 = r1, r = r)
  dd <- sm@.Data
  nm <- c(sprintf('Early Termination\n%.1f%% vs. %.1f%%\n', 1e2*dd[1L,1L], 1e2*dd[2L,1L]),
          sprintf('Fail\n%.1f%% vs. %.1f%%\n', 1e2*dd[1L,2L], 1e2*dd[2L,2L]), 
          #sprintf('Success\n\u03b1 = %.1f%%, 1-\u03b2 = %.1f%%\n', 1e2*dd[1L,3L], 1e2*dd[2L,3L]))
          sprintf('Success\nalpha = %.1f%%, power = %.1f%%\n', 1e2*dd[1L,3L], 1e2*dd[2L,3L]))
  
  list(
    # ?ggplot2::geom_rect wont work here
    geom_bar(mapping = aes(x = 2, y = dd[1L,], fill = nm), alpha = c(.3, .3, 1), stat = 'identity', color = 'white'),
    geom_bar(mapping = aes(x = 1, y = dd[2L,], fill = nm), alpha = c(.3, .3, 1), stat = 'identity', color = 'white'),
    coord_polar(theta = 'y', direction = -1),
    xlim(.3, 2.5),
    labs(fill = sprintf(
      fmt = 'Simon\'s 2-Stage Design\n%s\nResponse Rates: pu=%d%% vs. pa=%d%%\nExpected Total #: %.1f vs. %.1f', 
      switch(type, minimax = {
        'Minimum Total Sample Size'
      }, optimal = {
        'Minimum Expected Total Sample Size'
      }, n1 = {
        'Minimum Stage-1 Sample Size'
      }, maximax = {
        'Maximum Total Sample Size'
      }),
      1e2*pu, 1e2*pa, sm@eN[1L], sm@eN[2L]))
  )
}




#' @title Plot Simon's Two-Stage Design
#' 
#' @description 
#' Plot \link[clinfun]{ph2simon} object using \CRANpkg{ggplot2}.
#' 
#' @param object \link[clinfun]{ph2simon} object
#' 
#' @param ... potential parameters, currently not in use
#' 
#' @returns 
#' Function [autoplot.ph2simon()] returns a \link[ggplot2]{ggplot} object.
#' 
#' @examples
#' library(clinfun)
#' (x = ph2simon(pu = .2, pa = .4, ep1 = .05, ep2 = .1)) 
#' class(x)
#' autoplot(x, type = 'minimax')
#' autoplot(x, type = 'optimal')
#' autoplot(x, type = 'n1')
#' autoplot(x, type = 'maximax')
#' 
#' # example with r1 = 0
#' (des = ph2simon(pu = .05, pa = .3, ep1 = .05, ep2 = .2))
#' autoplot(des, type = 'optimal')
#' autoplot(des, type = 'minimax')
#' 
#' @importFrom ggplot2 autoplot ggplot theme_void
#' @export autoplot.ph2simon
#' @export
autoplot.ph2simon <- function(object, ...) {
  ggplot() + autolayer.ph2simon(object, ...) + theme_void()
}






#' @title Short Paragraph to Describe a \link[clinfun]{ph2simon} Object
#' 
#' @description
#' To create a short paragraph to describe a \link[clinfun]{ph2simon} object.
#' 
#' @param model \link[clinfun]{ph2simon} object
#' 
#' @param type \link[base]{character} scalar, type of Simon's two-stage design,
#' \describe{
#' \item{`'minimax'`}{(default) minimum total sample size}
#' \item{`'optimal'`}{minimum expected total sample size *under \eqn{p_0}*}
#' \item{`'n1'`}{minimum Stage-1 sample size}
#' \item{`'maximax'`}{maximum total sample size (as provided by end-user)}
#' }
#' 
#' @param ... additional parameters, currently not in use
#' 
#' @returns 
#' [Sprintf.ph2simon] returns a \link[base]{noquote} \link[base]{character} scalar.
#' 
#' @examples 
#' library(clinfun)
#' (x = ph2simon(pu = .2, pa = .4, ep1 = .05, ep2 = .1)) 
#' Sprintf.ph2simon(x, type = 'minimax')
#' Sprintf.ph2simon(x, type = 'optimal')
#' Sprintf.ph2simon(x, type = 'n1')
#' Sprintf.ph2simon(x, type = 'maximax')
#' 
#' @export Sprintf.ph2simon
#' @export
Sprintf.ph2simon <- function(model, type = c('minimax', 'optimal', 'n1', 'maximax'), ...) {
  # PASS output too long
  type <- match.arg(type)
  msum <- summary.ph2simon(model)
  design <- msum$design[type, ]
  ps <- msum$p[type, ]
  eN <- msum$EN[type, ]
  n1 <- design['n1'] 
  n <- design['n']
  noquote(sprintf(
    fmt = 'Simon\'s %s (i.e., %s) two-stage design for testing the null hypothesis p\u207A\u2264%.f%% versus the alternative hypothesis p\u207A>%.f%% with type-I-error rate %.1f%%, as described below, will achieve %.1f%% power at true p\u207A=%.0f%%. The drug will be tested on %d patients in the first stage. The trial will be terminate early if %d or fewer patients respond (early termination probability %.1f%% under the null p\u207A=%.f%%). Otherwise another %d patients will be enrolled in the second stage and the drug will be rejected if the total number of patients responding is %d or fewer. This design requires a maximum sample size of %d patients, with an expected sample size of %.1f patients under the null p\u207A=%.f%%. This design is provided by [R] package [clinfun].',
    type, switch(type, minimax = {
      'minimum total sample size'
    }, optimal = {
      'minimum expected total sample size'
    }, n1 = {
      'minimum Stage-1 sample size'
    }, maximax = {
      'maximum total sample size'
    }), 
    1e2*model$pu, 1e2*model$pu, 1e2*model$alpha, 1e2*(1-model$beta), 1e2*model$pa,
    n1, design['r1'], 1e2*ps['PET(pu)'], 1e2*model$pu, 
    n - n1, design['r'],
    n, eN['EN(pa)'], 1e2*model$pu))
}
