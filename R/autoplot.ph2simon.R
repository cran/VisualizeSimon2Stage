
#' @title Plot Simon's Two-Stage Design
#' 
#' @description 
#' Plot \link[clinfun]{ph2simon} object using \CRANpkg{ggplot2}.
#' 
#' @param object \link[clinfun]{ph2simon} object
#' 
#' @param type \link[base]{character} scalar, one of 
#' `'minimax'`, `'optimal'`, `'n1'` and `'maximax'`
#' 
#' @param n1,n (optional) \link[base]{integer} scalars, Stage-1 sample size \eqn{n_1} 
#' and total sample size \eqn{n}.  Overridden if `object` is given
#' 
#' @param r1,r (optional) \link[base]{integer} scalars, number of response
#' in Stage-1 \eqn{r_1} and overall \eqn{r} required *exclusively*,
#' i.e., passing Stage-1 means observing \eqn{>r_1} response.
#' Overridden if `object` is given
#' 
#' @param pu,pa \link[base]{double} scalars, see function \link[clinfun]{ph2simon}
#' 
#' @param ... potential parameters, currently not in use
#' 
#' @returns 
#' Function [autoplot.ph2simon] returns a \link[ggplot2]{ggplot} object.
#' 
#' Function [autolayer.ph2simon] returns a \link[base]{list} of \link[ggplot2]{ggproto} and labels.
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
#' @importFrom ggplot2 autoplot ggplot theme theme_void
#' @importFrom grid unit
#' @name autoplot.ph2simon
#' @export autoplot.ph2simon
#' @export
autoplot.ph2simon <- function(object, ...) {
  ggplot() + autolayer.ph2simon(object, ...) + 
    theme_void() +
    theme(
      legend.key.spacing.y = unit(.01, units = 'npc')
    )
}



#' @importFrom ggplot2 autolayer aes geom_bar coord_polar xlim labs
#' @rdname autoplot.ph2simon
#' @export autolayer.ph2simon
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
  nm <- c(sprintf('Early Termination\n%.1f%% vs. %.1f%%', 1e2*dd[1L,1L], 1e2*dd[2L,1L]),
          sprintf('Fail\n%.1f%% vs. %.1f%%', 1e2*dd[1L,2L], 1e2*dd[2L,2L]), 
          #sprintf('Success\n\u03b1 = %.1f%%, 1-\u03b2 = %.1f%%', 1e2*dd[1L,3L], 1e2*dd[2L,3L]))
          sprintf('Success\nalpha = %.1f%%, power = %.1f%%', 1e2*dd[1L,3L], 1e2*dd[2L,3L]))
  
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


