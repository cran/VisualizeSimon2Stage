
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
#' @name autoplot.ph2simon
#' @export autoplot.ph2simon
#' @export
autoplot.ph2simon <- function(object, ...) {
  ggplot() + autolayer.ph2simon(object, ...) + 
    theme_void() +
    theme(
      legend.position = 'inside'
    )
}



#' @importFrom ggplot2 autolayer aes geom_rect coord_polar ylim labs
#' @importFrom geomtextpath geom_textpath
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
  du <- sm@.Data[1L,]
  da <- sm@.Data[2L,]
  nm <- structure(1:3, levels = c('Frail', 'Fail', 'Success'), class = 'factor')
  
  max_a <- cumsum(da)
  min_a <- c(0, max_a[-3L])
  max_u <- cumsum(du)
  min_u <- c(0, max_u[-3L])
  
  list(
    geom_rect(mapping = aes(xmax = max_a, xmin = min_a, ymax = 1, ymin = .7, fill = nm), alpha = .2, colour = 'white'),
    geom_textpath(mapping = aes(x = (min_a+max_a)/2, y = .9, label = sprintf('%.1f%%', 1e2*da), color = nm), size = 3, show.legend = FALSE),
    geom_rect(mapping = aes(xmax = max_u, xmin = min_u, ymax = .65, ymin = .35, fill = nm), alpha = .2, colour = 'white', show.legend = FALSE),
    geom_textpath(mapping = aes(x = (min_u+max_u)/2, y = .55, label = sprintf('%.1f%%', 1e2*du), color = nm), size = 3, show.legend = FALSE),
    coord_polar(theta = 'x'),
    ylim(c(0,1.2)),
    labs(fill = NULL, caption = sprintf(
      fmt = 'Simon\'s 2-Stage Design\n%s\nResponse Rates: pu=%d%% vs. pa=%d%%\nE(N): %.1f vs. %.1f', 
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


