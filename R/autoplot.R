
#' @title Plot Simon's Two-Stage Design
#' 
#' @description 
#' Plot \link[clinfun]{ph2simon} object using \CRANpkg{ggplot2}.
#' 
#' @param object a \link[clinfun]{ph2simon} or \linkS4class{ph2simon4} object
#' 
#' @param ... parameters of function [ph2simon4()], most importantly `type`
#' 
#' @param r1,n1,r,n (optional) \link[base]{integer} scalars, see \linkS4class{ph2simon4}.
#' 
#' @param pu,pa \link[base]{double} scalars, see function \link[clinfun]{ph2simon}
#' 
#' @returns 
#' Function [autoplot.ph2simon()] returns a \link[ggplot2]{ggplot} object.
#' 
#' Function [autolayer.ph2simon()] returns a \link[base]{list} of \link[ggplot2]{ggproto} and labels.
#' 
#' @keywords internal
#' @name gg_ph2simon
#' @importFrom ggplot2 autoplot ggplot theme theme_void
#' @importFrom grid unit
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
#' @rdname gg_ph2simon
#' @export autolayer.ph2simon4
#' @export
autolayer.ph2simon4 <- function(
    object,
    r1 = object@r1, n1 = object@n1, r = object@r, n = object@n,
    pu = object@pu, pa = object@pa,
    ...
) {
  
  sm <- simon_pr.ph2simon4(prob = c(pu, pa), n1 = n1, n = n, r1 = r1, r = r)
  nm <- c(
    do.call(sprintf, args = c(list(fmt = 'Early Termination\n%.1f%% vs. %.1f%%'), as.list(1e2*sm@frail))),
    do.call(sprintf, args = c(list(fmt = 'Fail\n%.1f%% vs. %.1f%%'), as.list(1e2*(1 - sm@frail - sm@reject)))), 
    #do.call(sprintf, args = c(list(fmt = 'Success\n\u03b1 = %.1f%%, 1-\u03b2 = %.1f%%'), as.list(1e2*(sm@reject))))
    do.call(sprintf, args = c(list(fmt = 'Success\nalpha = %.1f%%, power = %.1f%%'), as.list(1e2*sm@reject)))
  )
  
  dd <- cbind(sm@frail, 1 - sm@frail - sm@reject, sm@reject)
  list(
    # ?ggplot2::geom_rect wont work here
    geom_bar(mapping = aes(x = 2, y = dd[1L,], fill = nm), alpha = c(.3, .3, 1), stat = 'identity', color = 'white'),
    geom_bar(mapping = aes(x = 1, y = dd[2L,], fill = nm), alpha = c(.3, .3, 1), stat = 'identity', color = 'white'),
    coord_polar(theta = 'y', direction = -1),
    xlim(.3, 2.5),
    labs(fill = sprintf(
      fmt = 'Simon\'s 2-Stage Design\n%s\nResponse Rates: pu=%d%% vs. pa=%d%%\nExpected Total #: %.1f vs. %.1f', 
      switch(object@type, minimax = {
        'Minimum Total Sample Size'
      }, optimal = {
        'Minimum Expected Total Sample Size'
      }, n1 = {
        'Minimum Stage-1 Sample Size'
      }, maximax = {
        'Maximum Total Sample Size'
      }, '(Customized)'),
      1e2*pu, 1e2*pa, sm@eN[1L], sm@eN[2L]))
  )
}





#' @rdname gg_ph2simon
#' @export autolayer.ph2simon
#' @export
autolayer.ph2simon <- function(object, ...) object |> ph2simon4(...) |> autolayer.ph2simon4()

