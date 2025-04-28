
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
#' @param type see slot `@type` of \linkS4class{ph2simon4} object
#' 
#' @returns 
#' Function [autoplot.ph2simon()] returns a \link[ggplot2]{ggplot} object.
#' 
#' Function [autolayer.ph2simon()] returns a \link[base]{list} of \link[ggplot2]{ggproto} and labels.
#' 
#' @keywords internal
#' @name gg_ph2simon
#' @importFrom ggplot2 autoplot
#' @importFrom grid unit
#' @export autoplot.ph2simon
#' @export
autoplot.ph2simon <- function(object, ...) {
  object |> 
    ph2simon4(...) |>
    autoplot.ph2simon4()
}

#' @rdname gg_ph2simon
#' @importFrom ggplot2 autoplot ggplot theme_void
#' @export autoplot.ph2simon4
#' @export
autoplot.ph2simon4 <- function(object, ...) {
  ggplot() + autolayer.ph2simon4(object, ...) + 
    coord_polar(theta = 'x') +
    theme_void()
}


#' @rdname gg_ph2simon
#' @importFrom ggplot2 autolayer aes geom_rect coord_polar ylim
#' @importFrom geomtextpath geom_textpath
#' @importFrom scales pal_hue
#' @export autolayer.ph2simon4
#' @export
autolayer.ph2simon4 <- function(
    object,
    r1 = object@r1, n1 = object@n1, r = object@r, n = object@n,
    pu = object@pu, pa = object@pa,
    type = object@type,
    ...
) {
  
  sm <- simon_pr.ph2simon4(prob = c(pu, pa), n1 = n1, n = n, r1 = r1, r = r)
  nm <- c(
    'Early\nTermination',
    'Fail', 
    'Success'
  )
  
  dd <- cbind(sm@frail, 1 - sm@frail - sm@reject, sm@reject)
    
  cdd <- cbind(0, dd) |>
    apply(MARGIN = 1L, FUN = cumsum, simplify = FALSE)
  xu <- cdd[[1L]]
  xa <- cdd[[2L]]
  
  dd[] <- sprintf(fmt = '%.1f%%', 1e2*dd)
  
  hue <- pal_hue()(n = 3L)
  
  list(
    geom_textpath(mapping = aes(x = c(0, 0, .5), y = c(1.4, .6, 1.4), label = c(
      sprintf(fmt = 'p\u1d64=%.0f%%', 1e2*pu),
      sprintf(fmt = 'p\u2090=%.0f%%', 1e2*pa),
      switch(type, minimax = {
        'Minimum Total Sample Size'
      }, optimal = {
        'Minimum Expected Total Sample Size'
      }, n1 = {
        'Minimum Stage-1 Sample Size'
      }, maximax = {
        'Maximum Total Sample Size'
      })
    ))),
    
    # outer circle: p_u
    geom_rect(mapping = aes(xmin = xu[-4L], xmax = xu[-1L], ymin = 1, ymax = 1.3), fill = hue, alpha = .15, color = 'white'),
    geom_rect(mapping = aes(xmin = xu[3L], xmax = xu[4L], ymin = 1, ymax = 1.3), fill = hue[3L], alpha = .7, color = 'white'),
    geom_textpath(
      mapping = aes(
        x = (xu[1:2] + xu[2:3])/2, 
        y = 1.2, 
        label = paste(c('Early Termination:', 'Fail:'), dd[1L, 1:2])
      ), color = hue[1:2], fontface = 2L
    ),
    
    # inner circle: p_a
    geom_rect(mapping = aes(xmin = xa[-4L], xmax = xa[-1L], ymin = .65, ymax = .95), fill = hue, alpha = .15, color = 'white'),
    geom_rect(mapping = aes(xmin = xa[3L], xmax = xa[4L], ymin = .65, ymax = .95), fill = hue[3L], alpha = .7, color = 'white'),
    geom_textpath(mapping = aes(x = (xa[1:2] + xa[2:3])/2, y = .85, label = dd[2L, 1:2]), color = hue[1:2], fontface = 2L),
    
    #geom_textpath(mapping = aes(x = sum(xa[3:4])/2, y = 1.2, label = dd[1L, 3L]), color = 'white'), # geom_textpath error, must be length-(2+). 
    geom_textpath(
      mapping = aes(
        x = c(sum(xu[3:4])/2, sum(xa[3:4])/2), 
        y = c(1.2, .85), 
        label = c(dd[1L,3L], dd[2L,3L] |> sprintf(fmt = 'Success: %s'))
      ), 
      color = 'white', fontface = 2L
    ),
    
    ylim(c(0, 1.5))
    
  )
  
}





#' @rdname gg_ph2simon
#' @export autolayer.ph2simon
#' @export
autolayer.ph2simon <- function(object, ...) object |> ph2simon4(...) |> autolayer.ph2simon4()

