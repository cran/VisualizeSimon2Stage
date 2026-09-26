

#' @title Power of \link[clinfun]{ph2simon} or \linkS4class{ph2simon4} Object
#' 
#' @param prob \link[base]{double} scalar, true response rate \eqn{p}
#' 
#' @param object a \link[clinfun]{ph2simon} or \linkS4class{ph2simon4} object
#' 
#' @param ... parameters of function [ph2simon4()], most importantly `type`
#' 
#' @param r1,n1,r,n (optional) \link[base]{integer} scalars, see \linkS4class{ph2simon4}.
#' 
#' @name power_ph2simon
#' @export
power.ph2simon4 <- function(
    prob, 
    object, ...,
    r1 = object@r1, n1 = object@n1, r = object@r, n = object@n
) {
  
  if (length(r1) != 1L || !is.integer(r1) || is.na(r1) || r1 < 0L) stop('`r1` must be non-negative integer scalar')
  if (length(n1) != 1L || !is.integer(n1) || is.na(n1) || n1 <= 0L) stop('`n1` must be positive integer scalar')
  if (length(r) != 1L || !is.integer(r) || is.na(r) || r < 0L) stop('`r` must be non-negative integer scalar')
  if (length(n) != 1L || !is.integer(n) || is.na(n) || n <= 0L) stop('`n` must be positive integer scalar')
  
  pr <- simon_pr.ph2simon4(prob = prob, r1 = r1, n1 = n1, r = r, n = n)
  return(pr@reject)
  
}



#' @rdname power_ph2simon
#' @export
power.ph2simon <- function(prob, object, ...) {
  object |> 
    ph2simon4(...) |> 
    power.ph2simon4(prob = prob, object = _)
}


#' @rdname power_ph2simon
#' @importFrom ggplot2 ggplot aes stat_function geom_point scale_x_continuous scale_y_continuous
#' @importFrom ggrepel geom_label_repel
#' @importFrom scales pal_hue label_percent
#' @export
powerCurve.ph2simon4 <- function(object, ...) {
  
  p <- c(object@pu, object@pa)
  pwr <- power.ph2simon4(prob = p, object = object)
  col <- pal_hue()(2L)
  
  ggplot() + 
    stat_function(
      fun = power.ph2simon4, 
      args = list(object = object), 
      xlim = c(0, 1),
      n = 501L) +
    geom_point(mapping = aes(x = p, y = pwr), colour = col) +
    geom_label_repel(mapping = aes(
      x = p, y = pwr,
      label = pwr |> label_percent(accuracy = .1)()
    ), colour = col) +
    scale_y_continuous(name = 'Power', labels = label_percent()) +
    scale_x_continuous(name = 'True Response Rate', labels = label_percent())
  
}


#' @rdname power_ph2simon
#' @export
powerCurve.ph2simon <- function(object, ...) {
  object |> 
    ph2simon4(...) |> 
    powerCurve.ph2simon4()
}

