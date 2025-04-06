

#' @rdname simon_pr
#' 
#' @param prob \link[base]{double} scalar or \link[base]{vector}, 
#' true response rate(s) \eqn{p}
#' 
#' @param object a \link[clinfun]{ph2simon} or \linkS4class{ph2simon4} object
#' 
#' @param r1,n1,r,n (optional) \link[base]{integer} scalars, see \linkS4class{ph2simon4}.
#' 
#' @param ... parameters of function [ph2simon4()], most importantly `type`
#' 
#' @returns 
#' 
#' Function [simon_pr()] returns \linkS4class{simon_pr} object.
#' 
#' @keywords internal
#' @export
simon_pr <- function(prob, object, ...) UseMethod(generic = 'simon_pr', object = object)


#' @rdname simon_pr
#' @export simon_pr.ph2simon
#' @export
simon_pr.ph2simon <- function(prob, object, ...) {
  object |>
    ph2simon4(...) |>
    simon_pr.ph2simon4(prob = prob)
}


#' @rdname simon_pr
#' @importFrom stats dbinom pbinom
#' @export simon_pr.ph2simon4
#' @export
simon_pr.ph2simon4 <- function(
    prob, 
    object,
    r1 = object@r1, n1 = object@n1, r = object@r, n = object@n,
    ...
) {
  
  if (length(r1) != 1L || !is.integer(r1) || is.na(r1) || r1 < 0L) stop('`r1` must be non-negative integer scalar')
  if (length(n1) != 1L || !is.integer(n1) || is.na(n1) || n1 <= 0L) stop('`n1` must be positive integer scalar')
  if (length(r) != 1L || !is.integer(r) || is.na(r) || r < 0L) stop('`r` must be non-negative integer scalar')
  if (length(n) != 1L || !is.integer(n) || is.na(n) || n <= 0L) stop('`n` must be positive integer scalar')
  
  if (!length(prob) || !is.double(prob) || anyNA(prob) || any(prob < 0, prob > 1)) stop('`prob` must be (0,1) vector')
  
  s1 <- (r1 + 1L):n1 # responses needed in Stage-1 to continue
  
  frail <- pbinom(q = r1, size = n1, prob = prob, lower.tail = TRUE) 
  # early termination; Prob(X1 <= r1), where X1 ~ Binom(n1, p)
  # vectorize on `prob`
  
  reject <- vapply(prob, FUN = \(p) {
    # `p`: response rate
    d1s <- dbinom(x = s1, size = n1, prob = p)
    sum(d1s * pbinom(q = r - s1, size = n - n1, prob = p, lower.tail = FALSE))
    # success in rejection H0; \sum_{s1 = r1+1}^{n1} Prob(X1 = s1) * Prob(X2 > (r-s1))
    # .. where X1 ~ Binom(n1, p), X2 ~ Binom(n-n1, p); X1 and X2 are independent
    # ?stats::pbinom is okay with negative `q`
  }, FUN.VALUE = NA_real_) # *cannot* vectorize on `prob`!!!
  
  new(Class = 'simon_pr',
      frail = frail,
      reject = reject,
      eN = frail * n1 + (1 - frail) * n, 
      prob = prob)
  
}






#' @export
`[.simon_pr` <- function(x, i) {
  x@frail <- x@frail[i]
  x@reject <- x@reject[i]
  x@eN <- x@eN[i]
  x@prob <- x@prob[i]
  return(x)
}



#' @title Convert \linkS4class{ph2simon4} to \link[flextable]{flextable}
#' 
#' @param x a \linkS4class{ph2simon4}
#' 
#' @param ... additional parameters, currently not in use
#' 
#' @returns 
#' Function [as_flextable.simon_pr()] returns a \link[flextable]{flextable}
#' 
#' @keywords internal
#' @importFrom officer fp_border
#' @importFrom flextable as_flextable flextable autofit add_header_row align vline
#' @export as_flextable.simon_pr
#' @export
as_flextable.simon_pr <- function(x, ...) {
  p <- cbind(
    'Early Termination' = x@frail, 
    'Fail' = 1-x@frail-x@reject, 
    'Success' = x@reject)
  p[] <- sprintf(fmt = '%.1f%%', 1e2*p)
  
  border_hard_ <- fp_border(width = 1.5, color = 'gray40')
  # *looks* like default border used in ?flextable::flextable
  # tzh does *not* know how to find out for sure, for now..
  # ?flextable:::print.flextable
  # ?flextable::htmltools_value
  
  data.frame(
    'Response Rate' = sprintf(fmt = 'p = %.0f%%', 1e2*x@prob),
    'E(N)' = sprintf(fmt = '%.1f', x@eN),
    p,
    check.names = FALSE
  ) |>
    flextable() |>
    autofit() |>
    add_header_row(values = c(' ', ' ', 'Probabilities'), colwidths = c(1L, 1L, 3L), top = TRUE) |>
    align(align = 'center', part = 'all') |>
    vline(j = 2L, border = border_hard_, part = 'all')
  
}




#' @importFrom ggplot2 autolayer aes geom_rect coord_polar labs
#' @export
autolayer.simon_pr <- function(object, ...) {
  object <- object[1L] # `[.simon_pr`
  x <- c(object@frail, 1-object@frail-object@reject, object@reject)
  nm <- paste(c('Early Termination', 'Fail', 'Success'), sprintf(fmt = '%.1f%%', 1e2*x))
  ymax <- cumsum(x)
  ymin <- c(0, ymax[-3L])
  list(
    geom_rect(mapping = aes(ymax = ymax, ymin = ymin, xmax = 1, xmin = 0, fill = nm), colour = 'white'),
    coord_polar(theta = 'y'),
    labs(fill = sprintf(fmt = 'Simon\'s 2-Stage\nResponse Rate %.0f%%\nExpected Total # = %.1f', 1e2*object@prob, object@eN))
  )
}


#' @importFrom ggplot2 autoplot ggplot theme_void
#' @export
autoplot.simon_pr <- function(object, ...) {
  ggplot() + autolayer.simon_pr(object, ...) + theme_void()
}


