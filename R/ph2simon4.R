
#' @rdname ph2simon4
#' 
#' @param object a \link[clinfun]{ph2simon} object
#' 
#' @param type \link[base]{character} scalar or \link[base]{vector}, see **Slots**
#' 
#' @param ... additional parameters, currently not in use
#' 
#' @examples
#' (x = clinfun::ph2simon(pu = .2, pa = .4, ep1 = .05, ep2 = .1)) 
#' ph2simon4(x, type = 'opt')
#' ph2simon4(x) # default type = 'minimax'
#' @export
ph2simon4 <- function(object, type = 'minimax', ...) {
  
  # [print_ph2simon()] uses multiple `type`
  
  type <- type |>
    vapply(FUN = match.arg, choices = c('minimax', 'optimal', 'n1', 'maximax'), FUN.VALUE = '') |> 
    unname()
  
  x <- object$out # sorted by column 'n'
  # colnames(x)
  # `PET` indicates probability of early termination (i.e. frail),
  # `EN` indicates expected sample size.
  
  rid <- c(
    minimax = 1L, # 'minimum total sample size'
    optimal = which.min(x[,'EN(p0)']), # 'minimum expected total sample size'
    n1 = which.min(x[,'n1']), # 'minimum Stage 1 sample size'
    maximax = dim(x)[1L] # 'maximum total sample size'
  )
  
  x_ <- x[rid[type], c('r1', 'n1', 'r', 'n'), drop = FALSE]
  storage.mode(x_) <- 'integer'
  
  new(Class = 'ph2simon4',
      r1 = x_[,'r1'],
      n1 = x_[,'n1'],
      r = x_[,'r'],
      n = x_[,'n'],
      pu = object$pu, pa = object$pa, alpha = object$alpha, beta = object$beta,
      nmax = as.integer(object$nmax),
      type = type)
  
}
