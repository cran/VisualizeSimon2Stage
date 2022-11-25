

# very difficult to overwrite \code{clinfun:::print.ph2simon} when publish to CRAN ..
# print.ph2simon <- function(x, ...) {
#  out <- summary.ph2simon(x)
#  out$EN[] <- sprintf('%.1f', out$EN)
#  out$p[] <- sprintf('%.1f%%', 1e2*out$p) 
#  cat('\n Simon\'s 2-Stage Design\n\n')
#  cat(sprintf('Unacceptable Response Rate: %.1f%%\n', 1e2*x[['pu']]))
#  cat(sprintf('Desirable Response Rate: %.1f%%\n', 1e2*x[['pa']]))
#  cat(sprintf('Controlled Error Rates: \u03B1 \u2264 %.f%%, \u03B2 \u2264 %.f%%\n', 1e2*x[['alpha']], 1e2*x[['beta']]))
#  cat(sprintf('Maximum Sample Size Allowed: %d\n\n', x$nmax))
#  print.default(cbind(out$design, out$EN[, 1L, drop = FALSE]), right = TRUE, quote = FALSE, ...)
#  cat('\n')
#}





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
  
  ret_EN <- do.call(rbind, args = lapply(sm, FUN = function(i) {
    y <- i@eN
    names(y) <- paste0('EN(', c('pu', 'pa'), ')')
    return(y)
  }))
  ret_p <- do.call(rbind, args = lapply(sm, FUN = function(i) {
    i <- unclass(i)
    c('PET(pu)' = i[1L,1L], 'PET(pa)' = i[2L,1L], 'alpha' = i[1L,3L], 'beta' = 1-i[2L,3L])
  }))
  
  # In \link[clinfun]{ph2simon} output, 
  # \code{PET} means probability of early termination (i.e. frail),
  # \code{EN} means expected sample size.
  rownames(ret_design) <- rownames(ret_EN) <- rownames(ret_p) <- names(rid)
  return(list(
    design = ret_design,
    EN = ret_EN,
    p = ret_p
  ))
}








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
          sprintf('Success\n\u03B1 = %.1f%%, 1-\u03B2 = %.1f%%\n', 1e2*dd[1L,3L], 1e2*dd[2L,3L]))
  
  list(
    # ?ggplot2::geom_rect wont work here
    geom_bar(mapping = aes(x = 2, y = dd[1L,], fill = nm), alpha = c(.3, .3, 1), stat = 'identity', color = 'white'),
    geom_bar(mapping = aes(x = 1, y = dd[2L,], fill = nm), alpha = c(.3, .3, 1), stat = 'identity', color = 'white'),
    coord_polar(theta = 'y', direction = -1),
    xlim(.3, 2.5),
    scale_fill_discrete(name = sprintf(
      fmt = 'Simon\'s %s 2-Stage\nResponse Rates: %d%% vs. %d%%\nExpected Total #: %.1f vs. %.1f', 
      type, 1e2*pu, 1e2*pa, sm@eN[1L], sm@eN[2L]))
  )
}


#' @title Plot A Simon's Two-Stage Design
#' 
#' @description Plot \link[clinfun]{ph2simon} object.
#' 
#' @param object \link[clinfun]{ph2simon} object
#' 
#' @param ... potential parameters, currently not in use
#' 
#' @return 
#' \link{autoplot.ph2simon} returns a \link[ggplot2]{ggplot} object
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
#' # example from user feedback
#' (des = ph2simon(pu = .05, pa = .3, ep1 = .05, ep2 = .2))
#' autoplot(des, type = 'optimal')
#' autoplot(des, type = 'minimax')
#' 
#' @export
autoplot.ph2simon <- function(object, ...) {
  ggplot() + autolayer.ph2simon(object, ...) + theme_void()
}



