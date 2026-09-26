## -----------------------------------------------------------------------------
library(VisualizeSimon2Stage)
library(clinfun)
library(flextable)
library(ggplot2)


## -----------------------------------------------------------------------------
#| echo: false
theme_minimal() |>
  theme_set()


## -----------------------------------------------------------------------------
(x = ph2simon(pu = .2, pa = .4, ep1 = .05, ep2 = .1)) 


## -----------------------------------------------------------------------------
x |> 
  ph2simon4(type = 'all')


## -----------------------------------------------------------------------------
x |> 
  simon_pr(prob = c(.2, .3, .4)) |> 
  as_flextable()


## -----------------------------------------------------------------------------
x |> 
  ph2simon4() |> 
  simon_pr(prob = c(.2, .3, .4)) |> 
  as_flextable()


## -----------------------------------------------------------------------------
simon_pr.ph2simon4(prob = c(.2, .3, .4), r1 = 5L, n1 = 24L, r = 13L, n = 45L) |>
  as_flextable()


## -----------------------------------------------------------------------------
#| fig-width: 4
#| fig-height: 4
x |> 
  autoplot(type = 'optimal')


## -----------------------------------------------------------------------------
#| fig-width: 4
#| fig-height: 4
x |> 
  ph2simon4(type = 'optimal') |> 
  autoplot()


## -----------------------------------------------------------------------------
#| fig-width: 4
#| fig-height: 4
autoplot.ph2simon4(pu = .2, pa = .4, r1 = 4L, n1 = 19L, r = 15L, n = 54L, type = 'optimal')


## -----------------------------------------------------------------------------
set.seed(15); s = x |> 
  r_simon(R = 1e4L, prob = .3, type = 'optimal')


## -----------------------------------------------------------------------------
set.seed(15); s1 = x |> 
  ph2simon4(type = 'optimal') |> 
  r_simon(R = 1e4L, prob = .3)
stopifnot(identical(s, s1))


## -----------------------------------------------------------------------------
set.seed(15); s2 = r_simon.ph2simon4(R = 1e4L, prob = .3, r1 = 4L, n1 = 19L, r = 15L, n = 54L)
stopifnot(identical(s, s2))


## -----------------------------------------------------------------------------
#| code-fold: true
#| code-summary: 'R code: Type-I-error rate at $p_u$, `1e4L` simulated copies'
set.seed(31); x |> 
  r_simon(R = 1e4L, prob = .2) |> 
  attr(which = 'dx', exact = TRUE) |> 
  table(Decision = _) |>
  as_flextable() |> 
  highlight(i = 3L, j = 3L) |>
  add_header_lines(values = 'pu = .2')


## -----------------------------------------------------------------------------
#| code-fold: true
#| code-summary: 'R code: power at $p_a$, `1e4L` simulated copies'
set.seed(24); x |> 
  r_simon(R = 1e4L, prob = .4) |>
  attr(which = 'dx', exact = TRUE) |>
  table(Decision = _) |>
  as_flextable() |> 
  highlight(i = 3L, j = 3L) |>
  add_header_lines(values = 'pa = .4')


## -----------------------------------------------------------------------------
#| eval: false
#| echo: false
#| results: hide
# summary(x)
# summary(x, type = 'all')


## -----------------------------------------------------------------------------
#| fig-width: 4
#| fig-height: 3
x |> 
  powerCurve.ph2simon()


## -----------------------------------------------------------------------------
#| fig-width: 4
#| fig-height: 3
x |> 
  ph2simon4() |>
  powerCurve.ph2simon4()


## -----------------------------------------------------------------------------
p = c(A = .3, B = .2, C = .15)


## -----------------------------------------------------------------------------
#| fig-width: 4
#| fig-height: 4
set.seed(52); x |> 
  simon_oc(prob = p, R = 1e4L, type = 'optimal')


## -----------------------------------------------------------------------------
#| fig-width: 4
#| fig-height: 4
set.seed(52); x |> 
  ph2simon4(type = 'optimal') |> 
  simon_oc(prob = p, R = 1e4L)


## -----------------------------------------------------------------------------
#| fig-width: 4
#| fig-height: 4
set.seed(52); simon_oc.ph2simon4(prob = p, R = 1e4L, r1 = 4L, n1 = 19L, r = 15L, n = 54L)

