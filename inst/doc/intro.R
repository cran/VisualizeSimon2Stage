## -----------------------------------------------------------------------------
#| message: false
library(VisualizeSimon2Stage)
library(clinfun)
library(flextable)
library(ggplot2)


## -----------------------------------------------------------------------------
#| echo: false
library(knitr) # for tables in this vignette
#options(mc.cores = 1L) # for CRAN submission


## -----------------------------------------------------------------------------
#| echo: false
#| results: asis
c(
  '[|>](https://search.r-project.org/R/refmans/base/html/pipeOp.html)', 'Forward pipe operator introduced in `R` 4.1.0', 
  '`$`', '[Extract](https://search.r-project.org/R/refmans/base/html/Extract.html) parts of an object',
  '[`binom`](https://search.r-project.org/R/refmans/stats/html/Binomial.html)', '[Binomial](https://en.wikipedia.org/wiki/Binomial_distribution) density and distribution', 
  '`CRAN`, `R`', '[The Comprehensive R Archive Network](https://cran.r-project.org)',
  '[`class`](https://search.r-project.org/R/refmans/base/html/class.html)', 'Object class',
  '[`flextable`](https://search.r-project.org/CRAN/refmans/flextable/html/flextable.html)', 'Flexible tables',
  '`PASS`', 'Power Analysis & Sample Size, <https://www.ncss.com/software/pass/>', 
  '`PET`', 'Probability of early termination',
  '[`ph2simon`](https://search.r-project.org/CRAN/refmans/clinfun/html/ph2simon.html)', 'Simon\'s 2-stage Phase II design',
  '`S3`, `generic`, [`methods`](https://search.r-project.org/R/refmans/utils/html/methods.html)', '`S3` object oriented system, [`UseMethod`](https://search.r-project.org/R/refmans/base/html/UseMethod.html); [`getS3method`](https://search.r-project.org/R/refmans/utils/html/getS3method.html); <https://adv-r.hadley.nz/s3.html>',  
  '`S4`, `generic`, `methods`', '`S4` object oriented system, [`isS4`](https://search.r-project.org/R/refmans/base/html/isS4.html); [`setClass`](https://search.r-project.org/R/refmans/methods/html/setClass.html); [`getMethod`](https://search.r-project.org/R/refmans/methods/html/getMethod.html); <https://adv-r.hadley.nz/s4.html>',
  '[`search`](https://search.r-project.org/R/refmans/base/html/search.html)', 'Search path',
  '[`seed`](https://search.r-project.org/R/refmans/base/html/Random.html)', 'Random number generation seed'
) |>
  matrix(nrow = 2L, dimnames = list(c('Term / Abbreviation', 'Description'), NULL)) |>
  t.default() |>
  as.data.frame.matrix() |> 
  kable() 


## -----------------------------------------------------------------------------
(x = ph2simon(pu = .2, pa = .4, ep1 = .05, ep2 = .1)) 


## -----------------------------------------------------------------------------
x |> ph2simon4(type = 'all')


## -----------------------------------------------------------------------------
x |> simon_pr(prob = c(.2, .3, .4)) |> as_flextable()


## -----------------------------------------------------------------------------
x |> ph2simon4() |> simon_pr(prob = c(.2, .3, .4)) |> as_flextable()


## -----------------------------------------------------------------------------
simon_pr.ph2simon4(prob = c(.2, .3, .4), r1 = 5L, n1 = 24L, r = 13L, n = 45L) |>
  as_flextable()


## -----------------------------------------------------------------------------
#| fig-width: 5
x |> autoplot(type = 'optimal')


## -----------------------------------------------------------------------------
#| fig-width: 5
x |> ph2simon4(type = 'optimal') |> autoplot()


## -----------------------------------------------------------------------------
#| fig-width: 5
autoplot.ph2simon4(pu = .2, pa = .4, r1 = 4L, n1 = 19L, r = 15L, n = 54L, type = 'optimal')


## -----------------------------------------------------------------------------
set.seed(15); s = x |> r_simon(R = 1e4L, prob = .3, type = 'optimal')


## -----------------------------------------------------------------------------
set.seed(15); s1 = x |> ph2simon4(type = 'optimal') |> r_simon(R = 1e4L, prob = .3)
stopifnot(identical(s, s1))


## -----------------------------------------------------------------------------
set.seed(15); s2 = r_simon.ph2simon4(R = 1e4L, prob = .3, r1 = 4L, n1 = 19L, r = 15L, n = 54L)
stopifnot(identical(s, s2))


## -----------------------------------------------------------------------------
set.seed(31); x |> r_simon(R = 1e4L, prob = .2) |> 
  attr(which = 'dx', exact = TRUE) |> 
  as_flextable() |> set_caption(caption = 'pu = .2')


## -----------------------------------------------------------------------------
set.seed(24); x |> r_simon(R = 1e4L, prob = .4) |>
  attr(which = 'dx', exact = TRUE) |>
  as_flextable() |> set_caption(caption = 'pa = .4')


## -----------------------------------------------------------------------------
#| eval: false
#| echo: false
#| results: hide
# summary(x)
# summary(x, type = 'all')


## -----------------------------------------------------------------------------
p = c(A = .3, B = .2, C = .15)


## -----------------------------------------------------------------------------
#| fig-width: 5
set.seed(52); x |> simon_oc(prob = p, R = 1e4L, type = 'optimal')


## -----------------------------------------------------------------------------
#| fig-width: 5
set.seed(52); x |> ph2simon4(type = 'optimal') |> simon_oc(prob = p, R = 1e4L)


## -----------------------------------------------------------------------------
#| fig-width: 5
set.seed(52); simon_oc.ph2simon4(prob = p, R = 1e4L, r1 = 4L, n1 = 19L, r = 15L, n = 54L)

