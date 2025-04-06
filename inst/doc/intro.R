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
  '', 'Forward pipe operator', '`?base::pipeOp` introduced in `R` 4.1.0', 
  '`$`', 'Extract parts of an object', '`` ?base::`$` ``',
  '`binom`', 'Binomial density and distribution', '`?stats::dbinom`, <https://en.wikipedia.org/wiki/Binomial_distribution>',
  '`CRAN`, `R`', 'The Comprehensive R Archive Network', '<https://cran.r-project.org>',
  '`class`', 'Object class', '`?base::class`',
  '`flextable`', 'Flexible tables', '`?flextable::flextable`',
  '`PASS`', 'Power Analysis & Sample Size', '<https://www.ncss.com/software/pass/>', 
  '`PET`', 'Probability of early termination', '',
  '`ph2simon`', 'Simon\'s 2-stage Phase II design', '`clinfun::ph2simon`',
  '`S3`, `generic`, `methods`', '`S3` object oriented system',  '`base::UseMethod`; `utils::methods`; `utils::getS3method`; <https://adv-r.hadley.nz/s3.html>',
  '`S4`, `generic`, `methods`', '`S4` object oriented system',  '`base::isS4`; `methods::setClass`; `methods::getMethod`; <https://adv-r.hadley.nz/s4.html>',
  '`search`', 'Search path', '`?base::search`',
  '`seed`', 'Random number generation seed', '`?base::set.seed`',
  '`table`', 'Cross tabulation', '`?base::table`'
) |>
  matrix(nrow = 3L, dimnames = list(c('Term / Abbreviation', 'Description', 'Reference'), NULL)) |>
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
x |> autoplot(type = 'optimal')


## -----------------------------------------------------------------------------
x |> ph2simon4(type = 'optimal') |> autoplot()


## -----------------------------------------------------------------------------
set.seed(15); r_u = x |> r_simon(R = 1e4L, prob = .2, type = 'optimal')


## -----------------------------------------------------------------------------
set.seed(15); r_u1 = x |> ph2simon4(type = 'optimal') |> r_simon(R = 1e4L, prob = .2)
stopifnot(identical(r_u, r_u1))


## -----------------------------------------------------------------------------
set.seed(15); r_u2 = r_simon.ph2simon4(R = 1e4L, prob = .2, r1 = 4L, n1 = 19L, r = 15L, n = 54L)
stopifnot(identical(r_u, r_u2))


## -----------------------------------------------------------------------------
r_u |> attr(which = 'dx', exact = TRUE) |> table() |> 
  as_flextable() |> set_caption(caption = 'pu = .2')


## -----------------------------------------------------------------------------
set.seed(24); x |> r_simon(R = 1e4L, prob = .4) |>
  attr(which = 'dx', exact = TRUE) |> table() |> 
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
set.seed(52); x |> simon_oc(prob = p, R = 1e4L, type = 'optimal')


## -----------------------------------------------------------------------------
set.seed(52); x |> ph2simon4(type = 'optimal') |> simon_oc(prob = p, R = 1e4L)


## -----------------------------------------------------------------------------
set.seed(52); simon_oc.ph2simon4(prob = p, R = 1e4L, r1 = 4L, n1 = 19L, r = 15L, n = 54L)

