

#' @title \link[utils]{bibentry} of \CRANpkg{VisualizeSimon2Stage}
#' 
#' @param key,... additional parameters of function \link[utils]{bibentry}
#' 
#' @keywords internal
#' @name simon_bib
#' @importFrom utils bibentry person
#' @export
.simon89 <- \(key = 'Simon89', ...) {
  bibentry(
    bibtype = 'article', key = key, ...,
    title = 'Optimal two-stage designs for phase {II} clinical trials',
    journal = 'Controlled Clinical Trials',
    volume = '10',
    number = '1',
    pages = '1-10',
    year = '1989',
    issn = '0197-2456',
    doi = '10.1016/0197-2456(89)90015-9',
    author = person(given = 'Richard', family = 'Simon')
  )
}