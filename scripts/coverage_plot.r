zip = function (...) {
    as.vector(do.call(rbind, list(...)))
}

#' Convert an Rle object to a ribbon surface
#'
#' @param rle the Rle object
rle_to_surface = function (rle) {
    modules::import_package('S4Vectors', attach = TRUE)
    tibble(Position = zip(start(rle), end(rle)),
           Coverage = rep(runValue(rle), each = 2)) %>%
        filter(! duplicated(Position))
}

case_of = function (x, ...) {
    cases = match.call(expand.dots = FALSE)$...

    for (case in cases) {
        if (case[[1]] == quote(`~`)) {
            if (x == eval.parent(case[[2]])) {
                return(eval.parent(case[[3]]))
            }
        } else {
            return(eval.parent(case))
        }
    }
    stop('No matching case found')
}

#' Scale for genomic positions
#'
#' @param multiplier the numeric multiplier (1, 1e3, 1e6 or 1e9) that
#' determines which unit prefix to use (\code{bp}, \code{kb}, \code{Mb} or
#' \code{Gb}). Defaults to 1e6 (\code{Gb}).
#' @param ... remaining arguments, forwarded to \code{scale_x_continuous}.
scale_x_genomic = function (multiplier = 1e6, ...) {
    unit = case_of(multiplier,
        1 ~ 'bp',
        1e3 ~ 'kb',
        1e6 ~ 'Mb',
        1e9 ~ 'Gb',
        stop('Invalid multiplier (not one of 1, 1e3, 1e6 or 1e9)'))

    labeller = function (x) sprintf('%s %s', x / multiplier, unit)
    ggplot2::scale_x_continuous(labels = labeller, ...)
}

#' Plot a coverage track
#'
#' \code{ggcoverage} plots a \code{GRanges} coverage track.
#'
#' @param data \code{GRanges} object with associated \code{score}s.
#' @return \code{ggcoverage} returns a \code{ggplot} object with associated
#' aesthetics for x and y coordinates.
ggcoverage = function (data) {
    ggplot2::ggplot(data = coverage_layout(data), envir = parent.frame()) +
        ggplot2::aes(Position, Coverage)
}

#' \code{coverage_layout} transforms a coverage track into a \code{data.frame}
#' suitable for plotting with ggplot2.
#'
#' @return \code{coverage_layout} returns a \code{data.frame} with columns
#' \code{Position} and \code{Coverage}, suitable for plotting via ggplot2.
#' @rdname ggcoverage
coverage_layout = function (data) {
    lapply(gr$coverage(data, weight = 'score'), rle_to_surface) %>%
        Map(function (name, x) mutate(x, seqname = name), names(.), .) %>%
        bind_rows()
}
