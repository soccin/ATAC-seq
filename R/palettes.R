#' Color palettes and scales for the R scripts.
#'
#' Every color a figure uses is defined here, once, as hex values, so a
#' figure does not change when a palette package is upgraded. The rules:
#'
#' - Sequential (positive values): from near white to dark, clipped to
#'   the 5th and 95th percentiles. PAL_EXPR for biology, PAL_QC for
#'   values read for their upper tail.
#' - Diverging (signed values): white at 0, blue negative, red positive,
#'   with symmetric limits +/- the 95th percentile of |value|; values
#'   outside are drawn in the end colors.
#' - Categorical: PAL_CAT up to 23 levels, PAL_CAT32 for 24-32, in a
#'   fixed order. Never viridis and its relatives, never ggplot's
#'   default hue palette.
#' - Grey means "none" only: NA, not significant, unassigned.
#'
#' So far only R/pathwayGSEA.R uses this file.

## Sequential, positive biology: grey -> yellow -> orange -> brown-red
PAL_EXPR <- c("#EDEDED", "#FEE391", "#FEC44F", "#FE9929",
              "#EC7014", "#CC4C02", "#8C2D04")

## Sequential, diagnostics: grey -> gold -> red -> deep purple
PAL_QC <- c("#F2F2F2", "#FFE0A3", "#FDA43A", "#F03B20",
            "#BD0026", "#7A0177", "#49006A")

## Diverging: ColorBrewer RdBu reversed (blue = negative)
PAL_SIGNED <- c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0",
                "#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B",
                "#67001F")

## Categorical, up to 23 levels: pals::cols25() without its two greys
PAL_CAT <- c("#1F78C8", "#FF0000", "#33A02C", "#6A33C2", "#FF7F00",
             "#FFD700", "#A6CEE3", "#FB6496", "#B2DF8A", "#CAB2D6",
             "#FDBF6F", "#EEE685", "#C8308C", "#FF83FA", "#C814FA",
             "#0000FF", "#36648B", "#00E2E5", "#00FF00", "#778B00",
             "#BEBE00", "#8B3B00", "#A52A3C")

## Categorical, 24-32 levels: pals::glasbey(32). Replaces PAL_CAT
## entirely; do not append it to PAL_CAT.
PAL_CAT32 <- c("#0000FF", "#FF0000", "#00FF00", "#000033", "#FF00B6",
               "#005300", "#FFD300", "#009FFF", "#9A4D42", "#00FFBE",
               "#783FC1", "#1F9698", "#FFACFD", "#B1CC71", "#F1085C",
               "#FE8F42", "#DD00FF", "#201A01", "#720055", "#766C95",
               "#02AD24", "#C8FF00", "#886C00", "#FFB79F", "#858567",
               "#A10300", "#14F9FF", "#00479E", "#DC5E93", "#93D4FF",
               "#004CFF", "#F2F318")

## "None" / unassigned
NA_GREY <- "grey70"

## The two sides of a signed value drawn as a category (bars, curves):
## the dark blue and red of PAL_SIGNED, not its darkest ends
SIDE_COLORS <- c(down = PAL_SIGNED[2], up = PAL_SIGNED[10])

## Lines that are not data: reference lines at 0 or a cutoff, outlines
## of points, borders between tiles
REF_LINE <- "grey40"
POINT_OUTLINE <- "grey20"
TILE_BORDER <- "white"

#' Limits of a sequential scale: the 5th and 95th percentiles
#'
#' @param x Numeric values the scale draws
#' @return Numeric vector of 2; the range of x if the percentiles are
#'   equal, and c(0, x) for a single value
seq_limits <- function(x) {
  x <- x[is.finite(x)]
  lims <- unname(quantile(x, c(0.05, 0.95)))
  if (lims[1] == lims[2]) {
    lims <- range(x)
  }
  if (lims[1] == lims[2]) {
    lims <- c(0, lims[2])
  }
  lims
}

#' Sequential color scale, clipped to the 5th and 95th percentiles
#'
#' @param x Numeric values the scale draws
#' @param palette Character hex vector, PAL_EXPR or PAL_QC
#' @param aesthetics Character "fill" or "colour"
#' @param ... Passed to ggplot2::scale_fill_gradientn (name, ...)
#' @return ggplot2 scale
scale_sequential <- function(x, palette, aesthetics = "fill", ...) {
  ggplot2::scale_fill_gradientn(colours = palette, limits = seq_limits(x),
                                oob = scales::squish, na.value = NA_GREY,
                                aesthetics = aesthetics, ...)
}

#' Diverging color scale, symmetric at +/- the 95th percentile of |x|
#'
#' @param x Numeric values the scale draws
#' @param aesthetics Character "fill" or "colour"
#' @param ... Passed to ggplot2::scale_fill_gradientn (name, ...)
#' @return ggplot2 scale
scale_signed <- function(x, aesthetics = "fill", ...) {
  lim <- unname(quantile(abs(x[is.finite(x)]), 0.95))
  if (lim == 0) {
    lim <- 1
  }
  ggplot2::scale_fill_gradientn(colours = PAL_SIGNED,
                                limits = c(-lim, lim),
                                oob = scales::squish, na.value = NA_GREY,
                                aesthetics = aesthetics, ...)
}
