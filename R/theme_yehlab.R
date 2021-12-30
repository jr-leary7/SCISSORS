#' A clean \code{ggplot2} theme for dimension reduction plots.
#'
#' @name theme_yehlab
#' @author Jack Leary
#' @description This is the Yeh Lab's default \code{ggplot2} theme for dimension reduction scatterplots, and was used throughout the SCISSORS manuscript.
#' Like all \code{ggplot2}-based themes, you can add more themes to it, or override them if you wish.
#' @importFrom ggplot2 theme element_blank element_rect
#' @export
#' @examples
#' \dontrun{DimPlot(pbmc, reduction = "umap") + theme_yehlab()}

theme_yehlab <- function() {
  ggplot2::theme(legend.position = "bottom",
                 axis.text = ggplot2::element_blank(),
                 axis.line = ggplot2::element_blank(),
                 panel.grid = ggplot2::element_blank(),
                 axis.ticks = ggplot2::element_blank(),
                 legend.direction = "horizontal",
                 legend.justification = "center",
                 panel.background = ggplot2::element_blank(),
                 panel.border = ggplot2::element_rect(colour = "black", linetype = 1, size = 1))
}
