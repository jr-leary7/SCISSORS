#' A clean `ggplot2` theme for dimension reduction plots.
#'
#' This is the Yeh Lab's default `ggplot2` theme for dimension reduction scatterplots, and was used throughout
#' the SCISSORS manuscript. Like all `ggplot2`-based themes, you can add more themes to it, or override them if you wish.
#' @import ggplot2
#' @export
#' @examples
#' DimPlot(pbmc, reduction = "tsne") + theme_yehlab()

theme_yehlab <- function() {
  theme(legend.position = "bottom",
        axis.text = element_blank(),
        axis.line = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        legend.direction = "horizontal",
        legend.justification = "center",
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", linetype = 1, size = 1))
}
