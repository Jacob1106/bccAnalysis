
#' @title Width Comparison of 2 GRanges Objects
#' @param Obj1 GRanges Object 
#' @param Obj2 GRanges Object
#' @return GGplot | Returns a plot that shows the width of Object 1 and Object 2. Plot = (x = Width of Peak) (y = Density of a particular width (How many peaks have that width))
#' @importFrom GenomicRanges width
#' @importFrom ggplot2 ggplot aes geom_density scale_x_log10 theme_classic labs
#' @author Jacob Martin
#' @export
compareGRangeWidth <- function(Obj1, Obj2){
    widthAnalysis <- data.frame(
        width = c(width(Obj1), width(Obj2)), #nolint
        condition = rep(c(deparse(substitute(Obj1)), deparse(substitute(Obj2))), c(length(Obj1), length(Obj2))) #nolint
    )
    widthPlot <- ggplot(widthAnalysis, aes(x = width, fill = condition))+ #nolint
        geom_density(alpha = 0.4)+ #nolint 
        scale_x_log10()+ #nolint
        theme_classic()+ # nolint: object_usage_linter.
        labs(title = "Peak width comparison", x = "Peak width (bp)", y = "Density") #nolint
    return(widthPlot)
}
