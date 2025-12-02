
#' @title Comparing Overlap of two GRange Objects
#' @param conditionA String with the name of the DataA
#' @param conditionB String with the name of DataB
#' @param DataA GRanges Object 
#' @param DataB GRanges Object
#' @return GGplot | Venn Diagram showing overlaps between two GRanges Objects. 
#' @export
custom_venn_diagram <- function(conditionA, conditionB, DataA, DataB){
        venn <- makeVennDiagram(list(conditionA = DataA, conditionB = DataB), #nolint
            fill = c("skyblue", "salmon"), 
            cat.col = c("blue", "red"),
            main = paste("Peak Overlap between", conditionA, "and", conditionB),
            plot = TRUE
        )
        return(venn)
}