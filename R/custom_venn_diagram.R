source("R/SourceFiles.R")

custom_venn_diagram <- function(conditionA, conditionB, DataA, DataB){
        venn <- makeVennDiagram(list(conditionA = DataA, conditionB = DataB), #nolint
            fill = c("skyblue", "salmon"), 
            cat.col = c("blue", "red"),
            main = paste("Peak Overlap between", conditionA, "and", conditionB),
            plot = TRUE
        )
        return(venn)
}