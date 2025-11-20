source("R/SourceFiles.R")


filter_by_pvalue <- function(list, pvalue, chromosomes = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")){
    list <- list
    for (x in seq_along(list)){
        #Convert log 10 p values to normal P values
        individual <- list[[x]]
        individual$pValue_norm <- 10^(-individual$pValue)

        #Filtering P value: 
        individual <- individual[which(individual$pValue_norm < pvalue), ]

        individual <- individual[seqnames(individual) %in% chromosomes] # nolint

        list[[x]] <- individual 


    }
    return(list)
}