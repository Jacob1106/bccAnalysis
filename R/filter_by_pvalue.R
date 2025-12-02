#' @title Filtering Function that filter by a custom p value and chromosome list
#' @param list List of GRanges Objects
#' @param pvalue float value p value
#' @param chromosomes list of characters e.g. c("chr1", "chr2")
#' @return Returns a list of GRange Objects with a new column for normal p values. The list is sorted and edited based on the p value and chromosome inputs.  
#' @importFrom GenomicRanges seqnames
#' @author Jacob Martin
#' @export

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