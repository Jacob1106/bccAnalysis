source("R/SourceFiles.R")

overlapGRangeList <- function(grange_list, num_of_overlaps_required = 1){
        unlist <- unlist(grange_list)
        index <- which(countOverlaps(unlist) > num_of_overlaps_required) #nolint
        filtered_GRange <- unlist[index]
        returnGRange <- GenomicRanges::reduce(filtered_GRange)
        return(returnGRange)
}