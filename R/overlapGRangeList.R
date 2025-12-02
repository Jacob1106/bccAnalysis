#' @title Filtering GRange Objects by number of overlaps with other samples. This function is used to merge multiple GRange Objects into one. 
#' @param grange_list List of GRanges Objects that are similar or part of the same condition 
#' @param num_of_overlaps_required Integer value that determines how many overlaps a peak must have with the other GRanges objects in order to be outputed by the function
#' @return Returns a single GRange Object that is merged and filtered by overlaps. 
#' @export
overlapGRangeList <- function(grange_list, num_of_overlaps_required = 1){
        unlist <- unlist(grange_list)
        index <- which(countOverlaps(unlist) > num_of_overlaps_required) #nolint
        filtered_GRange <- unlist[index]
        returnGRange <- GenomicRanges::reduce(filtered_GRange)
        return(returnGRange)
}