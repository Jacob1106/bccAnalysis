
#' @title Filtering Function that filter by a custom p value and chromosome list
#' @param object1 GRanges Object
#' @param object2 GRanges Object 
#' @param overlapobject1 The overlap G ranges object of object1 
#' @param overlapobject2 The overlap G ranges object of object2 
#' @param object1_unique The unique peaks as a G ranges object of object1 
#' @param object2_unique The unique peaks as a G ranges object of object2 
#' @return Returns a tibble that shows: Proportion of Overlap, TotalPeak coverage, Total number of unique bases, number of overlap peaks and percentage overlap of peaks for object 1 and 2.   
#' @export
proportionOverlapTibble <- function(object1, object2, overlapobject1, overlapobject2, object1_unique, object2_unique){
    total_base_overlap <- 0 # nolint
    totalwidthObj1 <- sum(width(object1)) #nolint
    totalwidthObj2 <- sum(width(object2)) #nolint

    for (x in seq_along(overlapobject1)){ #nolint
        over_lap_start <- max(start(overlapobject1[x]), start(overlapobject2[x])) #nolint
        over_lap_end <- min(end(overlapobject1[x]), end(overlapobject2[x])) #nolint
        overlap_length <- max(0, over_lap_end - over_lap_start)
        total_base_overlap <- total_base_overlap + overlap_length
    }

    object1_prop_Overlap <- total_base_overlap / totalwidthObj1 #nolint
    object2_prop_Overlap <- total_base_overlap / totalwidthObj2  #nolint

    total_base_unique_Obj1 <- sum(width(object1_unique)) #nolint
    total_base_unique_Obj2 <- sum(width(object2_unique)) #nolint

    cond_name1 <- deparse(substitute(object1))
    cond_name2 <- deparse(substitute(object2))

    propOverlap_Object1 <- length(overlapobject1) / length(object1)
    percentageOverlap_Obj1 <- propOverlap_Object1 * 100

    propOverlap_Object2 <- length(overlapobject2) / length(object2)
    percentageOverlap_Obj2 <- propOverlap_Object2 * 100

    totalPeaksObj1 <- length(object1) 
    totalPeaksObj2 <- length(object2)

    num_overlap_peaksobj1 <- length(overlapobject1)
    num_overlap_peaksobj2 <- length(overlapobject2)




    resultsTableBaseOverlap <- tibble( #nolint
    sample = c(cond_name1, cond_name2),
    ProportionOverlap = c(object1_prop_Overlap, object2_prop_Overlap),
    TotalPeakCoverage = c(totalwidthObj1, totalwidthObj2),
    UniqueBase = c(total_base_unique_Obj1, total_base_unique_Obj2), #nolint
    TotalPeaks = c(totalPeaksObj1, totalPeaksObj2),
    NumberofOverlapPeaks = c(num_overlap_peaksobj1, num_overlap_peaksobj2),
    PercentageOverlapPeaks = c(percentageOverlap_Obj1, percentageOverlap_Obj2)
    )
    return(resultsTableBaseOverlap)


}

