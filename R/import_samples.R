source("R/SourceFiles.R")
#' @title Importing narrowPeak files and converting into a GRanges Object 
#' @param named_file_list A named list with file paths for narrowPeak files. These files will be imported as GRanges Objects. 
#' @return Returns a list of GRange Objects with a new column for normal p values. The list is sorted and edited based on the p value and chromosome inputs.  
#' @export
import_samples <- function(named_file_list){
    GRANGE_List <- list()
    for (x in seq_along(named_file_list)){
        import_data <- rtracklayer::import(named_file_list[[x]], format = "narrowPeak")
        GRANGE_List <- append(GRANGE_List, setNames(list(import_data), names(named_file_list)[x]))
    }
    return(GRANGE_List)

}