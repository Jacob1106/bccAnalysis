source("R/SourceFiles.R")

import_samples <- function(named_file_list){
    GRANGE_List <- list()
    for (x in seq_along(named_file_list)){
        import_data <- rtracklayer::import(named_file_list[[x]], format = "narrowPeak")
        GRANGE_List <- append(GRANGE_List, setNames(list(import_data), names(named_file_list)[x]))
    }
    return(GRANGE_List)

}