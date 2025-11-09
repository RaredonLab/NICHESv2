#' Save metadata to disk, formatted for import to Xenium Explorer
#' 
#' @param object A Seurat object with cell barcodes (does not need to contain all cell barcodes present in Xenium Explorer)
#' @param metadata.to.export Vector of metadata column names to export and view in Xenium Explorer. One file will be saved per element.
#' @param tag String to make the export filenames unique. Optional.
#' @export

ExportForXeniumExplorer <- function(object,
                                    metadata.to.export,
                                    tag = NULL){
  output.list <- list()
  for(i in 1:length(metadata.to.export)){
    output.list[[i]] <- data.frame(cell_id = colnames(object),
                                   group = object[[metadata.to.export[i]]][,metadata.to.export[i]])
  }
  for(i in 1:length(output.list)){
    data.table::fwrite(output.list[[i]],
                       paste(metadata.to.export[i],tag,'exported',Sys.Date(),'csv',sep='.'))
  }
}

# Usage example
# tag = 'test'
# object = PPLR.master
# metadata.to.export <- c('Class.6','SubClass.6','Type.6')
# ExportForXeniumExplorer(object = object,
#                         metadata.to.export = metadata.to.export
#                         tag = 'End')