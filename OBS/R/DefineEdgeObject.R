DefineEdgeObject <- function(connect.obj,
                             feature,
                             assay = 'CellToCell'){
  # check & initialize
  if(!(feature %in% rownames(connect.obj))){
    stop('feature not present in connect.obj')
  }else{}

  # extract & store information
  rownames(connect.obj[[assay]]@layers$data) <- rownames(connect.obj)
  colnames(connect.obj[[assay]]@layers$data) <- colnames(connect.obj)
  edge.object <- data.frame(feature.value = connect.obj@assays[[assay]]@layers$data[feature,])
  edge.object$barcode.pair <- rownames(edge.object)

  # extract metadata
  edge.object <- cbind(edge.object,connect.obj@meta.data)
  return(edge.object)
}
