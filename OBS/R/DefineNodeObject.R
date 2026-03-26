DefineNodeObject <- function(transcr.obj,
                             feature,
                             assay = 'RNA'){
  # check & initialize
  if(!(feature %in% rownames(transcr.obj))){
    stop('feature not present in transcr.obj')
  }else{}

  # extract & store information
  rownames(transcr.obj[[assay]]@layers$data) <- rownames(transcr.obj)
  colnames(transcr.obj[[assay]]@layers$data) <- colnames(transcr.obj)
  node.object <- data.frame(feature.value = transcr.obj[[assay]]@layers$data[feature,])
  node.object$barcode <- rownames(node.object)

  # extract metadata
  node.object <- cbind(node.object,transcr.obj@meta.data)
  return(node.object)
}
