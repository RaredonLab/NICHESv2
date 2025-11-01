#' AggregateEdgeData
#'
#' @param edge.object A connectomic object
#' @param group.by The grouping variable to use to aggregate columns
#' @return A dataframe with aggregated edge data

AggregateEdgeData <- function(edge.object,
                              group.by){
  # check
  if(!(group.by %in% colnames(edge.object))){
    stop('group.by is not in column names of edge.object')
  }else{}
  edge.aggregate <- aggregate(x = edge.object,
                              by = list(edge.label = edge.object[[group.by]]),
                              FUN = mean,drop = F)
  # CHECK What is up with this separator being an EN-dash and not an EM-dash?? Why??
  edge.aggregate$sending.label <- stringr::str_split_fixed( edge.aggregate$edge.label , " - ",n=2)[,1]
  edge.aggregate$receiving.label <- stringr::str_split_fixed( edge.aggregate$edge.label , " - ",n=2)[,2]
  # remove columns that are all NA (character aggregation)
  edge.aggregate <- edge.aggregate[,-which(colSums(is.na(edge.aggregate))>0)]
  return(edge.aggregate)
}

