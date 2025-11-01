AggregateNodeData <- function(node.object,
                              group.by,
                              global.node.list){
  # check
  if(!(group.by %in% colnames(node.object))){
    stop('group.by is not in column names of node.object')
  }else{}
  node.aggregate <- data.frame(prop.table(table(node.object[[group.by]])))
  colnames(node.aggregate) <- c('node.label','system.fraction')
  rownames(node.aggregate) <- node.aggregate$node.label
  
  # add global nodes
  if(length(global.node.list)>length(unique(node.aggregate$node.label))){
    additional.nodes <- data.frame(node.label = global.node.list[which(!(global.node.list%in%node.aggregate$node.label))],
                                   system.fraction = 0)
    rownames(additional.nodes) <- additional.nodes$node.label
    node.aggregate <- rbind(node.aggregate,additional.nodes)
  }
  
  # order nodes based on global.node.list
  node.aggregate <- node.aggregate[global.node.list,]
  
  return(node.aggregate)
}
