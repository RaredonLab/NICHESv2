#' CircuitPlot
#'
#' @param transcr.obj A transcriptomic object
#' @param connect.obj A connectomic object, ideally paired with / derived from the transcriptomic object.
#' @param feature A connectomic feature (signaling mechanism), with and 'em dash' ("—") as separator
#' @param plot.function The base function to use. Currently only option is 'ggCircuit', which is a custom sub-function in NICHESMethods, based on trigonometry + ggplot2. May also be able to build an iGraph-based version in the future.
#' @param group.by Character string. Name of the grouping variable to use. Must be in the metadata of the transcriptomic object.
#' @param graph.angle The orientation angle of the whole CircuitPlot graph, rotational.
#' @param h Arrow start and end offset from node center (the amount the arrows are shortened on center.)
#' @param offset Spacing between opposite direction paracrine arrows, from common center line.
#' @param autocrine.offset #### MSBR not clear yet what this does, or if it is optimal approach
#' @param edge.scale.factor #### Bigger makes arrows thinner, only if edge.fixed.size == F #### Need to polish our approach to edge thickness.
#' @param arrow.head.angle Angle of the arrow head segments relative to arrow path
#' @param arrow.head.length Length of the arrow head lines
#' @param autocrine.arrow.curvature #### MSBR unclear how this operates, in need of polish to approach here
#' @param cols.use Named color palette / character vector. Colors to use for the nodes.
#' @param edge.fixed.size Default TRUE. Whether to fixe edge thickness, or have be proportional to connectivity strength, along with alpha.
#' @param global.node.list A list of nodes to layout within the graph. Useful for placing dummy nodes for cross-system comparison.
#' @param min.edge.value Default 0. The minimum represented value for connectivity strength.
#' @param max.edge.value Default NULL. The maximum represented value for connectivity strength. Allows user-defined global scale for multi-plot comparisons. ### MSBR may want to update how we handle this, may want to generalize better.
#' @param unity.normalize If TRUE, will scale edge values (centered) and unity-normalize so that the numbers range from 0 to 1. If FALSE, no change is made, and the numbers represented by the edge weight in the graph are the raw aggregated connectivity values between nodes. The graph will look the same in both cases; this only influences the actual numbers being represented.
#' @return A circuit plot (ggplot object if plot.function == 'ggCircuit')

CircuitPlot <- function(transcr.obj,
                        connect.obj,
                        feature,
                        plot.function = 'ggCircuit',
                        group.by = 'CellClass',
                        graph.angle = 45,
                        h = 0.2,
                        offset = 0.05,
                        autocrine.offset = 0.03,
                        edge.scale.factor = 20,
                        arrow.head.angle = 15,
                        arrow.head.length = 0.03,
                        autocrine.arrow.curvature = 8,
                        cols.use = RColorBrewer::brewer.pal(4, 'Set2'),
                        edge.fixed.size = TRUE,
                        global.node.list = c('Endothelium', 'Epithelium', 'Mesenchyme', 'Immune'),
                        min.edge.value = 0,
                        max.edge.value = NULL,
                        unity.normalize = F,
                        ...) {
  # Define Objects
  ligand <- strsplit(feature, split = "—")[[1]][1]
  receptor <- strsplit(feature, split = "—")[[1]][2]
  node.object <- DefineNodeObject(transcr.obj = transcr.obj,
                                  feature = ligand)
  edge.object <- DefineEdgeObject(connect.obj = connect.obj,
                                  feature = feature)

  # Aggregate Objects
  group.by.edge <- paste(group.by, 'Joint', sep = '.')
  node.aggregate <- AggregateNodeData(node.object = node.object,
                                      group.by = group.by,
                                      global.node.list = global.node.list)
  edge.aggregate <- AggregateEdgeData(edge.object = edge.object,
                                      group.by = group.by.edge)

  # Adding in an edge-aggregate scaling feature for more pronounced/specific signals to be highlighted (SEE, 12/03/2024)
 if (unity.normalize == T) {
 edge.aggregate$feature.value <- scale(edge.aggregate$feature.value, center = TRUE, scale = TRUE)
 edge.aggregate$feature.value <- scales::rescale(edge.aggregate$feature.value, to = c(0, 1))
}

  # Plot with desired plot.function
  if (plot.function == 'ggCircuit') {
    circuit.plot <- ggCircuit(node.aggregate = node.aggregate,
                              edge.aggregate = edge.aggregate,
                              graph.angle = graph.angle,
                              h = h,
                              offset = offset,
                              autocrine.offset = autocrine.offset,
                              edge.scale.factor = edge.scale.factor,
                              arrow.head.angle = arrow.head.angle,
                              arrow.head.length = arrow.head.length,
                              autocrine.arrow.curvature = autocrine.arrow.curvature,
                              cols.use = cols.use,
                              edge.fixed.size = edge.fixed.size,
                              min.edge.value = min.edge.value,
                              max.edge.value = max.edge.value,
                              title = feature)
  }

  return(circuit.plot)
}
