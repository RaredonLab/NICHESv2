CustomHeatmapDendrogram <- function(object,
                               data.type = 'CellToCell',
                               primary = 'seurat_clusters' ,
                               secondary = 'SendingType' ,
                               tertiary = 'ReceivingType' ,
                               #quarternary = 'orig.ident' ,
                               primary.cols = NULL,
                               secondary.cols = NULL, # Need to be a named list of colors
                               tertiary.cols = NULL,
                               #quarternary.cols = NULL,
                               features = NULL,
                               labels = NULL,
                               selected.row.anotations=NULL,
                               selected.label.size = 10,
                               use.scale.data = T,
                               range.frac = 1,
                               row.dendrogram = T){
  # Packages
  require(tidyverse)
  require(RColorBrewer)
  require(ggplot2)
  require(ComplexHeatmap)
  
  # Local functions
  source("~/GitHub/Engineered-Lung-Analysis/SamCode/General Functions/gg_color_hue.R")
  
  # add randomization column
  object$random <- sample(ncol(object))
  
  # Stash
  focus <- object
  
  # MOI
  MOI <- features
  
  # Get and organize data for fancy heatmap
  meta.data <- focus@meta.data
  meta.data$barcode <- rownames(meta.data)
  meta.data <- meta.data[order(
    meta.data[[primary]],
    meta.data[[secondary]], # TEMPORARY EXPERIMENT FOR VISUALS - THIS CONTROLS THE ORDER OF THE BARS
    meta.data[[tertiary]],
    meta.data[['random']]),] # Experiment 2024-04-18 to randomize order
  if(use.scale.data == T){
    to.plot <- as.matrix(focus@assays[[data.type]]@scale.data[MOI,meta.data$barcode])
  }else{
    to.plot <- as.matrix(focus@assays[[data.type]]@data[MOI,meta.data$barcode])
  }
  
  # Rownames and Column names
  rownames(to.plot) <- MOI
  colnames(to.plot) <- NULL
  
  # Colors
  # 1
  if(is.null(primary.cols)){
    primary.colors <- gg_color_hue(length(unique(meta.data[[primary]])))
    names(primary.colors) <- unique(meta.data[[primary]])
  }else{
    primary.colors <- primary.cols[unique(meta.data[[primary]])]
    names(primary.colors) <- unique(meta.data[[primary]])
  }
  # 2
  if(is.null(secondary.cols)){
    cols.2 <- RColorBrewer::brewer.pal(n=9,name='Set1')
    secondary.colors <- colorRampPalette(cols.2)(length(unique(meta.data[[secondary]])))
    names(secondary.colors) <- unique(meta.data[[secondary]])
  }else{
    secondary.colors <- secondary.cols[unique(meta.data[[secondary]])]
    names(secondary.colors) <- unique(meta.data[[secondary]])
  }
  # 3
  if(is.null(tertiary.cols)){
    cols.3 <- RColorBrewer::brewer.pal(n=8,name='Accent')
    tertiary.colors <- colorRampPalette(cols.3)(length(unique(meta.data[[tertiary]])))
    names(tertiary.colors) <- unique(meta.data[[tertiary]])
  }else{
    tertiary.colors <- tertiary.cols[unique(meta.data[[tertiary]])]
    names(tertiary.colors) <- unique(meta.data[[tertiary]])
  }
  # # 4
  # if(is.null(quarternary.cols)){
  #   cols.4 <- RColorBrewer::brewer.pal(n=12,name='Paired')
  #   quarternary.colors <- colorRampPalette(cols.4)(length(unique(meta.data[[quarternary]])))
  #   names(quarternary.colors) <- unique(meta.data[[quarternary]])
  # }else{
  #   quarternary.colors <- quarternary.cols[unique(meta.data[[quarternary]])]
  #   names(quarternary.colors) <- unique(meta.data[[quarternary]])
  # }
  
  # Annotation [HOW DO I GENERALIZE?]
  # Define annotations
  stuff <- data.frame(primary = meta.data[[primary]],
                      secondary = meta.data[[secondary]],
                      tertiary = meta.data[[tertiary]],
                      #quartenary = meta.data[[quarternary]],
                      check.names = F)
  names(stuff) <- labels
  
  # Define colors
  colors <- list(
    primary = primary.colors,
    secondary = secondary.colors,
    tertiary = tertiary.colors)#,
  #quarternary = quarternary.colors)
  names(colors) <- labels
  
  column_ha <- ComplexHeatmap::HeatmapAnnotation(
    df = stuff,
    col = colors,
    annotation_name_side = "left")
  
  # Value colors
  # col_fun = circlize::colorRamp2(c(-2, 0, 2), c("grey", "white", "blue"))
  if(use.scale.data == T){
    col_fun = circlize::colorRamp2(c(min(to.plot)*range.frac,0, max(to.plot)*range.frac), c("#440154",'#21918c',"#fde725"))
    viridis::scale_fill_viridis() 
  }else{
    col_fun = circlize::colorRamp2(c(0,max(to.plot)*range.frac), c("#440154", "#fde725"))
  }
  
  # Legend title
  if(use.scale.data == T){
    legend.title = 'Scaled Connectivity'
  }else{
    legend.title = 'Connectivity'
  }
  # Heatmap with selected rows annotated
  if(!is.null(selected.row.anotations)){
    # Set up row annotations
    #htmp_input <- nlmeganoNA.small@assays$RNA@scale.data[nlmega.markers_res0.8$gene, ] #create matrix for heatmap input
    #ShowGenes <- nlmega.markers_res0.8$gene[c(1, 3, 4, X, etc...)] #select genes that you want to emphasize out of total gene list
    HAleft <- rowAnnotation(foo = anno_mark(at = which(rownames(to.plot) %in% selected.row.anotations), side = 'left',
                                            labels = rownames(to.plot)[rownames(to.plot) %in% selected.row.anotations],
                                            labels_gp = gpar(fontsize=selected.label.size)))#,
    # extend = unit(0, "mm"), #not sure if this is bringing anything to the table
    # link_width = unit(5, "mm"), #this actually works and sets the width of the lines
    # link_gp = gpar(lineheight = X))) #not sure if this works, but I tried to set line height like this.. X = multiple of labels font size??
    row.names.stash <- row.names(to.plot)
    row.names(to.plot) <- NULL
    # Make heatmap
    heatmap.object <- Heatmap(to.plot,
                              column_split = stuff[,1], # breaks on primary
                              column_gap = unit(0.3,'mm'),
                              col = col_fun,
                              use_raster = F,
                              cluster_rows = row.dendrogram,
                              cluster_columns = F,
                              show_column_dend = F,
                              show_row_dend = row.dendrogram,
                              top_annotation = column_ha,
                              #row_names_side = 'left',
                              #row_names_gp = gpar(fontsize = 10),
                              color_space = colors.inferno,
                              name = legend.title,
                              heatmap_legend_param = list(title_position = 'leftcenter-rot'),
                              row_title = 'Signaling Mechanisms',
                              column_title = NULL,
                              left_annotation = HAleft)
    draw(heatmap.object,
         heatmap_legend_side='right',
         annotation_legend_side='right')
  }else{
    # OR Heatmap with all rows annotated
    heatmap.object <- Heatmap(to.plot,
                              column_split = stuff[,1], # breaks on primary
                              column_gap = unit(0.3,'mm'),
                              col = col_fun,
                              use_raster = F,
                              cluster_rows = row.dendrogram,
                              cluster_columns = F,
                              show_column_dend = F,
                              show_row_dend = row.dendrogram,
                              top_annotation = column_ha,
                              row_names_side = 'left',
                              row_names_gp = gpar(fontsize = 10),
                              color_space = colors.inferno,
                              name = legend.title,
                              heatmap_legend_param = list(title_position = 'leftcenter-rot'),
                              row_title = 'Signaling Mechanisms',
                              column_title = NULL)
    draw(heatmap.object,
         heatmap_legend_side='right',
         annotation_legend_side='right')
  }
  
  row.order.output <<- row_order(heatmap.object)
}
