NetworkPlot <- function(transcriptome.object, # scRNAseq object containing all populations within which edges will be plotted. The embedding of this object will determine the shape of the graph.
                               connectome.object, # CellToCell signaling object containing the edges of interest to be visualized as a graph.
                               mechanism.of.interest, # An individual L-R signaling mechanism of interest. Character vector of length 1.
                               legends.to.plot = NULL, # Metadata to use for legend plot(s). Character vector of length > 1. Default NULL. Determines number of rows (or columns if transposed) for joint output plot.
                               split.by = NULL, # Metadata handle to make a split plot. Character vector of length = 1. Default NULL. Determines number of columns (or rows if transposed) for joint output plot.
                               transpose = FALSE, # Whether to transpose plot organization or not. Default FALSE.
                               reduction = 'umap.rpca',
                               min.connectivity.thresh = 0,
                               point.size = 0.1,
                               alpha.min = 0,
                               alpha.max = 0.85,
                               line.thickness = 0.2,
                               arrowhead.size = 0,
                               connectivity.color.min = 0, # Sets lower bound on the connectivity color range
                               connectivity.color.max = NULL, # Sets upper bound on the connectivity color range
                               ligand.color.min = 0, # Sets lower bound on the ligand color range
                               receptor.color.min = 0, # Sets lower bound on the receptor color range
                               legend.palettes = NULL, # A named list of named palettes for the legend plots to use. The names must be the same as the elements of legends.to.plot
                               legend.alpha = 0.25, # Alpha for legend coloring
                               black.points = T,
                               node.assay = 'RNA',
                              edge.assay = 'CellToCell'
                               ){
  #### SETUP ####
  # Break mechanism of interest into ligand and receptor
  ligand <- strsplit(mechanism.of.interest,split='—')[[1]][[1]]
  receptor <- strsplit(mechanism.of.interest,split='—')[[1]][[2]]

  # Check if transcriptome.object contains the desired embedding
  if(reduction %in% names(transcriptome.object@reductions)){
    message(paste('Using',reduction,'reduction for transcriptome embedding'))
  }else{
    stop('transcriptome object does not contain the requested reduction for graph embedding. Please address.')
  }

  # Check if mechanism of interest is in the connectome.object
  if(!(mechanism.of.interest %in% rownames(connectome.object))){
    stop('mechanism of interest is not in the rownames of the connectome object. Please address.')
  }

  # Check if ligand an receptor information are in the transcriptome.object
  if(!(ligand %in% rownames(transcriptome.object))){
    stop('ligand of interest is not in the rownames of the transcriptome object. Please address.')
  }
  if(!(receptor %in% rownames(transcriptome.object))){
    stop('receptor of interest is not in the rownames of the transcriptome object. Please address.')
  }

  # Sketch output plot structure
  if(transpose == FALSE){
    message(paste('Planning',length(legends.to.plot),'legend row(s) and 1 connectivity row'))
    if(!is.null(split.by)){
      message(paste('split.by = ',split.by,', which contains ',length(names(table(transcriptome.object[[split.by]]))),' categories',sep = ''))
    }
    message(paste('Total output plot grid will be',length(legends.to.plot)+1,'rows by 1 column'))
  }else{
    message(paste('Output plot will have',length(legends.to.plot),'legend columns and 1 connectivity column'))
    if(!is.null(split.by)){
      message(paste('split.by = ',split.by,', which contains ',length(names(table(transcriptome.object[[split.by]]))),' categories',sep = ''))
    }
    message(paste('Total output plot grid will be 1 row by',length(legends.to.plot)+1,'columns'))
  }

  # Gather data for transcriptomic plot layers
  message('Extracting transcriptomic information for plotting...')
    # Coordinates
    umap.coords <- data.frame(transcriptome.object@reductions[[reduction]]@cell.embeddings)
    # Convert rownames to row identifiers
    umap.coords$barcode <- rownames(umap.coords)
    # Metadata
    point.meta.data <- transcriptome.object@meta.data
    # LR.info
    rownames(transcriptome.object[[node.assay]]@layers$data) <- rownames(transcriptome.object)
    colnames(transcriptome.object[[node.assay]]@layers$data) <- colnames(transcriptome.object)
    LR.info <- data.frame(ligand.info = transcriptome.object@assays[[node.assay]]@layers$data[ligand,],
                          receptor.info = transcriptome.object@assays[[node.assay]]@layers$data[receptor,])
    # Assemble into dataframe for ggplotting
    if(sum(rownames(point.meta.data) != rownames(umap.coords))==0 &
       sum(rownames(LR.info) != rownames(umap.coords))==0 &
       sum(rownames(point.meta.data) != rownames(LR.info))==0){
      transcriptome.information <- cbind(umap.coords,LR.info,point.meta.data)
    }else{
      stop('the rownames for UMAP coords are different from the rownames for point metdata or LR.info, please address')
    }

  # Gather data for connectomic plot layers
    message('Extracting connectomic information for plotting...')
    # Connectivity
    rownames(connectome.object[[edge.assay]]@layers$data) <- rownames(connectome.object)
    colnames(connectome.object[[edge.assay]]@layers$data) <- colnames(connectome.object)
    connectivity <- data.frame(connectivity = connectome.object@assays$CellToCell@layers$data[mechanism.of.interest,])
    # Convert rownames to row identifiers
    connectivity$edge.ident <- rownames(connectivity)
    # Split the barcodes into sending vs. receiving:
    split.barcodes <- stringr::str_split_fixed(connectivity$edge.ident,"—",2)
    connectivity$SendingBarcode <- split.barcodes[,1]
    connectivity$ReceivingBarcode <- split.barcodes[,2]
    # Metadata
    edge.meta.data <- connectome.object@meta.data
    # Assemble
    connectome.information <- cbind(connectivity,edge.meta.data)

  # Compute the starting and ending points for each edge
  message('Computing starting and ending points for each edge...')
  sending.barcode.umap.coords <- data.frame(transcriptome.information[connectome.information$SendingBarcode,])
  receiving.barcode.umap.coords <- data.frame(transcriptome.information[connectome.information$ReceivingBarcode,])
  connectome.information$sending.barcode.umap1 <- sending.barcode.umap.coords[,1] # bad idea, works?
  connectome.information$sending.barcode.umap2 <- sending.barcode.umap.coords[,2] # bad idea, works?
  connectome.information$receiving.barcode.umap1 <- receiving.barcode.umap.coords[,1] # bad idea, works?
  connectome.information$receiving.barcode.umap2 <- receiving.barcode.umap.coords[,2] # bad idea, works?

  # Threshold edges
  message('Will only plot edges with >',min.connectivity.thresh,' connectivity, per user input...')
  downsampled <- connectome.information[connectome.information$connectivity>min.connectivity.thresh,]

  # Log transform connectivity values to be plotted
  message('Log transforming connectivity values for plotting...')
  downsampled$connectivity.to.plot <- log1p(downsampled$connectivity)

  # Set connectivity color range based on the data and/or user input
  if(!is.null(connectivity.color.max)){
    connectivity.color.max <- connectivity.color.max
  }else{
  connectivity.color.max <- max(downsampled$connectivity.to.plot)
  }
  if(!is.null(connectivity.color.min)){
    connectivity.color.min <- connectivity.color.min
  }else{
    connectivity.color.min <- min(downsampled$connectivity.to.plot)
  }

  # Set receptor color range based on data range or user input
  receptor.info.to.plot <- transcriptome.information[downsampled$ReceivingBarcode,] # not sure why this sometimes has NA in it?
  if(!is.null(receptor.color.min)){
    receptor.color.min <- receptor.color.min
  }else{
    receptor.color.min <- min(receptor.info.to.plot$receptor.info)
  }
  receptor.color.max <- max(receptor.info.to.plot$receptor.info,na.rm = T)

  # Set ligand color range based on data range or user input
  ligand.info.to.plot <- transcriptome.information[downsampled$SendingBarcode,] # not sure why this sometimes has NA in it?
  if(!is.null(ligand.color.min)){
    ligand.color.min <- ligand.color.min
  }else{
    ligand.color.min <- min(ligand.info.to.plot$ligand.info)
  }
  ligand.color.max <- max(ligand.info.to.plot$ligand.info,na.rm = T)

  #### BASE PLOT CONSTRUCTION ####
    # Base plot
    message('Plotting base plot...')
    base.plot <- ggplot(data.frame(transcriptome.information),
                   aes_q(x=sym(colnames(transcriptome.information)[1]),
                       y=sym(colnames(transcriptome.information)[2])))+ # RISKY GENERALIZATION
      theme_classic()

    # ggplot(data.frame(transcriptome.information),
    #        aes(x=umaprpca_1,
    #            y=umaprpca_2))+ # RISKY GENERALIZATION
    #   theme_classic()

  #### SPLIT.BY == NULL ####
    if(is.null(split.by)){

      # Connectivity plot
    message('Plotting edges...')
    connectivity.plot <- base.plot +
      geom_point(size=point.size,alpha=0.1,color='grey')+ # put a light-grey point base layer down first
      scale_alpha_continuous(range=c(alpha.min,alpha.max),limits=c(connectivity.color.min,connectivity.color.max),name='Connectivity')+ # sets the alpha range for the segments
      scale_colour_gradientn(colours = c('#4C1E4F','#348AA7','#FAA916','#EF233C'),limits=c(connectivity.color.min,connectivity.color.max),name='Connectivity')+
      geom_segment(data = downsampled[sample(rownames(downsampled),size=nrow(downsampled)),], # randomizes the segment plotting order
                   aes(x = sending.barcode.umap1,
                       y = sending.barcode.umap2,
                       xend = receiving.barcode.umap1,
                       yend = receiving.barcode.umap2,
                       color = connectivity.to.plot, # edges colored by connectivity
                       alpha = connectivity.to.plot), # lower value edges are more translucent
                   linewidth = line.thickness,
                   arrow = grid::arrow( # Adding this here for the first time, adding little tiny arrowheads to be able to better tell signaling direction
                     angle = 10, # how 'fat' the arrowheads are
                     length = unit(arrowhead.size, "npc"), # size of the arrowheads, including 0
                     ends = 'last',
                     type = 'closed'))

    # Transcriptomic plot
    message('Plotting transcriptome information...')

    # define ligand and receptor information
    receptor.info.to.plot <- transcriptome.information[downsampled$ReceivingBarcode,] # not sure why this sometimes has NA in it?
    ligand.info.to.plot <- transcriptome.information[downsampled$SendingBarcode,] # not sure why this sometimes has NA in it?

    # plot colored points
    if(black.points==F){
      # Add receptor expressivity

      output.plot <- connectivity.plot +
        ggnewscale::new_scale_color() +
        scale_colour_gradientn(colours = c('#FF000000','#FF000050'),limits=c(receptor.color.min,receptor.color.max),name = "Receptor Expression")+
        geom_point(data=receptor.info.to.plot,
                   aes(color = receptor.info.to.plot$receptor.info),#,
                       #alpha = receptor.info.to.plot$receptor.info),
                   size = point.size)

      # Add ligand expressivity

      output.plot <- output.plot +
        ggnewscale::new_scale_color() +
        scale_colour_gradientn(colours = c('#0000FF00','#0000FF50'),limits=c(ligand.color.min,ligand.color.max),name = "Ligand Expression")+
        geom_point(data=ligand.info.to.plot,
                   aes(color = ligand.info.to.plot$ligand.info),#,
                       #alpha = ligand.info.to.plot$ligand.info),
                   size = point.size)
    }
      # if black points ==T
      if(black.points==T){
        output.plot <- connectivity.plot +
          geom_point(data=transcriptome.information,
                     size = point.size)# use aes_() here
      }

      # Add connectivity plot title
      output.plot <- output.plot+ ggtitle(paste(mechanism.of.interest,'Connectivity'))+NoLegend()

      # Add legend(s)
      # 1. Set default legend colors if required
      if(sum(!(legends.to.plot %in% names(legend.palettes)))){
        warning('a legend palette has not been provided for every legend requested. Default colors will be used as needed.')
      }
      # 2. Make legend plot(s)
      if(!is.null(legends.to.plot)){
        legend.plot.list <- list()
        for(i in 1:length(legends.to.plot)){
          if(legends.to.plot[i] %in% names(legend.palettes)){
          legend.plot.list[[i]] <- Seurat::DimPlot(transcriptome.object,
                                                    group.by = legends.to.plot[i],
                                                    cols = legend.palettes[[legends.to.plot[i]]],
                                                    label=T)+NoLegend()+ggtitle(paste('Legend:',legends.to.plot[i]))
          }else{
            legend.plot.list[[i]] <- Seurat::DimPlot(transcriptome.object,
                                                     group.by = legends.to.plot[i],
                                                     #cols = legend.palettes[[legends.to.plot[i]]],
                                                     label=T)+NoLegend()+ggtitle(paste('Legend:',legends.to.plot[i]))
          }
        }

      }else{
        legend.plot.list <- list()
      }

      # Concatenate legend and connectivity plots together into output plot
      plot.list <- legend.plot.list
      plot.list[[length(plot.list)+1]] <- output.plot
      if(!is.null(legends.to.plot)){
        if(transpose==FALSE){
        output.plot.total <- cowplot::plot_grid(plotlist = plot.list,ncol = 1)
        }else{
          output.plot.total <- cowplot::plot_grid(plotlist = plot.list,nrow = 1)
        }
      }else{
        output.plot.total <- cowplot::plot_grid(plotlist = plot.list,nrow = 1)
      }

      # Test output

      output.plot.total
      # Return output total plot
    return(output.plot.total)
    }

  #### SPLIT.BY == TRUE ####
    if(!is.null(split.by)){

      # Test that user input makes sense
      if(length(split.by)>1){
        stop('This function only accepts x1 split.by parameter. Please revise input.')
      }

      # Define split.by.names
      split.by.names <- names(table(connectome.object[[split.by]]))

      # Split connectivity information
      downsampled.list <- list()
      for(i in 1:length(split.by.names)){
        downsampled.list[[i]] <- downsampled[downsampled[[split.by]]==split.by.names[i],]
      }
      names(downsampled.list) <- split.by.names

      # Split transcriptomic information
      transcriptomic.list <- list()
      for(i in 1:length(split.by.names)){
        transcriptomic.list[[i]] <- transcriptome.information[transcriptome.information[[split.by]]==split.by.names[i],]
      }
      names(transcriptomic.list) <- split.by.names

      # For each element of downsampled.list, make a connectivity plot with associated legends:
      # Initialize the strip-wise output
      plot.strip.list <- list()
      # Build each strip-wise output
      for(i in 1:length(split.by.names)){

      # Edge plotting on base
      connectivity.plot <- base.plot +
        geom_point(size=point.size,alpha=0.1,color='grey')+ # put a light-grey point base layer down first
        scale_alpha_continuous(range=c(alpha.min,alpha.max),limits=c(connectivity.color.min,connectivity.color.max),name='Connectivity')+ # sets the alpha range for the segments
        scale_colour_gradientn(colours = c('#4C1E4F','#348AA7','#FAA916','#EF233C'),limits=c(connectivity.color.min,connectivity.color.max),name='Connectivity')+
        geom_segment(data = downsampled.list[[i]][sample(rownames(downsampled.list[[i]]),size=nrow(downsampled.list[[i]])),], # randomizes the segment plotting order
                     aes(x = sending.barcode.umap1,
                         y = sending.barcode.umap2,
                         xend = receiving.barcode.umap1,
                         yend = receiving.barcode.umap2,
                         color = connectivity.to.plot, # edges colored by connectivity
                         alpha = connectivity.to.plot), # lower value edges are more translucent
                     linewidth = line.thickness,
                     arrow = grid::arrow( # Adding this here for the first time, adding little tiny arrowheads to be able to better tell signaling direction
                       angle = 10, # how 'fat' the arrowheads are
                       length = unit(arrowhead.size, "npc"), # size of the arrowheads, including 0
                       ends = 'last',
                       type = 'closed'))

      # Define ligand and receptor info
      receptor.info.to.plot <- transcriptomic.list[[i]][downsampled.list[[i]]$ReceivingBarcode,] # not sure why this sometimes has NA in it?
      ligand.info.to.plot <- transcriptomic.list[[i]][downsampled.list[[i]]$SendingBarcode,] # not sure why this sometimes has NA in it?

      # Plot black points or colored points

      if(black.points==F){
      # Add receptor expressivity
      output.plot <- connectivity.plot +
        ggnewscale::new_scale_color() +
        scale_colour_gradientn(colours = c('#FF000000','#FF000050'),limits=c(receptor.color.min,receptor.color.max),name = "Receptor Expression")+
        geom_point(data=receptor.info.to.plot,
                   aes(color = receptor.info.to.plot$receptor.info),#,
                   #alpha = receptor.info.to.plot$receptor.info),
                   size = point.size)

      # Add ligand expressivity
      output.plot <- output.plot +
        ggnewscale::new_scale_color() +
        scale_colour_gradientn(colours = c('#0000FF00','#0000FF50'),limits=c(ligand.color.min,ligand.color.max),name = "Ligand Expression")+
        geom_point(data=ligand.info.to.plot,
                   aes(color = ligand.info.to.plot$ligand.info),#,
                   #alpha = ligand.info.to.plot$ligand.info),
                   size = point.size)
      }
      # if black points ==T
      if(black.points==T){
        output.plot <- connectivity.plot +
          geom_point(data=transcriptomic.list[[i]],
                      size = point.size)# use aes_() here
      }

      # Add connectivity plot title
      output.plot <- output.plot+ ggtitle(paste(mechanism.of.interest,'Connectivity'))+NoLegend()

      # Add legend(s)
      # 1. Set default legend colors if required
      if(sum(!(legends.to.plot %in% names(legend.palettes)))){
        warning('a legend palette has not been provided for every legend requested. Default colors will be used as needed.')
      }
      # 2. Make legend plot(s) (needs to be done custom here with ggplot, cannot use Seurat this time)
      if(!is.null(legends.to.plot)){
        legend.plot.list <- list()
        for(j in 1:length(legends.to.plot)){
          if(legends.to.plot[j] %in% names(legend.palettes)){
            legend.plot.list[[j]] <- base.plot +
              geom_point(size=point.size,alpha=0.1,color='grey')+ # put a light-grey point base layer down first
              geom_point(data=transcriptomic.list[[i]],
                         aes_(color = transcriptomic.list[[i]][[legends.to.plot[j]]]), # use aes_() here
                         size = point.size,alpha = legend.alpha)+
              scale_color_manual(values = legend.palettes[[legends.to.plot[j]]])+
              ggtitle(paste('Legend:',legends.to.plot[j]))+NoLegend()
          }else{
            legend.plot.list[[j]] <- base.plot +
              geom_point(size=point.size,alpha=0.1,color='grey')+ # put a light-grey point base layer down first
              geom_point(data=transcriptomic.list[[i]],
                         aes_(color = transcriptomic.list[[i]][[legends.to.plot[j]]]),
                         size = point.size,alpha = legend.alpha)+
              ggtitle(paste('Legend:',legends.to.plot[j]))+NoLegend()
          }
        }

      }

      # Concatenate legend and connectivity plots together into sub output plots
      plot.list <- legend.plot.list
      plot.list[[length(plot.list)+1]] <- output.plot
      if(!is.null(legends.to.plot)){ # if legends requested
        if(transpose==FALSE){ # if transpose not requested
          plot.strip.list[[i]] <- cowplot::plot_grid(plotlist = plot.list,ncol = 1)
        }else{ # if transpose requested
          plot.strip.list[[i]] <- cowplot::plot_grid(plotlist = plot.list,nrow = 1)
        }
      }else{ # if no legends requested
        plot.strip.list[[i]] <- output.plot
        }
      } # end of single interation for a single split.by.name

      # Name the strips
      names(plot.strip.list) <- split.by.names

      # Concatenate all split.by strips together
      if(transpose==FALSE){ # if transpose not requested
      global.output <- cowplot::plot_grid(plotlist = plot.strip.list,nrow = 1)
      }else{# if transpose requested
        global.output <- cowplot::plot_grid(plotlist = plot.strip.list,ncol = 1)
      }
      # Return output total plot
      return(global.output)
    }

}

