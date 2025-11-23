# creation of test data for NICHESv2
setwd("~/GitHub/NICHESv2-msbr/tests/testthat/fixtures")

# source data file
load("~/Large Files/Pneumonectomy_Project_transfer/MSBR_polishing_2025-11-01/outputs/PPLR.polished.with.territory.and.neighborhood.2025-11-17.Robj")

# subset to one puck
sub <- subset(PPLR.polished, subset = TMA == 'PNX3_M2_Sp')
rm(PPLR.polished)

# define small crop
cropped.coords <- Crop(sub[["PNX3_M2_Sp"]], 
                       x = c(6900, 7075), 
                       y = c(4450, 4600), 
                       coords = "plot")

# set a new field of view (fov)
sub[["demo"]] <- cropped.coords

# inspect
ImageDimPlot(sub, 
             fov = 'demo',
             border.color = NA,
             axes=T,
             group.by = 'Type.6')+
  NoLegend()

# isolate demo into solitary object
sub$demo.cells <- NA
cropped.cells <- sub@images$demo$centroids@cells
sub$demo.cells[cropped.cells] <- 'demo'
demo <- subset(sub, subset = demo.cells == 'demo')
demo

ImageDimPlot(demo,
             fov = 'demo',
             border.color = NA,
             axes=T,
             group.by = 'Type.6')+
  NoLegend()

demo <- UpdateSeuratObject(demo)
DefaultAssay(demo) <- 'Xenium'

# clean up
# remove extraneous assays
names(demo@assays)
demo[['N2C_rad_raw']] <- NULL
demo[['C2N_rad_raw']] <- NULL
demo[['Neighborhood_Type5_Rad25']] <- NULL
names(demo@assays)

# remove extraneous graphs
names(demo@graphs)
demo[['N2C_rad_norm_nn']] <- NULL
demo[['N2C_rad_norm_snn']] <- NULL
demo[['N2C_knn_norm_nn']] <- NULL
demo[['N2C_knn_norm_snn']] <- NULL
demo[['N2C_rad_raw_nn']] <- NULL
demo[['N2C_rad_raw_snn']] <- NULL
demo[['N2C_knn_raw_nn']] <- NULL
demo[['N2C_knn_raw_snn']] <- NULL
demo[['C2N_rad_norm_nn']] <- NULL
demo[['C2N_rad_norm_snn']] <- NULL
demo[['C2N_knn_norm_nn']] <- NULL
demo[['C2N_knn_raw_nn']] <- NULL
demo[['C2N_rad_raw_nn']] <- NULL
demo[["C2N_knn_raw_snn"]] <- NULL    
demo[['C2N_rad_raw_snn']] <- NULL
demo[['N2C2N_rad_norm_nn']] <- NULL
demo[["Xenium_nn"]] <- NULL           
demo[["Xenium_snn"]] <- NULL          
demo[["C2N_knn_norm_snn"]] <- NULL    
demo[["N2C2N_rad_norm_snn"]] <- NULL  
demo[["N2C2N_knn_norm_nn"]] <- NULL   
demo[["N2C2N_knn_norm_snn"]] <- NULL 
demo[["N2C2N_rad_raw_nn"]] <- NULL
demo[["N2C2N_rad_raw_snn"]] <- NULL   
demo[["N2C2N_knn_raw_nn"]] <- NULL    
demo[["N2C2N_knn_raw_snn"]] <- NULL   
demo[["neighborhoodPCA.KNN"]] <- NULL 
demo[["neighborhoodPCA.SNN"]] <- NULL
names(demo@graphs)


# remove extraneous fovs
demo[['fov']] <- NULL
demo[['PNX3_M2_Sp']] <- NULL

# remove extraneous embeddings
demo[['pca']] <- NULL
demo[['umap']] <- NULL
demo[['pca_N2C_rad_raw_pc30_res0.6']] <- NULL
demo[['umap_N2C_rad_raw_pc30_res0.6']] <- NULL
demo[['pca_C2N_rad_raw_pc30_res0.6']] <- NULL
demo[['umap_C2N_rad_raw_pc30_res0.6']] <- NULL
demo[['neighborhoodPCA']] <- NULL
demo[['neighborhoodUMAP']] <- NULL

DefaultAssay(demo) <- 'Xenium'
demo <- UpdateSeuratObject(demo)

ImageDimPlot(demo,
             fov = 'demo',
             border.color = NA,
             axes=T,
             group.by = 'Class.6')

# save as seurat object
save(demo,file = 'demo.Seurat.Robj')

