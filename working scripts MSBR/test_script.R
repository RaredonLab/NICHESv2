# This is a test script for ToolBuilding

# Set WD
setwd("/Users/msbr/GitHub/NICHESMethods/R Scripts/ToolBuilding")

# Load transcriptomic data
load("~/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/My Drive/Tuft_Sox9_Pneumonectomy_Project/Homeostatic_Single_Cell/lung.combined.clean.classed.annotated.final.2023-01-29.Robj")
lung.combined

# Load connectomic data
load("~/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/My Drive/Tuft_Sox9_Pneumonectomy_Project/Homeostatic_Single_Cell/cell.to.cell.Robj")
cell.to.cell

# test DefineNodeObject [CHECKED]
node.object <- DefineNodeObject(transcr.obj = lung.combined,
                 feature = 'Tgfb1')
class(node.object)
dim(node.object)
names(node.object)
head(node.object)

# test DefineEdgeObject [CHECKED]
edge.object <- DefineEdgeObject(connect.obj = cell.to.cell,
                                feature = 'Tgfb1—Tgfbr1')
class(edge.object)
dim(edge.object)
names(edge.object)
head(edge.object)

# test AggregateNodeData [CHECKED]
node.aggregate <- AggregateNodeData(node.object = node.object,
                                    group.by = 'CellClass')
# test AggregateEdgeData [CHECKED]
edge.aggregate <- AggregateEdgeData(edge.object = edge.object,
                                    group.by = 'CellClass.Joint')

# test ggCircuit
ggCircuit(node.aggregate = node.aggregate,
          edge.aggregate = edge.aggregate,
          offset = 0.05,
          edge.fixed = T,
          edge.scale.factor = 20)

CircuitPlot(transcr.obj = lung.combined,
            connect.obj = cell.to.cell,
            feature = 'Vegfa—Kdr',
            group.by = 'CellClass',
            graph.angle = 45,
            h = 0.2,
            offset = 0.05,
            autocrine.offset = 0.02,
            edge.scale.factor = 20,
            arrow.head.angle = 15,
            arrow.head.length = 0.03,
            autocrine.arrow.curvature = 10,
            cols.use = RColorBrewer::brewer.pal(4,'Set2'),
            edge.fixed = T)

CircuitPlot(transcr.obj = lung.combined,
            connect.obj = cell.to.cell,
            feature = 'Osm—Osmr',
            group.by = 'CellClass',
            edge.fixed = F)

CircuitPlot(transcr.obj = lung.combined,
            connect.obj = cell.to.cell,
            feature = 'Angpt1—Tie1',
            group.by = 'CellClass',
            edge.fixed = T)

# test on allie's dataset
load("~/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/.shortcut-targets-by-id/1VLqBlyzO-Qad5O2kbwXkBRh_1cQho39t/Engineered Lung Paper/Global_Connectomics/global.connectomics.2023-11-11.Robj")
load("~/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/.shortcut-targets-by-id/1VLqBlyzO-Qad5O2kbwXkBRh_1cQho39t/Engineered Lung Paper/Global_Connectomics/downsampled.engOnly.customEmbed.2023-12-14.Robj") # 'sub'
# First we have to scale the connectivity data, just in case we haven't before
connectivity <- global.connectomics[['CellToCell']][['alra']]
connectivity <- ScaleData(connectivity)
# Build a new metadata slot ('Dataset4') to split based on M/C/T/Q
Idents(sub) <- sub$Dataset2
Idents(connectivity) <- connectivity$Dataset2.Sending
table(Idents(sub))
table(Idents(connectivity))
sub <- RenameIdents(sub,
                    'Tri_E'='Tri',
                    'Tri_L'='Tri',
                    'Quad_E'='Quad',
                    'Quad_L'='Quad')
connectivity <- RenameIdents(connectivity,
                             'Tri_E'='Tri',
                             'Tri_L'='Tri',
                             'Quad_E'='Quad',
                             'Quad_L'='Quad')
sub$Dataset4 <- Idents(sub)
connectivity$Dataset4 <- Idents(connectivity)

sub$Dataset4 <- factor(sub$Dataset4,levels=c('Mono','Co','Tri','Quad'))
connectivity$Dataset4 <- factor(connectivity$Dataset4,levels=c('Mono','Co','Tri','Quad'))

# Divide up into 4 systems
Idents(sub) <- sub$Dataset4
tran.mono <- subset(sub,idents = 'Mono')
tran.co <- subset(sub,idents = 'Co')
tran.tri <- subset(sub,idents = 'Tri')
tran.quad <- subset(sub,idents = 'Quad')
Idents(connectivity) <- connectivity$Dataset4
con.mono <- subset(connectivity,idents = 'Mono')
con.co <- subset(connectivity,idents = 'Co')
con.tri <- subset(connectivity,idents = 'Tri')
con.quad <- subset(connectivity,idents = 'Quad')

feature = 'Vegfa—Kdr'
feature = 'Tgfb1—Tgfbr1'
feature = 'Sct—Sctr'
feature = 'Fgf1—Egfr'

TempFunc <- function(feature,
                     max.edge.value = 0.1){
  p4 <- CircuitPlot(transcr.obj = tran.quad,
              connect.obj = con.quad,
              feature = feature,
              group.by = 'class',
              edge.fixed.size = T,min.edge.value = 0,max.edge.value = max.edge.value,cols.use = cols.use)
  p3 <- CircuitPlot(transcr.obj = tran.tri,
              connect.obj = con.tri,
              feature = feature,
              group.by = 'class',
              edge.fixed.size = T,min.edge.value = 0,max.edge.value = max.edge.value,cols.use = cols.use)
  p2 <- CircuitPlot(transcr.obj = tran.co,
              connect.obj = con.co,
              feature = feature,
              group.by = 'class',
              edge.fixed.size = T,min.edge.value = 0,max.edge.value = max.edge.value,cols.use = cols.use)
  p1 <- CircuitPlot(transcr.obj = tran.mono,
              connect.obj = con.mono,
              feature = feature,
              group.by = 'class',
              edge.fixed.size = T,min.edge.value = 0,max.edge.value = max.edge.value,cols.use = cols.use)
  print(cowplot::plot_grid(p1,p2,p3,p4,nrow = 1))
}

# For allie 2024-02-16
cols.use <- RColorBrewer::brewer.pal(4,'Set2')
names(cols.use) <- global.node.list
scales::show_col(cols.use,ncol=4)
# create data that references the palette
colorKey = data.frame(colorName=names(cols.use))
colorKey$colorName <- factor(colorKey$colorName,levels = global.node.list)
# plot with ggplot, referencing the palette
ggplot(data=colorKey, aes(x=1, y = nrow(colorKey):1, fill=colorName, label=colorName)) +
  geom_tile() +
  scale_fill_manual(values = cols.use) +
  theme_void()+
  theme(legend.position="none") + 
  geom_text(size = 10)


setwd("/Users/msbr/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/.shortcut-targets-by-id/1VLqBlyzO-Qad5O2kbwXkBRh_1cQho39t/Engineered Lung Paper/Global_Connectomics")

# Epi autocrine
feature.list <- c('Tnf—Tnfrsf1b','Sct—Sctr','Wnt7a—Lrp6',
                  'Npnt—Itgb1','Ereg—Egfr','Wnt7b—Lrp5',
                  'Hbegf—Egfr','Lamb3—Itga3','Areg—Egfr')
max.edge.value.list <- c(0.25,0.05,0.1,
                         0.1,0.1,0.1,
                         0.25,1,0.1)
for(i in 1:length(feature.list)){
  png(filename = paste(feature.list[i],'.png',sep=''),width = 10,height = 2.5,units = 'in',res=300)
  TempFunc(feature = feature.list[i],
           max.edge.value = max.edge.value.list[i])
  dev.off()
}

# Mes to Epi
feature.list <- c('Fgf1—Fgfr3','Fgf1—Egfr','Lama1—Itga2','Csf2—Csf2ra',
                  'Ngf—Maged1','Wnt16—Fzd6','Fgf7—Fgfr3','Fn1—Itgb6',
                  'Vegfa—Ephb2','Wnt5a—Lrp5','Col4a1—Itga2','Tnc—Egfr')
max.edge.value.list <- c(0.05,0.05,0.5,0.5,
                         0.1,0.25,0.1,1,
                         0.5,0.5,0.5,0.1)
for(i in 1:length(feature.list)){
  png(filename = paste(feature.list[i],'.png',sep=''),width = 8,height = 2,units = 'in',res=300)
  TempFunc(feature = feature.list[i],
           max.edge.value = max.edge.value.list[i])
  dev.off()
}

