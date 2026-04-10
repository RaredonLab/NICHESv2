neighborhoods <- readRDS("~/Desktop/TEMP/my_spatial_network_neighborhoods_v3.rds")



source("/Users/nw426/Documents/GitHub/NICHESv2-NW/working-scripts-NW/Claude/extract_neighborhoods.R")
source("/Users/nw426/Documents/GitHub/NICHESv2-NW/working-scripts-NW/Claude/validate_neighborhoods.R")

# load("/Volumes/T7/Pneum/Xenium/PPLR.polished.with.territory.and.neighborhood.2025-11-17.Robj")
# count.mtx <- PPLR.polished@assays$Xenium$counts
# meta.data <- PPLR.polished@meta.data
# 
# new.PPLR.niches <- create_NICHESObject(count.mtx,
#                                        meta.data,
#                                        LRM.db = 'fantom5',
#                                        mode = "spatial",
#                                        species = 'rat',
#                                        sample.col = 'TMA',
#                                        rad = 25,
#                                        n.cores = 4)

# saveRDS(new.PPLR.niches, file = 'new.PPLR.niches.rds')

new.PPLR.niches <- readRDS("~/Desktop/NICHESv2-NW/new.PPLR.niches.rds")

my.result <- extract_neighborhoods(new.PPLR.niches$edge.list, n.cores = 8)

canonical.advisor <- advisor_to_canonical(neighborhoods)
canonical.mine <- niches_to_canonical(my.result)

compare_neighborhoods(canonical.advisor, canonical.mine)


