# load demo data
load("~/GitHub/NICHESv2-msbr/tests/testthat/fixtures/demo.Seurat.Robj")

# niches.obj

## node.object

### node.list
node.list <- colnames(demo)

### node.meta.data
node.meta.data <- demo@meta.data

### node.data.counts
node.data.counts <- demo@assays$Xenium@layers$counts
rownames(node.data.counts) <- rownames(demo)

### node.data.normalized
node.data.normalized <- demo@assays$Xenium@layers$data
rownames(node.data.normalized) <- rownames(demo)

### node.data.scaled
node.data.scaled <- demo@assays$Xenium@layers$scale.data
rownames(node.data.scaled) <- rownames(demo)

## mechanism.object
### ground.truth.raw
### ground.truth.filtered


## edge.object

### edge.list

### edge.meta.data

### edge.data.counts

### edge.data.normalized

### edge.data.scaled

### edge.list.params




## outgoing (edges aggregated by sender)
### outgoing.list
### outgoing.meta.data
### outgoing.data.counts
### outgoing.data.normalized
### outgoing.data.scaled

## incoming (edges aggregated by receiver)

## local (edges aggregated within neighborhood of each node)



