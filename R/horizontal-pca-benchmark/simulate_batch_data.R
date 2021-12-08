library("splatter")
library("scater")
library("ggplot2")
library(zellkonverter)
require(BiocSingular)

data_path<-'/home/anne/Documents/featurecloud/pca/horizontal-pca/data/simulated/'

dir.create(data_path)
params <- newSplatParams()
# Set the number of genes to 1000
params <- setParam(params, "nGenes", 1000)

bc=500
sim.groups <- splatSimulate(params, batchCells = rep(bc, 5), group.prob = c(0.5, 0.5),
                            method = "groups", verbose = FALSE,
                            out.prob = 0.01,
                            batch.facLoc = 0.001, batch.facScale = 0.001,
                            )
sim.groups <- logNormCounts(sim.groups)
sim.groups <- runPCA(sim.groups, ntop=1000, scale=TRUE)
gp1<- plotPCA(sim.groups, shape_by = "Batch", colour_by = "Group")
ggsave(gp1, file=file.path(data_path, 'major_group.pdf'))
writeH5AD(sim.groups, file=file.path(data_path, 'major_group.h5ad'), X_name = 'logcounts')
gp1


##############
bc=500
sim.groups <- splatSimulate(params,batchCells = rep(bc, 5), group.prob = c(0.5, 0.5),
                            method = "groups", verbose = FALSE,
                            out.prob = 0.01,
                            batch.facLoc = 0.1, batch.facScale = 0.25,
)
sim.groups <- logNormCounts(sim.groups)
sim.groups <- runPCA(sim.groups, ntop=1000, scale=TRUE)
gp2<- plotPCA(sim.groups, shape_by = "Batch", colour_by = "Group")
gp2
ggsave(gp2, file=file.path(data_path, 'major_site.pdf'))
writeH5AD(sim.groups, file=file.path(data_path, 'major_site.h5ad'), X_name = 'logcounts')


##############
bc=500
params <- newSplatParams()
# Set the number of genes to 1000
params <- setParam(params, "nGenes", 1000)
sim.groups <- splatSimulate(params,batchCells = rep(bc, 5), group.prob = c(0.5, 0.5),
                            method = "groups", verbose = FALSE,
                            out.prob = 0.5, out.facLoc=1, out.facScale=1,
                            #batch.facLoc = 0.1, batch.facScale = 0.1,
                            batch.facLoc = 0.001, batch.facScale = 0.001,
)
sim.groups <- logNormCounts(sim.groups)
sim.groups <- runPCA(sim.groups, ntop=1000, scale=TRUE,BSPARAM=IrlbaParam(), exprs_values='logcounts')
gp3<-plotPCA(sim.groups, shape_by = "Batch", colour_by = "Group")
gp3

ggsave(gp3, file=file.path(data_path, 'major_group_outliers.pdf'))
writeH5AD(sim.groups, file=file.path(data_path, 'major_group_outlier.h5ad'), X_name = 'logcounts')



