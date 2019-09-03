require(data.table)
require(ggplot2)
require(RColorBrewer)
require(viridisLite)
require(colorspace)
library(ggbiplot)
require(Matrix)
require(ggforce)

data <- fread('/home/anne/Documents/featurecloud/data/TCGA/htseq/simulation/eigenvalues/nonoise_projections.txt')
cn<-1:(ncol(data)-4)
cn<-paste0('PC', cn)


nrdat<-nrow(unique(data[,c(ncol(data)-3,ncol(data)-2,ncol(data)-1,ncol(data)), with=F]))
indx<-nrow(data)/nrdat

data$index <- rep(1:indx,nrdat)

colnames(data)<-c(cn, "split", "re", "epsilon", "delta", "index")
sliderSplit <- unique(data$split)
sliderEpsilon<- unique(data$epsilon)
sliderDelta <- unique(data$delta)