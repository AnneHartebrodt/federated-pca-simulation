require(data.table)
require(ggplot2)
require(RColorBrewer)
require(viridisLite)
require(colorspace)
library(ggbiplot)
require(Matrix)
require(ggforce)

data <- fread('/home/anne/Documents/featurecloud/results/simulation_breast/noise_pca.projection')
cn<-1:(ncol(data)-4)
cn<-paste0('PC', cn)

t<-table(data[,c(ncol(data)-3,ncol(data)-2,ncol(data)-1,ncol(data)), with=F])
nrdat<-nrow(unique(data[,c(ncol(data)-3,ncol(data)-2,ncol(data)-1,ncol(data)), with=F]))

id<-c()
for(el in rep(t[1:5],each=max(data$V12)+1)){
  id<-c(id, 1:el)
}

data$index<-id
colnames(data)<-c(cn, "split", "re", "epsilon", "delta", "index")
sliderSplit <- unique(data$split)
sliderEpsilon<- unique(data$epsilon)
sliderDelta <- unique(data$delta)
