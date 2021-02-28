require(data.table)
require(ggplot2)
require(dplyr)
require(tidyr)
require("ggrepel")

source('/home/anne/Documents/featurecloud/pca/federated_dp_pca/R/horizontal-pca-benchmark/library.R')

outdir <- '/home/anne/Documents/featurecloud/pca/horizontal-pca/results/accuracy/merged/5/'
data<- read_cancer_types(outdir)
data$splits<-5

outdir <- '/home/anne/Documents/featurecloud/pca/horizontal-pca/results/accuracy/merged/2/'
data2<- read_cancer_types(outdir)
data2$splits<-2

outdir <- '/home/anne/Documents/featurecloud/pca/horizontal-pca/results/accuracy/pre_split/'
data3<- read_cancer_types(outdir)
data3$splits<-'TCGA'

data<-rbind(data, data2, data3)


ggplot(data[algorithm=='balcan_proxy' & name %in% paste0('', 1:10)], aes(name, value, fill=as.factor(splits)))+
  geom_boxplot(outlier.shape = NA)+
  xlab('Eigenvalue rank')+
  scale_shape_manual(values=c(1, 8))+
  scale_color_manual(values = c('#000000','#FF0000'))+
  guides(color=F)+
  my_theme+
  scale_fill_manual('#Sites', values = palette_div[c(1,4,8)])+
  ylab('Angle w.r.t reference [degree]')
