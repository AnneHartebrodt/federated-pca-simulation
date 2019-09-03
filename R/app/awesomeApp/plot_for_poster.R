require(data.table)
require(ggplot2)
require(RColorBrewer)
require(viridisLite)
require(colorspace)
library(ggbiplot)
require(Matrix)
require(ggforce)
require(cowplot)

selection<-c( 148,   5,  13,  104,  16,  124,  28,  30,  49,  34,
              1,83,8,63,10,117,150,43,52,98,
              3,90,23,42,58,128,86,110,147, 122)

data <- fread('/home/anne/Documents/featurecloud/featurecloud-test-data/PCA01/IRIS_simulation_multiple.tsv')
labels <- fread('/home/anne/Documents/featurecloud/featurecloud-test-data/PCA01/IRIS_labels.tsv', header = F)
data$label <- rep(labels$V1,500)
data$index<-rep(1:150,500)
data<-data[index %in% selection]

neworder<-data[1:30]
neworder<-neworder[order(neworder$label)]

data<-data[order(data$label)]
colnames(data)<-c('V1', "V2", 'V3', "V4", "split", "iteration", "epsilon", "delta", "label","index")

epsilon1 = 10
delta1 = 0.01
owners1 = 1

pp<-c(rep(inferno(20)[1:10],1), rep(viridis(20)[11:20],1),rep(plasma(20)[11:20],1))
pp<-pp[order(neworder$index)]
g4 <- ggplot(data[epsilon %in% c(0.01,0.1,1,10) & split %in% c(1,5)], aes(x=V1, y=V2))+
  ggtitle('PCA with different levels of privacy')+
  scale_shape_manual(values = c(15,17,19))+
  ylab('PC2')+xlab('PC1')+theme(legend.title = element_blank(), plot.title = element_text(size = 28))+
  #geom_mark_ellipse(aes(color=as.factor(index),alpha=0.001), expand = unit( 2, "mm"))+
  geom_point(aes(shape= as.factor(label), color = as.factor(index)))+
  theme(legend.position = 'bottom')+
  facet_grid(vars(split), vars(epsilon))+
  scale_color_manual(values = pp)
g4


ggsave(filename='/home/anne/Documents/ECCB19/poster/PCA_plot.pdf', g4, width = 8, height = 5)
