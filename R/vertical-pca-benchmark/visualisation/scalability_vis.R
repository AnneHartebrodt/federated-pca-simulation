require(data.table)
require(ggplot2)
require(dplyr)
require(tidyr)

dat<-fread('/home/anne/Documents/featurecloud/pca/vertical-pca/results/results_for_david//wide.tikz.transmission.transmission_cost.tsv')

da<- dat %>% pivot_longer(-counter)
da<-as.data.table(da)
da$value <- da$value*4/(1000000)


ggplot(da, aes(name, value))+geom_boxplot()+theme(axis.text.x = element_text(angle=90))+
  theme_bw()+
  theme(axis.line=element_line(),
      legend.position = 'none',
      strip.background = element_blank(),
      panel.grid = element_blank(),
      axis.text = element_text(size=7))+
  ylab('transmitted data volume in MB')+
  xlab('percentage of samples')+
  ggtitle('Net volume of transmitted data')

scale_color_manual(values = viridis(4)[1:3])

data[, name_qr := paste0(matrix, eigenvector_update, qr_method, orthonormalisation_skip)]
data$rank <-as.factor(data$rank)
data$name_qr<-as.factor(data$name_qr)
data$name_qr<- recode(data$name_qr, matrixpowerfederated_qr1='FED-GS', matrixpowerfederated_qr100='NO-GS', vectorgradientcentral_qr1='GUO')
data <- as.data.table(data)
