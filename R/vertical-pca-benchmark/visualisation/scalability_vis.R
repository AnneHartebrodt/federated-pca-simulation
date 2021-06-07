require(data.table)
require(ggplot2)
require(dplyr)
require(tidyr)

data<-fread('/home/anne/Documents/featurecloud/pca/vertical-pca/results-new-tests/scalability/summary.transmission_cost.tsv')

data$cost<- data$nr_float*4/(1000000)
#data$facet_title <- ordered(data$facet_title, levels=levels(data$rank))
data[, name_qr := paste0(matrix, eigenvector_update, qr_method, orthonormalisation_skip)]
data$rank <-as.factor(data$rank)
data$name_qr<-as.factor(data$name_qr)
data$name_qr<- recode(data$name_qr, matrixpowerfederated_qr1='FED-GS', matrixpowerfederated_qr100='NO-GS', vectorgradientcentral_qr1='GUO')
data <- as.data.table(data)



ggplot(data, aes(as.factor(data_fraction), max_it, fill=name_qr))+
  geom_boxplot()+
  theme_bw()+
  theme(axis.line=element_line(),
      legend.position = 'top',
      strip.background = element_blank(),
      panel.grid = element_blank(),
      axis.text = element_text(size=7))+
  ylab('transmitted data volume in MB')+
  xlab('percentage of samples')+
  ggtitle('Net volume of transmitted data')+
  scale_fill_manual(values = viridis(4)[1:3])

data[, name_qr := paste0(matrix, eigenvector_update, qr_method, orthonormalisation_skip)]
data$rank <-as.factor(data$rank)
data$name_qr<-as.factor(data$name_qr)
data$name_qr<- recode(data$name_qr, matrixpowerfederated_qr1='FED-GS', matrixpowerfederated_qr100='NO-GS', vectorgradientcentral_qr1='GUO')
data <- as.data.table(data)
