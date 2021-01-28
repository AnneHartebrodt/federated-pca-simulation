require(data.table)
require(ggplot2)
require(tidyr)
require(cowplot)
require(optparse)
require(dplyr)
require(R.utils)
library(ggthemes)   
require(scales)
require(R.utils)
library(ggforce)
require(stringr)

palette_div<-c('#062865', '#2f497d', '#4c6d96', '#6793af', '#81bac8', '#ffc4b3', '#f68888', '#d75161', '#ab203f', '#720022')
palette_seq<-c('#062865', '#203a72', '#324d80', '#43618d', '#52759b', '#618aa9', '#70a0b7', '#7eb6c5', '#8dcdd4', '#9ce4e2')

my_theme <-
  theme_classic() + theme(
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 12),
    plot.title = element_text(size = 20, hjust = 0.5),
    axis.text = element_text(size = 8),
    legend.text = element_text(size = 8),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    legend.key.size = unit(1.5, 'lines'))

data <- fread('/home/anne/Documents/featurecloud/pca/vertical-pca/results/1000g/chr1/summaries/wide.vertical.angles_precomp.summary.tsv')
d <- data %>% pivot_longer(-c(iterations, rank))
d<-as.data.table(d)
selection<-c("matrix_2_power_central_qr", "matrix_2_power_federated_qr", "vector_2_power_central_qr", "vector_2_power_federated_qr","vector_2_gradient_central_qr")



p<-ggplot(d[name %in% selection], aes(iterations, value, col=as.factor(rank)))+geom_line()+facet_wrap(~name, nrow = 3)+
  xlab('Angle [degree]')+ylab('Eigenvector rank')+my_theme+scale_color_manual('Eigenvector\nrank',values = palette_div)

ggsave(p, file='/home/anne/Documents/featurecloud/pca/vertical-pca//paper/plink_comparison.pdf')

p
