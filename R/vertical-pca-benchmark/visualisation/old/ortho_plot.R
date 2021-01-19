require(data.table)
require(ggplot2)
require(tidyr)
require(cowplot)
require(optparse)
require(dplyr)

palette_div<-c('#062865', '#2f497d', '#4c6d96', '#6793af', '#81bac8', '#ffc4b3', '#f68888', '#d75161', '#ab203f', '#720022')
palette_seq<-c('#062865', '#203a72', '#324d80', '#43618d', '#52759b', '#618aa9', '#70a0b7', '#7eb6c5', '#8dcdd4', '#9ce4e2')


my_theme <-
  theme_classic() + theme(
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 12),
    plot.title = element_text(size = 20, hjust = 0.5),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 12),
    plot.subtitle = element_text(size = 12, hjust = 0.5)
  )
angles<-fread('/home/anne/Documents/featurecloud/pca/vertical-pca/results/local_vs_global/mnist/164157.244188193.mnist.angles', header=F)

angles <- angles %>% pivot_longer(-c(V1, V2))
angles$name <- sapply(angles$name, function(x)as.numeric(gsub('V','', x))-2)
angles.plot<-ggplot(angles, aes(factor(name), value))+geom_boxplot()+ylab('Angle [degree] between  reference\nand naively federated eigenvectors')+my_theme+
  xlab('Eigenvector rank')
ggsave(angles.plot, filename = '/home/anne/Documents/featurecloud/pca/vertical-pca/figures/angle_naive_fed_vert.pdf')
angles.plot
