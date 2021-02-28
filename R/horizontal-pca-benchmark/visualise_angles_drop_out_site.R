require(data.table)
require(ggplot2)
require(dplyr)
require(tidyr)
require("ggrepel")

source('/home/anne/Documents/featurecloud/pca/federated_dp_pca/R/horizontal-pca-benchmark/library.R')

outdir <- '/home/anne/Documents/featurecloud/pca/horizontal-pca/results/leave1out/random/'
data_random <- read_cancer_types(outdir)
data_random$experiment <-'random'
outdir <- '/home/anne/Documents/featurecloud/pca/horizontal-pca/results/leave1out/site'
data_site <- read_cancer_types(outdir = outdir)
data_site$experiment <- 'site'

data<-rbind(data_random, data_site)
leave1out<-ggplot(data, aes(name, value, fill=experiment))+geom_boxplot()+
  xlab('Eigenvalue rank')+
  scale_shape_manual(values=c(1, 8))+
  scale_color_manual(values = c('#000000','#FF0000'))+
  guides(color=F)+
  scale_fill_manual('Experiment', values = palette_div[c(1,4,8)])+
  ylab('Angle w.r.t reference [degree]')+
  my_theme
#+
  geom_label_repel(data = subset(data_lo, algorithm %in% ol & value>20), aes(label=algorithm), size=3)

leave1out

ggsave(leave1out, file='/home/anne/Documents/featurecloud/pca/horizontal-pca/figures/effect_of_leaving_one_site_out.pdf', width = 20, height = 15, units = 'cm')

