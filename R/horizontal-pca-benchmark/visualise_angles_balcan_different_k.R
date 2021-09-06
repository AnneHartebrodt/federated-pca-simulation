require(data.table)
require(ggplot2)
require(dplyr)
require(tidyr)
require("ggrepel")

source('/home/anne/Documents/featurecloud/pca/federated_dp_pca/R/horizontal-pca-benchmark/library.R')

outdir <- '/home/anne/Documents/featurecloud/pca/horizontal-pca/results/single-cell'

read_cancer_types<-function(outdir, folder){
  dirs<-list.files(outdir, no..=F)
  data_list<-list()
  for (d in dirs){
    data<-fread(file.path(outdir,folder, paste0(folder, '_angles.tsv')), sep='\t', fill=TRUE)
    
    colnames(data)<- c('algorithm', 'factor', 1:(ncol(data)-2))
    data$dataset<-d
    data<- data %>% pivot_longer(-c(algorithm,dataset, factor))
    data<-as.data.table(data)
    data_list[[d]]<-data
  }
  data <- rbindlist(data_list)
  data<-data[!is.na(value)]
  data$name<-as.factor(as.numeric(data$name))
  return(data)
}
data <- read_cancer_types(outdir = outdir, 'paga2')

gp<- ggplot(data[name %in% c('1', '2', '3', '19')], aes(name, value, fill=as.factor(factor)))+geom_boxplot(outlier.shape = NA)+
  xlab('Eigenvalue rank')+
  scale_shape_manual(values=c(1, 8))+
  #scale_color_manual(values = c('#000000','#FF0000'))+
  guides(color='none')+
  my_theme+
  #scale_fill_manual('#Factor', values = palette_div[c(1,4,8)])+
  ylab('Angle w.r.t reference [degree]')+facet_wrap(~algorithm)
  
gp
ggsave(gp, file='/home/anne/Documents/featurecloud/pca/horizontal-pca/figures/effect_of_increased_k.pdf', width = 20, height = 15, units = 'cm')

