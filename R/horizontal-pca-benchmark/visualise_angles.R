require(data.table)
require(ggplot2)
require(dplyr)
require(tidyr)
require("ggrepel")

outdir <- '/home/anne/Documents/featurecloud/pca/horizontal-pca/results/pre_split'

dirs<-list.files(outdir, no..=F)

data_list<-list()


for (d in dirs){
  data<-fread(file.path(outdir,d, paste0(d, '_angles.tsv')), sep='\t')

  colnames(data)<- c('algorithm', 1:(ncol(data)-1))
  data$dataset<-d
  data<- data %>% pivot_longer(-c(algorithm,dataset))
  data<-as.data.table(data)
  data_list[[d]]<-data
}

data <- rbindlist(data_list)
data<-data[!is.na(value)]
data$name<-as.factor(as.numeric(data$name))
data[, mnist:=as.factor(ifelse(dataset=='mnist', 2, 1))]
ggplot(data, aes(name, value, fill=algorithm, col=mnist))+geom_boxplot(outlier.shape = NA)+xlab('Eigenvalue rank')+
  scale_shape_manual(values=c(1, 8))+scale_color_manual(values = c('#000000','#FF0000'))+guides(color=F)

iteration_list<-list()
for (d in dirs){
  iterations<-fread(file.path(outdir,d, paste0(d, '_iterations.tsv')), sep='\t')
  colnames(iterations)<- c('algorithm', 'iterations')
  iterations$dataset<-d
  iteration_list[[d]]<-iterations
}
iterations<-rbindlist(iteration_list)
ggplot(iterations, aes(algorithm, iterations, label=dataset))+geom_boxplot(outlier.shape = NA)+geom_jitter()+
  geom_label_repel(data = subset(iterations, iterations>600 | dataset == 'mnist'),
                   aes(label=dataset), size=3)
