require(data.table)
require(ggplot2)
require(dplyr)
require(tidyr)
require("ggrepel")
require(gridExtra)
library(ggpubr)
require(GGally)

read_iterations<-function(outdir){
  iteration_list<-list()
  dirs<-dirs<-list.files(outdir, no..=F)
  for (d in dirs){
    iterations<-fread(file.path(outdir,d, paste0(d, '_iterations.tsv')), sep='\t')
    colnames(iterations)<- c('algorithm', 'iterations')
    iterations$dataset<-d
    iteration_list[[d]]<-iterations
  }
  iterations<-rbindlist(iteration_list)
  return(iterations)
}


source('/home/anne/Documents/featurecloud/pca/federated_dp_pca/R/horizontal-pca-benchmark/library.R')

outdir1 <- '/home/anne/Documents/featurecloud/pca/horizontal-pca/results/accuracy/pre_split//'
data1<- read_cancer_types(outdir1)
data1$experiment <- 'TCGA'
data1$algorithm<-as.factor(data1$algorithm)
#levels(data1$algorithm)<- c('Bai', 'Balcan', 'Power iteration', 'Proxy', 'Proxy naive', 'Proxy weighted', 'Vertical power iteration')
#data3<-data1



outdir2 <- '/home/anne/Documents/featurecloud/pca/horizontal-pca/results/accuracy/merged/5'
data2<- read_cancer_types(outdir2)
data2$experiment <- '5'

outdir3 <- '/home/anne/Documents/featurecloud/pca/horizontal-pca/results/accuracy/merged/2'
data3<- read_cancer_types(outdir3)
data3$experiment <- '2'

data3 <- rbind(data1,data2, data3)
data3$algorithm<-as.factor(data3$algorithm)
levels(data3$algorithm)<- c('Bai', 'Balcan', 'Power iteration', 'Proxy', 'Proxy naive', 'Vertical power iteration')
data3$experiment<-ordered(data3$experiment, levels = c('TCGA', '5', '2'))

data3<-data3[algorithm != 'Vertical power iteration']

comparison<- ggplot(data3[name %in% paste0('', 1:10)], aes(name, value, fill=as.factor(algorithm)))+
  geom_boxplot(outlier.shape = NA)+
  xlab('Eigenvalue rank')+
  scale_shape_manual(values=c(1, 8))+
  scale_color_manual(values = c('#000000','#FF0000'))+
  guides(color=F)+
  my_theme+
  scale_fill_manual('Algorithm', values =palette_div[c(2,9,4,8,6,7,5,3)])+
  ylab('Angle w.r.t reference [degree]')+facet_wrap(~experiment, ncol = 2)+
  guides(fill=guide_legend(ncol=2))

c1<- ggplot(data3[name %in% paste0('', 1:10) & experiment=='TCGA'], aes(name, value, fill=as.factor(algorithm)))+
  geom_boxplot(outlier.shape = NA)+
  xlab('Eigenvalue rank')+
  scale_shape_manual(values=c(1, 8))+
  scale_color_manual(values = c('#000000','#FF0000'))+
  guides(color=F)+
  my_theme+
  facet_wrap(~experiment)+
  scale_fill_manual('Algorithm', values =palette_div[c(2,9,4,8,6,7,5,3)])+
  ylab('Angle w.r.t reference [degree]')+
  theme(legend.position = 'bottom')
c2<- ggplot(data3[name %in% paste0('', 1:10) & experiment==5], aes(name, value, fill=as.factor(algorithm)))+
  geom_boxplot(outlier.shape = NA)+
  xlab('Eigenvalue rank')+
  scale_shape_manual(values=c(1, 8))+
  scale_color_manual(values = c('#000000','#FF0000'))+
  facet_wrap(~experiment)+
  my_theme+
  scale_fill_manual('Algorithm', values =palette_div[c(2,9,4,8,6,7,5,3)])+
  ylab('Angle w.r.t reference [degree]')
c3<- ggplot(data3[name %in% paste0('', 1:10) & experiment==2], aes(name, value, fill=as.factor(algorithm)))+
    geom_boxplot(outlier.shape = NA)+
    xlab('Eigenvalue rank')+
    scale_shape_manual(values=c(1, 8))+
    scale_color_manual(values = c('#000000','#FF0000'))+
    my_theme+
  facet_wrap(~experiment)+
    scale_fill_manual('Algorithm', values =palette_div[c(2,9,4,8,6,7,5,3)])+
    ylab('Angle w.r.t reference [degree]')
comparison
l<-grab_legend(comparison)

gm<-ggmatrix(list(c1,c2,c3, l), nrow=2, ncol=2, 
             xlab = 'Eigenvector rank', ylab = 'Angle [degree] w.r.t. reference',
            showStrips = TRUE)
  

ggsave(gm, file='/home/anne/Documents/featurecloud/pca/horizontal-pca/figures/comparison_angles_all_methods.pdf', width = 20, height = 15, units = 'cm')


mn<-data1[dataset=='mnist' & name %in% paste0('', 1:10)]
av <- mn %>% select(algorithm, name, value) %>% group_by(algorithm) %>% summarise(avg_angle=mean(value))

## Also visualise the iterations until convergence.
iterations<-read_iterations(outdir1)
iterations$experiment<-'TCGA'
iterations2<-read_iterations(outdir2)
iterations2$experiment <- '5'
iterations3<-read_iterations(outdir3)
iterations3$experiment <- '2'

iterations<-rbind(iterations, iterations2, iterations3)
iterations$experiment<-ordered(iterations$experiment, levels = c('TCGA', '5', '2'))
iterations$algorithm<-as.factor(iterations$algorithm)
levels(iterations$algorithm)<-c('Bai', 'Balcan', 'Power iteration', 'Proxy', 'Proxy naive', 'Vertical power iteration')

iterations<-iterations[algorithm != 'Vertical power iteration']

gp<-ggplot(iterations, aes(algorithm, iterations, fill=experiment))+geom_boxplot()+
  my_theme+
  ylab('Number of iterations')+
  scale_fill_manual('Experiment', values = palette_div[c(2,5,8)])+
  theme(axis.text.x = element_text(angle = 90), legend.position = 'None', axis.title.x = element_blank())

gp
ggsave(gp, file='/home/anne/Documents/featurecloud/pca/horizontal-pca/figures/comparison_iterations_all_methods.pdf', width = 20, height = 15, units = 'cm')



sample_sizes<-fread('/home/anne/Documents/featurecloud/data/tcga/metadata/sample_provenance_300_10_summary.tsv')
sample_sizes[, nr_sites := .N, by=primary_site]
sample_sizes<-sample_sizes[, .(primary_site, samples, nr_sites)]
sample_sizes$a <- '&'
sample_sizes$b <- '&'
sample_sizes$c <- '\\'
sample_sizes<-unique(sample_sizes)
sample_sizes[, .(primary_site, a, samples, b, nr_sites, c)]


d <- 20269
k <- 10
algorithms <- c('Bai', 'Balcan', 'Power iteration', 'Proxy', 'Proxy naive', 'Proxy weighted', 'Vertical power iteration')
bai<-0.5*d*d+d*k
balcan_proxy<-d*k*2
powerit<-d*k*2

vert<- d*k*2
communication_factor_d <-c (bai, balcan_proxy, powerit, balcan_proxy, bai, 1, vert)
communication_factor_n <-c (0,0,0,0,0,0,k*2)

communication_factor<-data.table(algorithm=algorithms, 
                                 communication=communication_factor_d, 
                                 communication_factor_n=communication_factor_n)
real_com <- merge(iterations, communication_factor, by.x = 'algorithm', by.y = 'algorithm')
real_com<-merge(real_com, sample_sizes, by.x = 'dataset', by.y = 'primary_site')
real_com[, cost:= iterations*communication+communication_factor_n*samples]
real_com[, cost_mb:=cost*4/(1000*1000)]

cost<- ggplot(real_com,  aes(as.factor(algorithm), cost_mb, fill=experiment))+geom_boxplot()+
  ylab('Transmission cost [MB]')+
  my_theme+
  scale_fill_manual('Experiment', values = palette_div[c(2,5,8)])+
  theme(axis.text.x = element_text(angle=90), legend.position = 'bottom', axis.title.x = element_blank())
cost

combined <- ggarrange(gp, cost, common.legend = TRUE, legend="bottom")
combined
ggsave(combined, file='/home/anne/Documents/featurecloud/pca/horizontal-pca/figures/comparison_iterations_and_cost_all_methods.pdf', width = 20, height = 15, units = 'cm')


