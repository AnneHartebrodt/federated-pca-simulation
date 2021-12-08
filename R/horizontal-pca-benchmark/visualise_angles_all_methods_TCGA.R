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

data_all <- rbind(data1,data2, data3)
data_all<-data_all[algorithm!='vertical_pca']
data_all$algorithm<-as.factor(data_all$algorithm)
algorithms <- c('QR-PCA', 'APSTACK', 'SIT', 'APCOV', 'PCOV', 'VSIT')
levels(data_all$algorithm)<- algorithms
data_all$experiment<-ordered(data_all$experiment, levels = c('TCGA', '5', '2'))
data_all$score <- 'Angle w.r.t reference'
#data_all<-data_all[algorithm != 'Vertical power iteration']

outdir1 <- '/home/anne/Documents/featurecloud/pca/horizontal-pca/results/accuracy/pre_split//'
data1<- read_sre(outdir1)
data1$experiment <- 'TCGA'
#data1$algorithm<-as.factor(data1$algorithm)
#levels(data1$algorithm)<- c('Bai', 'Balcan', 'Power iteration', 'Proxy', 'Proxy naive', 'Proxy weighted', 'Vertical power iteration')
#data3<-data1



outdir2 <- '/home/anne/Documents/featurecloud/pca/horizontal-pca/results/accuracy/merged/5'
data2<- read_sre(outdir2)
data2$experiment <- '5'

outdir3 <- '/home/anne/Documents/featurecloud/pca/horizontal-pca/results/accuracy/merged/2'
data3<- read_sre(outdir3)
data3$experiment <- '2'

data3 <- rbind(data1,data2, data3)

data3<-data3[algorithm!='vertical_pca']
data3$algorithm<-as.factor(data3$algorithm)
levels(data3$algorithm)<- algorithms
data3$experiment<-ordered(data3$experiment, levels = c('TCGA', '5', '2'))

#data3<-data3[algorithm != 'Vertical power iteration' & dataset == 'Liver_and_intrahepatic_bile_ducts']
baseline <- unique(data3[algorithm=='SIT', .(dataset, value, name)])
colnames(baseline)<- c('dataset', 'baseline', 'name')
data3 <- data3 %>% left_join(baseline)
data3<-as.data.table(data3)
data3$value <- data3$value/data3$baseline
data3$score<-'SRE (normalised)'
data3 <- data3[, .(algorithm, dataset, name, value, experiment, score)]

data_all<-rbind(data_all, data3)

my_theme <-
  theme_classic() + theme(
    axis.title = element_text(size = 9),
    legend.title = element_text(size = 12),
    plot.title = element_text(size = 22, hjust = 0.5),
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size = 9),
    strip.text = element_text(size = 9),
    legend.text = element_text(size = 9),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    legend.key.size = unit(1, 'lines'))

comparison<- ggplot(data_all[name %in% paste0('', 1:10) &algorithm!='Vertical power iteration'], aes(name, value, fill=as.factor(algorithm)))+
  geom_boxplot(outlier.size = 0.5)+
  xlab('Eigenvalue rank')+
  scale_shape_manual(values=c(1, 8))+
  #scale_color_manual(values = c('#000000','#FF0000'))+
  guides(color='none')+
  my_theme+
  scale_fill_manual('Algorithm', values =palette_div[c(2,9,4,8,6,7,5,3)])+
  facet_grid(c('score','experiment'), scales = 'free_y', labeller = labeller(score = label_wrap_gen(width = 25)))+
  guides(fill=guide_legend(nrow=1))+
  theme(axis.title.y = element_blank(), legend.position = 'bottom')

comparison

ggsave(comparison, file='/home/anne/Documents/manuscripts/horizontal-pca/figures/comparison_angles_all_methods.pdf', width = 20, height = 10, units = 'cm')


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
algorithms <- c('QR-PCA', 'APSTACK', 'SIT', 'APCOV', 'PCOV')
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

iterations$score<-'Number of iterations'
real_com$score <- 'Transmission cost [MB]'
runstats <- rbind(iterations, real_com[, .(algorithm, cost_mb, dataset, experiment, score)], use.names=F)
runstats$algorithm<-as.factor(runstats$algorithm)
levels(runstats$algorithm)<- c('QR-PCA', 'APSTACK', 'SIT', 'APCOV', 'PCOV', 'VSIT')

my_theme <-
  theme_classic() + theme(
    axis.title = element_text(size = 24),
    legend.title = element_text(size = 24),
    plot.title = element_text(size = 24, hjust = 0.5),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 16),
    strip.text.x = element_text(size = 16),
    legend.text = element_text(size = 16),
    plot.subtitle = element_text(size = 24, hjust = 0.5),
    legend.key.size = unit(2, 'lines'))

cost<- ggplot(runstats[algorithm!='VSIT'],  aes(iterations,as.factor(algorithm), fill=experiment))+
  geom_boxplot()+
  my_theme +
  scale_fill_manual('Experiment', values = palette_div[c(2,5,8)])+
  theme(legend.position = 'bottom', axis.title.x = element_blank(), axis.title.y = element_blank())+
  scale_x_log10()+
  facet_wrap(~score)
cost

ggsave(cost, file='/home/anne/Documents/manuscripts/horizontal-pca//figures/comparison_iterations_and_cost_all_methods.pdf', width = 20, height = 10, units = 'cm')


