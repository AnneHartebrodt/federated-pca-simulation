require(data.table)
require(ggplot2)
require(dplyr)
require(tidyr)
require("ggrepel")

source('/home/anne/Documents/featurecloud/pca/federated_dp_pca/R/horizontal-pca-benchmark/library.R')


read_and_pca<-function(type){
  data<-fread(paste0('/home/anne/Documents/featurecloud/data/tcga/cancer_type/',type,'/sites/data_all/',type,'.tsv'))
  names<-fread(paste0('/home/anne/Documents/featurecloud/data/tcga/cancer_type/',type,'/sites/data_all/',type,'.names.tsv'),
               header=F, sep='\t')
  
  names[,c('type', 'tss', 'id') := tstrsplit(V1, '-')]
  pca <- prcomp(data, scale. = F, center = F)
  
  x <- as.data.table(pca$x)[, 1:10]
  x$id<-names$id
  x$tss<-names$tss
  x$name<-names$V1
  return(x)
}

plot_save_map<-function(x, type, ol){

  ou<-unique(as.numeric(unlist(apply(x[,1:10], 2, function(a) robust.outlier(a, 6)))))
  ou.names<-x[ou]$name
  ou<-data.table(name = ou.names, method = 'mad')
  ol$method<-'angle'
  ol1<-right_join(ol, ou, by= c('algorithm'='name'))
  ol1$col<- 'both'
  ol1$col<- ifelse(!is.na(ol1$method.x) & is.na(ol1$method.y) , 'angle', 'both')
  ol1$col<- ifelse(is.na(ol1$method.x) & !is.na(ol1$method.y) , 'mad', 'both')
  x<-merge(x, ol1, by.x = 'name', by.y = 'algorithm', all.x = TRUE)
  # make this plot for the legend
  g<-ggplot(x, aes(PC1, PC2, col=as.factor(tss), shape=tss))+
    geom_point(aes(shape=tss))+
    theme_bw()+theme(legend.box = "horizontal")+
    scale_color_manual(name='TSS',values = c(palette_div, palette_div,palette_div, palette_div))+
    scale_shape_manual(name = 'TSS',values = c(1:25,1:25))+
    guides(color=guide_legend(ncol=4))+theme(legend.key.size = unit(5, 'mm'), 
                                             legend.title = element_text(size = 10))+
    geom_label_repel(data = subset(x, name %in% ol1$algorithm),
                     aes(label=id, fill=col), show.legend = FALSE, size=2)+
    scale_fill_manual("Outlier\nidentification", values = c('#fde600','#fffac5','#fcf27e' ))
  
  gd <- ggplot_build(g)
  x$color<-gd$data[[1]]$colour
  map <- x[, .(name, color)]
  
  l<- grab_legend(g)
  gp<- ggpairs(x, columns = c('PC1', 'PC2', 'PC3'), aes( shape=tss, col=tss),
               diag = NULL, 
               upper = list(continuous = "points"), 
               lower='blank')+
    scale_color_manual(name='TSS',values = c(palette_div, palette_div, palette_div, palette_div))+
    scale_shape_manual(name = 'TSS', values = c(1:25, 1:25))+
    geom_label_repel(data = subset(x, name %in% ol1$algorithm),
                     aes(label=id, fill=col), size=2,show.legend = FALSE)+
    scale_fill_manual("Outlier\nidentification", values = c('#fde600','#fffac5','#fcf27e' ))
  gp[2,2]<-l
  ggsave(gp, file = paste0('/home/anne/Documents/featurecloud/pca/horizontal-pca/figures/leave1out/', type, '_pcaplot.pdf'), 
         width=30, height=20,units = 'cm' )
  
  return(map)
}

#https://www.r-bloggers.com/detecting-outlier-samples-in-pca/
robust.outlier<-function(x, m) which( (abs(x - median(x)) / mad(x)) > m )

outdir <- '/home/anne/Documents/featurecloud/pca/horizontal-pca/results/leave1out/sample/'

data_lo <- read_cancer_types(outdir = outdir)
data_lo[,c('type', 'tss', 'id') := tstrsplit(algorithm, '-')]

leave1out<-ggplot(data_lo, aes(name, value))+geom_boxplot()+
  xlab('Eigenvalue rank')+
  guides(color=F)+
  ylab('Angle w.r.t reference [degree]')+
  my_theme
leave1out
ggsave(leave1out, file='/home/anne/Documents/featurecloud/pca/horizontal-pca/figures/effect_of_leaving_one_sample_out_all.pdf', width = 20, height = 15, units = 'cm')

# determine outlier
data_lo<- as.data.table(data_lo)
outlier_summarised<-data_lo[value>15 & name %in% paste0('', 1:10), .N, by=c('algorithm', 'dataset')]
ol<-outlier_summarised[N>1]
o<-ol[dataset == 'Corpus_uteri']


# plots per cancer type
# one boxplot and one
# pca plot
dal<-data_lo %>% group_by(dataset) %>% group_split()
d<-dal[[7]]
for(d in dal){
  d<-as.data.table(d)
  d<-d[name %in% paste0("", 1:10)]
  type <- d$dataset[1]
  pca <- read_and_pca(type)
  map <- plot_save_map(pca, type, ol)
  d<-merge(d, map, by.x = 'algorithm', by.y = 'name')
  
  leave1out<-ggplot(d, aes(name, value))+geom_boxplot()+
    xlab('Eigenvalue rank')+
    guides(color=F)+
    ylab('Angle w.r.t reference [degree]')+
    my_theme+
    geom_label_repel(data = subset(d, algorithm %in% ol$algorithm & value>15),
                     aes(label=id, col=color), size=2)+
  scale_color_identity()
  ggsave(leave1out, file=paste0('/home/anne/Documents/featurecloud/pca/horizontal-pca/figures/leave1out/l1o', d$dataset[1], '.pdf'), width = 20, height = 15, units = 'cm')
  
}

# same plot for proxy method
outdir <- '/home/anne/Documents/featurecloud/pca/horizontal-pca/results/leave1out/sample_balcan//'

data_lo <- read_cancer_types(outdir = outdir)
#outliers<-data_lo[value>20, .N, by=algorithm]
#ol<- outliers[N>2]$algorithm
leave1out<-ggplot(data_lo[name %in% paste0('', 1:10)], aes(name, value))+geom_boxplot()+
  xlab('Eigenvalue rank')+
  scale_shape_manual(values=c(1, 8))+
  scale_color_manual(values = c('#000000','#FF0000'))+
  guides(color=F)+
  scale_fill_manual('#Factor', values = palette_div[c(1,4,8)])+
  ylab('Angle w.r.t reference [degree]')+
  my_theme

#geom_label_repel(data = subset(data_lo, algorithm %in% ol & value>20), aes(label=algorithm), size=3)
leave1out
ggsave(leave1out, file='/home/anne/Documents/featurecloud/pca/horizontal-pca/figures/effect_of_leaving_one_sample_out_balcan.pdf', width = 20, height = 15, units = 'cm')


