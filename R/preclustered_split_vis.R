require(data.table)
require(ggplot2)
require(dplyr)
require(ggpubr)
require(cowplot)
require(gridExtra)
require(RColorBrewer)

plotdir <-'/home/anne/Documents/featurecloud/results_final/plots/'

#setwd('/home/anne/Documents/featurecloud/results_final/splits_fake/')
#study<-'fake'

setwd('/home/anne/Documents/featurecloud/results_final/splits_preclustered/')
study<-'preclustered'


files<-c('proxy_angles_unequal_splits_weighted_1.0.tsv', 
         'proxy_angles_unequal_splits_balcan_1.0.tsv',
         'proxy_cluster_weighted_angles_unequal_splits.tsv',
         'power_cluster_weighted_angles_unequal_splits.tsv',
         'power_subspace_angles_unequal_splits.tsv',
         'power_subspace_weighted_angles_unequal_splits.tsv')

run <-c('Proxy unequal weighted', 
        'Proxy unequal balcan', 
        'Proxy cluster weighted',
        'Power cluster weighted', 
        'Power unequal', 
        'Power unequal weighted')

angs<-list()
i = 1
for(set in list.dirs(recursive = F)){
  for(fi in 1:length(files)){
    for(d in list.dirs(path = set, recursive = F)){
      if(!file.exists(file.path(d, 'meta_splits_perc.tsv'))){
        next
      }
      if(!file.exists(file.path(d, files[fi]))){
        next
      }
      print(set)
      print(files[fi])
      meta<-fread(file.path(d, 'meta_splits_perc.tsv'), fill = T)
      
      meta<-meta[, -9]
      colnames(meta)<-c('dataset', "i", "nr_splits", paste0('split_', 1:5))
      meta$sp<-apply(meta[,.(split_1, split_2, split_3, split_4, split_5)],1, function(x) paste(na.omit(x), collapse=' '))
      ang1<-fread(file.path(d, files[fi]), fill = T)
      ang1<-ang1[, which(unlist(apply(ang1, 2, function(x) !all(is.na(x))))), with=F]
      if(nrow(ang1)<7 & nrow(ang1)>1){
        next
      }
      colnames(ang1)<-c('dataset', paste('', 1:(ncol(ang1)-1), sep = ''))
      
      ang1<-cbind(ang1, meta[,.(i, nr_splits, sp)])
      ang1$run<-run[fi]
      angs[[i]]<-ang1
      i = i+1
    }
  }
 } 
ang<-rbindlist(angs, fill = T, use.names = T)
ang<-melt(ang, variable.name = 'rank', value.name = 'angle',
            id.vars = c('dataset', 'i', 'nr_splits', 'sp', 'run'))
ang$angle<-as.numeric(ang$angle)
ang$rank<-as.factor(ang$rank)
  
col<-c(brewer.pal(name='Blues', n = 9),brewer.pal(name='Blues', n = 9),brewer.pal(name='Blues', n = 9))
  
#subti <- ifelse(w =='weighted', 'Weighted Proxy Covariance Matrix', 'Stacked Principal Components')
  
p<-ggplot(ang[rank %in% c(1,2,3,4,5,6,7,8,9,10)  & sp == '0.5 0.5'], aes(run, angle, fill=rank))+
    geom_boxplot()+
    theme(plot.title = element_text(size = 25, hjust = 0.5), 
          legend.key.size = unit(1, "cm"),
          legend.title = element_text(size = 25),
          strip.text = element_text(size = 35), 
          axis.text = element_text(size=20),
          axis.text.x = element_text(angle = 90),
          axis.title = element_text(size=35),
          plot.subtitle = element_text(size = 25, hjust = 0.5),
          legend.text=element_text(size=20),
          legend.position = 'top')+
    xlab('')+
    scale_fill_manual('Rank', values = col)+ 
    guides(fill=guide_legend(nrow=1, byrow=TRUE))+
    ylab('Angle [degree]')+
  ggtitle('Angles between leading eigenvectors for different split')


ggsave(file=file.path(plotdir, paste0(study, '_angles_eigenvectors.pdf')), p, dpi = 'print', width = 20, height = 15)


