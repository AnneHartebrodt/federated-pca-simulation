require(data.table)
require(ggplot2)
require(dplyr)
require(ggpubr)
require(cowplot)
require(gridExtra)

mm<-data.table(i = rep(0:9, 7), nr_splits=rep(c(2,2,2,5,5,5,5), each=10),
               split_1=rep(c(0.1, 0.3, 0.5, 0.2, 0.1, 0.1, 0.2375) , each=10),
               split_2=rep(c(0.9, 0.7, 0.5, 0.2, 0.1, 0.1, 0.2375), each=10), 
               split_3=rep(c(NA, NA, NA, 0.2, 0.2, 0.1, 0.2375), each=10),
               split_4=rep(c(NA,NA, NA, 0.2, 0.2, 0.1, 0.2375), each=10),
               split_5=rep(c(NA,NA,NA, 0.2, 0.4,0.6, 0.05), each=10))
mm[is.na(mm)]<-0
mm$sp<-apply(mm[,3:7],1, function(x) paste(na.omit(x), collapse=' '))
mm2<-melt(mm, id.vars = c('sp', 'i', 'nr_splits'))
setwd('/home/anne/Documents/featurecloud/results/sand/splits/')
angs = list()
i  = 1
for(w in c('weighted', 'balcan')){
for(set in list.dirs(recursive = F)){
  for (var in c('0.2', '0.5')){
    if(!file.exists(file.path(set, var, 'meta_splits.tsv'))){
      next
    }
for( d in c('0.5', '1.0', '2.0')){
  t<-fread(file.path(set, var, 'time_log.tsv'))
  meta<-fread(file.path(set, var, 'meta_splits.tsv'), fill = T)
  meta<-meta[, which(unlist(apply(meta, 2, function(x) !all(is.na(x))))), with=F]
  colnames(meta)<-c('datatset', 'i', 'nr_splits', paste('split_', 1:5, sep = ''))
  meta[,4:8]<-round(meta[,4:8]/max(na.omit(meta[, 4:8])),2)
  meta<-meta[mm, roll='nearest', on=c('i', 'nr_splits', 'split_1', 'split_2', 'split_3', 'split_4', 'split_5')]
  meta$sp<-apply(meta[,.(split_1, split_2, split_3, split_4, split_5)],1, function(x) paste(na.omit(x), collapse=' '))
  meta<-meta[-which(t[V1 %in% c('Unequal split proxy', 'Time excpetion')]=='Time excpetion')]
  ang1<-fread(file.path(set, var, paste0('prox_angles_unequal_splits_', w, '_',d,'.tsv')), fill = T)
  ang1<-ang1[, which(unlist(apply(ang1, 2, function(x) !all(is.na(x))))), with=F]
  colnames(ang1)<-c('dataset', paste('', 1:(ncol(ang1)-1), sep = ''))
  ang1$intermediated<-d
  ang1$var<-var
  ang1<-cbind(ang1, meta[,.(i, nr_splits, sp)])
  angs[[i]]<-ang1
  i = i+1
}
}
}

ang<-rbindlist(angs, fill = T)
ang<-melt(ang, variable.name = 'rank', value.name = 'angle',
          id.vars = c('dataset', 'intermediated', 'i', 'nr_splits', 'sp', 'var'))
ang$angle<-as.numeric(ang$angle)
ang$rank<-as.factor(ang$rank)
p<-ggplot(ang[rank %in% c('1', '2', '3', '4', '5', '6', '7', '8', '9')], aes(sp, angle, fill=rank))+
  geom_boxplot()+facet_grid(c('var', 'intermediated'))+
  theme(axis.text.x=element_blank(),legend.position = 'top', legend.key.size = unit(0.3, "cm"))+xlab('')+
  scale_fill_brewer(palette = 'Blues')+ guides(fill=guide_legend(nrow=1,byrow=TRUE))+
  ylab('Angle (degree)')

annot<-ggplot(mm2[value!=0], aes(sp, value, fill=variable))+
  geom_col(position = position_dodge(preserve = 'total'))+theme(legend.position = 'none')+
  theme(axis.line = element_blank(),
       axis.text = element_blank(),
       axis.ticks = element_blank(), 
       axis.title = element_blank(), 
  plot.margin= unit(c(0, 0, 0, 0), "lines"))+ 
  scale_fill_grey(start=0.2, end=0.8)+scale_y_reverse()

a<-ggarrange(annot,annot, annot, ncol=3)+theme(plot.margin= unit(c(-1.5, 1.1, 0, 3.1), "lines"))
gg<-grid.arrange(p, a, heights=c(10,1))
ggsave(file=paste0('/home/anne/Documents/paper_fed_PCA/figures/angles_eigenvectors_', w, '.pdf'), gg, dpi = 'print', width = 20, height = 15)
}

ang<-list()
for(set in list.dirs(recursive = F)){
  if(!file.exists(file.path(set, '0.5', 'meta_splits.tsv'))){
    print(set)
  }else{
  print(set)
  meta<-fread(file.path(set, '0.5', 'meta_splits.tsv'), fill = T)
  meta<-meta[, which(unlist(apply(meta, 2, function(x) !all(is.na(x))))), with=F]
  
  colnames(meta)<-c('datatset', 'i', 'nr_splits', paste('split_', 1:5, sep = ''))
  meta[,4:8]<-round(meta[,4:8]/max(na.omit(meta[, 4:8])),2)
  meta<-meta[mm, roll='nearest', on=c('i', 'nr_splits', 'split_1', 'split_2', 'split_3', 'split_4', 'split_5')]
  meta$sp<-apply(meta[,.(split_1, split_2, split_3, split_4, split_5)],1, function(x) paste(na.omit(x), collapse=' '))
  ang1<-fread(file.path(set, '0.5','dpit_subspace_angles_unequal_splits.tsv'), fill = T)
  ang1<-ang1[, which(unlist(apply(ang1, 2, function(x) !all(is.na(x))))), with=F]
  print(nrow(ang1))
  ang1<-cbind(ang1, meta)
  ang[[set]]<-ang1
  }
}
ang1<-rbindlist(ang)
ang1<-melt(ang1, variable.name = 'rank', value.name = 'angle', id.vars = c('V1', 'i', 'nr_splits', 'sp'))
ang1$angle<-as.numeric(ang1$angle)
#ang1$sp<-factor(ang1$sp, levels = c("0.5 0.5 0 0 0","0.2 0.2 0.2 0.2 0.2",  
#                                    "0.2375 0.2375 0.2375 0.2375 0.05", "0.3 0.7 0 0 0",  "0.1 0.1 0.1 0.1 0.6", "0.1 0.1 0.2 0.2 0.4","0.1 0.9 0 0 0" ))
powerp<-ggplot(ang1[rank %in% c('V2', 'V3', 'V4', 'V5', 'V6', 'V7')], aes(x= rank, y=angle, fill=sp))+geom_boxplot()+scale_fill_brewer(palette = 'Blues')+
  ggtitle('Angles between leading eigenvectors for different splits')+xlab('Rank of eigenvector')+ylab('Angle (degree)')+
  theme(legend.position = 'bottom', plot.title = element_text(size = 25), 
        legend.key.size = unit(1.5, "cm"),
        legend.title = element_text(size = 20),
        strip.text = element_text(size = 25), axis.text = element_text(size=20),
        axis.title = element_text(size=25))+guides(fill=guide_legend(nrow=1,byrow=TRUE))
powerp
ggsave(powerp, file = '/home/anne/Documents/paper_fed_PCA/figures/angles_eigenvectors_power_iteration.pdf', height = 15, width = 20)
a<-ggarrange(annot,annot, annot, ncol=3)+theme(plot.margin= unit(c(0, 1.1, 0, 3.1), "lines"))
