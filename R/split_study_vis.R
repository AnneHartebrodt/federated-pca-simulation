require(data.table)
require(ggplot2)
require(dplyr)
require(ggpubr)
require(cowplot)
require(gridExtra)


setwd('/home/anne/Documents/featurecloud/results_final/split/')
meta<-fread('TCGA-COAD/0.5/meta_splits_perc.tsv', fill = T)
meta<-meta[, -9]
colnames(meta)<-c('dataset', "i", "nr_splits", paste0('split_', 1:5))
meta$sp<-apply(meta[,.(split_1, split_2, split_3, split_4, split_5)],1, function(x) paste(na.omit(x), collapse=' '))

res<-list()
mm<-melt(meta, c('dataset', 'i', 'nr_splits', 'sp'))
mm<-mm[!is.na(value)]
mm2<-mm
mm$f<-1
mm2$f<-2
mm<-rbind(mm, mm2)

annot<-ggplot(mm, aes(sp, value))+
  geom_col(aes(fill=variable),position = position_dodge(preserve = 'total'))+
  theme(legend.position = 'none',
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.text = element_text(size=20),
        axis.title = element_text(size=20, hjust = 0, vjust = 0.5),
        plot.margin = unit(c(-1.0,0,0,0), "cm"),
        strip.background = element_blank(),
        strip.text.x = element_blank())+ facet_grid(~f)+
  scale_fill_grey(start=0.2, end=0.8)+
  scale_y_reverse()+ylab('%Samples')



for(w in c('weighted', 'balcan')){
  angs = list()
  i  = 1
for(set in list.dirs(recursive = F)){
  print (set)
  var_exp <- fread(file.path(set, var, 'nr_vars_explain_aor.tsv'))
  v = var_exp[V1==0.6]$V2
  for (var in c('0.2','0.5')){
    if(!file.exists(file.path(set, var, 'meta_splits_perc.tsv'))){
      next
    }
for( d in c('0.5', '1.0', '2.0')){
  meta<-fread(file.path(set, var, 'meta_splits_perc.tsv'), fill = T)
  meta<-meta[, -9]
  colnames(meta)<-c('dataset', "i", "nr_splits", paste0('split_', 1:5))
  meta$sp<-apply(meta[,.(split_1, split_2, split_3, split_4, split_5)],1, function(x) paste(na.omit(x), collapse=' '))
  ang1<-fread(file.path(set, var, paste0('proxy_angles_unequal_splits_', w, '_',d,'.tsv')), fill = T)
  ang1<-ang1[, which(unlist(apply(ang1, 2, function(x) !all(is.na(x))))), with=F]
  ang1<-ang1[, 1:(v+1), with=F]
  colnames(ang1)<-c('dataset', paste('', 1:(ncol(ang1)-1), sep = ''))
  ang1$intermediated<-d
  ang1$var<-var
  ang1<-cbind(ang1, meta[,.(i, nr_splits, sp)])
  angs[[i]]<-ang1
  i = i+1
}
}
}

ang<-rbindlist(angs, fill = T, use.names = T)
ang<-melt(ang, variable.name = 'rank', value.name = 'angle',
          id.vars = c('dataset', 'intermediated', 'i', 'nr_splits', 'sp', 'var'))
ang$angle<-as.numeric(ang$angle)
ang$rank<-as.factor(ang$rank)

res[[w]]<-ang
col<-c(brewer.pal(name='Blues', n = 9),brewer.pal(name='Blues', n = 3))

subti <- ifelse(w=='weighted', 'Weighted Proxy Covariance Matrix', 'Stacked Principal Components')

p<-ggplot(ang[var ==0.5 & intermediated %in% c('1.0', '2.0')], aes(sp, angle, fill=rank))+
  geom_boxplot()+facet_wrap(c('intermediated'))+
   theme(plot.title = element_text(size = 25, hjust = 0.5), 
               legend.key.size = unit(1, "cm"),
               legend.title = element_text(size = 25),
               strip.text = element_text(size = 35), 
               axis.text = element_text(size=20),
         axis.text.x = element_blank(),
               axis.title = element_text(size=35),
               plot.subtitle = element_text(size = 25, hjust = 0.5),
               legend.text=element_text(size=20),
         legend.position = 'top')+
  xlab('')+
  scale_fill_manual('Rank', values = col)+ guides(fill=guide_legend(nrow=1, byrow=TRUE))+
  ylab('Angle [degree]')+ggtitle('Angles between leading eigenvectors for different split')+
  labs(subtitle =subti )



gg<-grid.arrange(p, annot, heights=c(10,1))
gg

ggsave(file=paste0('/home/anne/Documents/paper_fed_PCA/figures/angles_eigenvectors_', w, '.pdf'), gg, dpi = 'print', width = 20, height = 15)
}





ang<-list()
for(set in list.dirs(recursive = F)){
  if(!file.exists(file.path(set, '0.5', 'meta_splits_perc.tsv'))){
    print(set)
  }else{
    var_exp <- fread(file.path(set, var, 'nr_vars_explain_aor.tsv'))
    v = var_exp[V1==0.6]$V2
  print(set)
    meta<-fread(file.path(set, var, 'meta_splits_perc.tsv'), fill = T)
    meta<-meta[, -9]
    colnames(meta)<-c('dataset', "i", "nr_splits", paste0('split_', 1:5))
  meta$sp<-apply(meta[,.(split_1, split_2, split_3, split_4, split_5)],1, function(x) paste(na.omit(x), collapse=' '))
  ang1<-fread(file.path(set, '0.5','power_subspace_angles_unequal_splits.tsv'), fill = T)
  ang1<-ang1[, which(unlist(apply(ang1, 2, function(x) !all(is.na(x))))), with=F]
  #ang1<-ang1[, 1:(v+1), with=F]
  print(nrow(ang1))
  ang1<-cbind(ang1, meta)
  ang[[set]]<-ang1
  }
}
ang1<-rbindlist(ang)
ang1<-melt(ang1, variable.name = 'rank', value.name = 'angle', id.vars = c('V1', 'i', 'nr_splits', 'sp'))
ang1$rank<-as.factor(as.numeric(gsub("V", "", ang1$rank))-1)
ang1$angle<-as.numeric(ang1$angle)
ang1<-ang1[!is.na(ang1$angle)]
ang1<-ang1[!is.na(ang1$rank)]
col<-c(brewer.pal(name='Blues', n = 9),brewer.pal(name='Blues', n = 9),brewer.pal(name='Blues', n = 3))
powerp<-ggplot(ang1, aes(x= sp, y=angle, fill=rank))+
  geom_boxplot()+scale_fill_manual("Rank", values = col)+ 
  ggtitle('Angles between leading eigenvectors for different splits (Subspace Iteration)')+
  xlab('')+ylab('Angle (degree)')+
  theme(plot.title = element_text(size = 25, hjust = 0.5), 
        legend.key.size = unit(1, "cm"),
        legend.title = element_text(size = 25),
        strip.text = element_text(size = 35), 
        axis.text = element_text(size=20),
        axis.text.x = element_blank(),
        axis.title = element_text(size=35),
        plot.subtitle = element_text(size = 25, hjust = 0.5),
        legend.text=element_text(size=20),
        legend.position = 'top') +
  guides(fill=guide_legend(nrow=1, byrow=TRUE))
powerp

annot<-ggplot(mm, aes(sp, value))+
  geom_col(aes(fill=variable),position = position_dodge(preserve = 'total'))+
  theme(legend.position = 'none',
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.text = element_text(size=20),
        axis.title = element_text(size=20, hjust = 0, vjust = 0.5),
        plot.margin = unit(c(-1.0,0,0,0), "cm"),
        strip.background = element_blank(),
        strip.text.x = element_blank())+
  scale_fill_grey(start=0.2, end=0.8)+
  scale_y_reverse()+ylab('%Samples')


gp<-grid.arrange(powerp, annot, heights=c(10,1))
gp
ggsave(gp, file = '/home/anne/Documents/paper_fed_PCA/figures/angles_eigenvectors_power_iteration.pdf', height = 15, width = 20)


