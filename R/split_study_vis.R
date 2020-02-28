require(data.table)
require(ggplot2)

setwd('/home/anne/Documents/featurecloud/results/sandbox/splits/')
angs = list()
i  = 1
for(w in c('weighted', 'balcan')){
for(set in list.dirs(recursive = F)){
  for (var in c('0.2', '0.5')){
for( d in c('0.5', '1.0', '2.0')){
  meta<-fread(file.path(set, var, 'meta_splits.tsv'), fill = T)
  meta<-meta[, which(unlist(apply(meta, 2, function(x) !all(is.na(x))))), with=F]
  colnames(meta)<-c('datatset', 'i', 'nr_splits', paste('split_', 1:5, sep = ''))
  
  meta[,sp:=paste(split_1, split_2, split_3, split_4, split_5)]
  m<-meta
  meta<-meta[,.( i, nr_splits, sp)]
  ang1<-fread(file.path(set, var, paste0('prox_angles_unequal_splits_', w, '_',d,'.tsv')), fill = T)
  ang1<-ang1[, which(unlist(apply(ang1, 2, function(x) !all(is.na(x))))), with=F]
  colnames(ang1)<-c('dataset', paste('', 1:(ncol(ang1)-1), sep = ''))
  ang1$intermediated<-d
  ang1$var<-var
  ang1<-cbind(ang1, meta)
  angs[[i]]<-ang1
  i = i+1
}
}
}

ang<-rbindlist(angs, fill = T)
ang$split<-as.factor(ang$split)
ang<-melt(ang, variable.name = 'rank', value.name = 'angle', id.vars = c('dataset', 'intermediated', 'i', 'nr_splits', 'sp', 'var'))
ang$angle<-as.numeric(ang$angle)
m<-melt(m[,2:ncol(m), with=T], id.vars = c('i', 'nr_splits', 'sp'))

p<-ggplot(ang[rank %in% c('1', '2', '3', '4', '5', '6', '7', '8', '9')], aes(sp, angle, fill=rank))+
  geom_boxplot()+facet_grid(c('var', 'intermediated'))+
  theme(axis.text.x = element_blank(), legend.position = 'top', legend.key.size = unit(0.3, "cm"))+xlab('')+
  scale_fill_brewer(palette = 'Blues')+ guides(fill=guide_legend(nrow=1,byrow=TRUE))

annot<-ggplot(m, aes(sp, value, fill=variable))+
  geom_bar(stat = 'identity',position = 'fill' )+theme(legend.position = 'none')+
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(), 
        axis.title = element_blank(), plot.margin= unit(c(-1.5, 0, 0, 0), "lines"))+
  scale_fill_grey(start=0.2, end=0.8)
a<-ggarrange(annot,annot, annot, ncol=3)+theme(plot.margin= unit(c(0, 1.1, 0, 3.1), "lines"))
gg<-grid.arrange(p, a, heights=c(10,1))
ggsave(paste0('/home/anne/Documents/paper_fed_PCA/figures/angles_eigenvectors', w, '.pdf'), gg, dpi = 'print', width = 20, height = 15)
}

ang<-list()

for(set in list.dirs(recursive = F)){
  meta<-fread(file.path(set, '0.5', 'meta_splits.tsv'), fill = T)
  meta<-meta[, which(unlist(apply(meta, 2, function(x) !all(is.na(x))))), with=F]
  colnames(meta)<-c('datatset', 'i', 'nr_splits', paste('split_', 1:5, sep = ''))
  meta[,sp:=paste(split_1, split_2, split_3, split_4, split_5)]
  m<-meta
  meta<-meta[,.( i, nr_splits, sp)]
  ang1<-fread(file.path(set, '0.5','single_site_subspace_dpit_angles_unequal_splits.tsv' ), fill = T)
  ang1<-ang1[, which(unlist(apply(ang1, 2, function(x) !all(is.na(x))))), with=F]
  ang1<-cbind(ang1, meta)
  ang[[set]]<-ang1
}
ang1<-rbindlist(ang)
ang1<-melt(ang1, variable.name = 'rank', value.name = 'angle', id.vars = c('V1', 'i', 'nr_splits', 'sp'))

ggplot(ang1[rank %in% c('V2', 'V3', 'V4')], aes(x = rank, y=angle, fill=sp))+geom_boxplot()+facet_wrap(~V1)
