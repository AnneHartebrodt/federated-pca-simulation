require(data.table)
require(ggplot2)
require(cowplot)

option_list = list(
  make_option(c("-r", "--resultfolder"), action="store", default=NA, type='character',
              help="annoation file in gtf format"),
  make_option(c("-p", "--plotdir"), action="store", default=NA, type='character',
              help="The directory for plot")
)
opt = parse_args(OptionParser(option_list=option_list))
theme_set(theme_cowplot())
plotdir <-opt$p
setwd(opt$r)

ang<-list()
for(set in list.dirs(recursive = F)){
  if(!file.exists(file.path(set,'meta_splits.tsv'))){
    print(set)
  }else{
    print(set)
    meta<-fread(file.path(set, 'meta_splits.tsv'), fill = T)
    meta<-meta[, which(unlist(apply(meta, 2, function(x) !all(is.na(x))))), with=F]
    
    colnames(meta)<-c('datatset', 'i', 'nr_splits', paste('split_', 1:5, sep = ''))
    meta[,4:8]<-round(meta[,4:8]/max(na.omit(meta[, 4:8])),2)
    meta<-meta[mm, roll='nearest', on=c('i', 'nr_splits', 'split_1', 'split_2', 'split_3', 'split_4', 'split_5')]
    meta$sp<-apply(meta[,.(split_1, split_2, split_3, split_4, split_5)],1, function(x) paste(na.omit(x), collapse=' '))
    ang1<-fread(file.path(set,'dpit_subspace_angles_unequal_splits.tsv'), fill = T)
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
ggsave(powerp, file = '/home/anne/Documents/paper_fed_PCA/figures/angles_eigenvectors_power_iteration_clu.pdf', height = 15, width = 20)
a<-ggarrange(annot,annot, annot, ncol=3)+theme(plot.margin= unit(c(0, 1.1, 0, 3.1), "lines"))


setwd('/home/anne/Documents/featurecloud/results/splits_clu/data_fake/')
ang<-list()
for(set in list.dirs(recursive = F)){
  if(!file.exists(file.path(set,'meta_splits.tsv'))){
    print(set)
  }else{
    print(set)
    #meta<-fread(file.path(set, 'meta_splits.tsv'), fill = T)
    #meta<-meta[, which(unlist(apply(meta, 2, function(x) !all(is.na(x))))), with=F]
    
    #colnames(meta)<-c('datatset', 'i', 'nr_splits', paste('split_', 1:5, sep = ''))
    #meta[,4:8]<-round(meta[,4:8]/max(na.omit(meta[, 4:8])),2)
    #meta<-meta[mm, roll='nearest', on=c('i', 'nr_splits', 'split_1', 'split_2', 'split_3', 'split_4', 'split_5')]
    #meta$sp<-apply(meta[,.(split_1, split_2, split_3, split_4, split_5)],1, function(x) paste(na.omit(x), collapse=' '))
    ang1<-fread(file.path(set,'power_cluster_weighted_angles_unequal_splits.tsv'), fill = T)
    ang1<-ang1[, which(unlist(apply(ang1, 2, function(x) !all(is.na(x))))), with=F]
    #print(nrow(ang1))
    #ang1<-cbind(ang1, meta)
    ang[[set]]<-ang1
  }
}
angall<-rbindlist(ang)
ang1<-melt(angall, variable.name = 'rank', value.name = 'angle', id.vars = c('V1'))
ang1$angle<-as.numeric(ang1$angle)
#ang1$sp<-factor(ang1$sp, levels = c("0.5 0.5 0 0 0","0.2 0.2 0.2 0.2 0.2",  
#                                    "0.2375 0.2375 0.2375 0.2375 0.05", "0.3 0.7 0 0 0",  "0.1 0.1 0.1 0.1 0.6", "0.1 0.1 0.2 0.2 0.4","0.1 0.9 0 0 0" ))
powerp<-ggplot(ang1[rank %in% c('V2', 'V3', 'V4', 'V5', 'V6', 'V7')], aes(x= rank, y=angle))+geom_boxplot()+scale_fill_brewer(palette = 'Blues')+
  ggtitle('Angles between leading eigenvectors for different splits')+xlab('Rank of eigenvector')+ylab('Angle (degree)')+
  theme(legend.position = 'bottom', plot.title = element_text(size = 25), 
        legend.key.size = unit(1.5, "cm"),
        legend.title = element_text(size = 20),
        strip.text = element_text(size = 25), axis.text = element_text(size=20),
        axis.title = element_text(size=25))+guides(fill=guide_legend(nrow=1,byrow=TRUE))
powerp
