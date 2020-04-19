require(data.table)
require(ggplot2)
require(dplyr)
require(ggpubr)
require(cowplot)
require(gridExtra)
require(RColorBrewer)

require(optparse)



option_list = list(
  make_option(c("-r", "--resultfolder"), action="store", default=NA, type='character',
              help="annoation file in gtf format"),
  make_option(c("-s", "--study.name"), action="store", default=NA, type='character',
              help="annoation file in gtf format"),
  make_option(c("-t", "--sub.title"), action="store", default=NA, type='character',
              help="annoation file in gtf format"),
  make_option(c("-p", "--plotdir"), action="store", default=NA, type='character',
              help="The directory for plot")
)
opt = parse_args(OptionParser(option_list=option_list))

plotdir <-opt$p
setwd(opt$resultfolder)


study<-opt$study.name
subtitle <- opt$sub.title

files<-c('proxy_angles_unequal_splits_weighted_2.0.tsv', 
         'proxy_cluster_weighted_angles_unequal_splits.tsv',
         'power_cluster_weighted_angles_unequal_splits.tsv',
         'power_subspace_weighted_angles_unequal_splits.tsv')

run <-c('Proxy covariance\nequal split\nweighted', 
        #'Proxy unequal balcan', 
        'Proxy covariance\npreclustered\nweighted',
        'Power iteration\npreclustered\nweighted', 
        #'Power unequal', 
        'Power iteration\nequal split\nweighted')

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
  
p<-ggplot(ang[rank %in% c(1,2,3,4,5,6)  & sp %in% c('0.5 0.5', '0.2 0.2 0.2 0.2 0.2')], aes(run, angle, fill=rank))+
    geom_boxplot()+
    theme(plot.title = element_text(size = 25, hjust = 0.5), 
          legend.key.size = unit(1, "cm"),
          legend.title = element_text(size = 25),
          strip.text = element_text(size = 35), 
          axis.text = element_text(size=20),
          axis.title = element_text(size=35),
          plot.subtitle = element_text(size = 25, hjust = 0.5),
          legend.text=element_text(size=20),
          legend.position = 'top')+
    xlab('')+
    scale_fill_manual('Rank', values = col)+ 
    guides(fill=guide_legend(nrow=1, byrow=TRUE))+
    ylab('Angle [degree]')+
  ggtitle('Angles between leading eigenvectors for biased splits')+labs(subtitle = subtitle)


ggsave(file=file.path(plotdir, paste0(study, '_angles_eigenvectors.eps')), p, width = 20, height = 15)
ggsave(file=file.path(plotdir, paste0(study, '_angles_eigenvectors.png')), p, dpi = 350, width = 15, height = 10)

