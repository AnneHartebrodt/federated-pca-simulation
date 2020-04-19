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
  make_option(c("-p", "--plotdir"), action="store", default=NA, type='character',
              help="The directory for plot")
)
opt = parse_args(OptionParser(option_list=option_list))
theme_set(theme_cowplot())
plotdir <-opt$p
setwd(opt$r)


myplots<-list()
metafile<-grep(x = list.files(recursive = T), pattern = 'meta_splits_perc.tsv', value = T)[1]
meta<-fread(metafile, fill = T)
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
        axis.title = element_text(size=20),
        plot.margin = unit(c(-1.0,0,0,0), "cm"),
        strip.background = element_blank(),
        strip.text.x = element_blank())+ facet_grid(~f)+
  scale_fill_grey(start=0.2, end=0.8)+
  scale_y_reverse()+ylab('%Samples')


angs <- list()
i  <- 1
for(w in c('weighted', 'balcan')){
ww<-ifelse(w == 'weighted', 'weighted', 'stacked')
for(set in list.dirs(recursive = F)){
  for (var in c('0.2','0.5')){
    if(!file.exists(file.path(set, var, 'meta_splits_perc.tsv'))){
      next
    }
    var_exp <- fread(file.path(set, var, 'nr_vars_explain_aor.tsv'))
    v = var_exp[V1==0.6]$V2
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
  ang1$method <- paste('Proxy covariance', ww)
  ang1$var<-var
  ang1<-cbind(ang1, meta[,.(i, nr_splits, sp)])
  angs[[i]]<-ang1
  i = i+1
}
}
}

ang<-rbindlist(angs, fill = T, use.names = T)
ang<-melt(ang, variable.name = 'rank', value.name = 'angle',
          id.vars = c('dataset', 'intermediated', 'method', 'i', 'nr_splits', 'sp', 'var'))
ang$angle<-as.numeric(ang$angle)
ang$rank<-as.factor(ang$rank)
}

aaa <- ang

#### Subspace
i = 1
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
    ang1$method <-'Power iteration regular'
    ang[[i]]<-ang1
    i = i+1
  }
}

for(set in list.dirs(recursive = F)){
  for(d in list.dirs(path = set, recursive = F)){
    print(d)
  if(!file.exists(file.path(d, 'meta_splits_perc.tsv'))){
    print('no meta')
  }else{
    var_exp <- fread(file.path(d, 'nr_vars_explain_aor.tsv'))
    v = var_exp[V1==0.6]$V2
    print(set)
    meta<-fread(file.path(d, 'meta_splits_perc.tsv'), fill = T)
    if(nrow(meta)!=70){
      print('flase')
      next
    }
    meta<-meta[, -9]
    colnames(meta)<-c('dataset', "i", "nr_splits", paste0('split_', 1:5))
    meta$sp<-apply(meta[,.(split_1, split_2, split_3, split_4, split_5)],1, function(x) paste(na.omit(x), collapse=' '))
    ang1<-fread(file.path(d,'power_subspace_weighted_angles_unequal_splits.tsv'), fill = T)
    ang1<-ang1[, which(unlist(apply(ang1, 2, function(x) !all(is.na(x))))), with=F]
    print(nrow(ang1))
    ang1<-cbind(ang1, meta)
    ang1$method = 'Power iteration weighted'
    ang[[i]]<-ang1
    i = i+1
  }
}
}
ang1<-rbindlist(ang)
ang1<-melt(ang1, variable.name = 'rank', value.name = 'angle', id.vars = c('V1', 'i', 'nr_splits', 'sp', 'method'))
ang1$rank<-as.factor(as.numeric(gsub("V", "", ang1$rank))-1)
ang1$angle<-as.numeric(ang1$angle)
ang1<-ang1[!is.na(ang1$angle)]
ang1<-ang1[!is.na(ang1$rank)]
colnames(ang1) <-c('dataset', 'i', 'nr_splits', 'sp', 'method', 'rank', 'angle')

all <- rbind(aaa, ang1, fill = T)
all$method <-factor(all$method, levels = c("Proxy covariance weighted", "Power iteration weighted","Proxy covariance stacked","Power iteration regular"))

col <-brewer.pal(6, 'Blues')

p<-ggplot(all[rank %in% c(1,2,3,4,5,6) & var %in% c(0.5, NA) & intermediated %in% c('2.0', NA)], aes(sp, angle, fill=rank))+
  geom_boxplot()+facet_wrap(c('method'), ncol = 2)+
  theme(plot.title = element_text(size = 35, hjust = 0.5), 
        legend.key.size = unit(0.75, "cm"),
        legend.title = element_text(size = 15),
        strip.text = element_text(size = 25), 
        axis.text = element_text(size=20),
        axis.text.x = element_blank(),
        axis.title = element_text(size=35),
        plot.subtitle = element_text(size = 30, hjust = 0.5),
        legend.text=element_text(size=20),
       legend.position = c(0.02, 0.90),
       legend.box.background = element_rect(colour = "black"),
       legend.margin = margin(1,1,1,1))+
  xlab('')+ylab('Angle [degree]')+
  scale_fill_manual('Rank', values = col)+ 
  guides(fill=guide_legend(nrow=2, byrow=TRUE))

p
pp<-plot_grid(p, annot, rel_heights = c(0.9, 0.1), align = 'v', nrow = 2, axis = 'l')
pp
ggsave(pp, file = file.path(plotdir, 'angles_all_methods.eps'), width = 20,height = 15)
ggsave(pp, file = file.path(plotdir, 'angles_all_methods.png'), width = 20,height = 15, dpi = 350)




# Supplement 0.5 and 1.0 intermediated
p<-ggplot(all[rank %in% c(1,2,3,4,5,6, 7,8,9,10) & var %in% c(0.5) & intermediated %in% c('0.5', '1.0', '2.0')], aes(sp, angle, fill=rank))+
  geom_boxplot()+facet_grid(c('intermediated','method'))+
  theme(plot.title = element_text(size = 35, hjust = 0.5), 
        legend.key.size = unit(0.75, "cm"),
        legend.title = element_text(size = 15),
        strip.text = element_text(size = 25), 
        axis.text = element_text(size=20),
        axis.text.x = element_blank(),
        axis.title = element_text(size=35),
        plot.subtitle = element_text(size = 30, hjust = 0.5),
        legend.text=element_text(size=20),
        legend.position = 'top',)+
  xlab('')+ylab('Angle [degree]')+
  scale_fill_manual('Rank', values = col)+ 
  guides(fill=guide_legend(nrow=1, byrow=TRUE))

p
pp<-plot_grid(p, annot, rel_heights = c(0.9, 0.1), align = 'v', nrow = 2, axis = 'l')
pp
ggsave(pp, file = file.path(plotdir, 'angles_proxy_all.png'), width = 20, height = 15, dpi = 350)
ggsave(pp, file = file.path(plotdir, 'angles_proxy_all.eps'), width = 20, height = 15)


