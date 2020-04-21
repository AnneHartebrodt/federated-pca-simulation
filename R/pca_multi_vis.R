require(data.table)
require(ggplot2)
require(gridExtra)
require(RColorBrewer)
require(cowplot)
library(ggplot2)
library(grid)
library(gtable)
require(optparse)


option_list = list(
  make_option(c("-r", "--resultfolder"), action="store", default=NA, type='character',
              help="annoation file in gtf format"),
  make_option(c("-p", "--plotdir"), action="store", default=NA, type='character',
              help="The directory for plot")
)
opt = parse_args(OptionParser(option_list=option_list))

plotdir <-opt$plotdir
setwd(opt$resultfolder)
theme_set(theme_cowplot())


eigenvalue<-function(A, v){
  Av<-A %*% v
  return(v %*% Av)
}

angle <- function(x,y){
  dot.prod <- x%*%y 
  norm.x <- norm(x,type="2")
  norm.y <- norm(y,type="2")
  t<-dot.prod / (norm.x * norm.y)
  #print(t)
  #Make the angle save to numerical inaccuracies
  if(t<=(-1)){
    #print('t<1')
    t<--1
  }
  if(t>1){
    #print('t>1')
    t<-1
  }
  theta <- acos(t)
  #print(theta)
  #print(as.numeric(theta))
  a<-as.numeric(theta)*180/pi
  #print(a)
  if(is.nan(a)){
    #print(a)
    return(1)
  }
  #print(a)
  a
}

canonize<-function(t, a){
if(angle(t,a)>=90){
  return(-t)
}else{
    return(t)
  }
}


if(!dir.exists(file.path(plotdir,'example_vis'))){
  dir.create(file.path(plotdir,'example_vis'))
}

setwd(opt$resultfolder)
#import splits
ma<-data.table(i = rep(0:9, 7), nr_splits=rep(c(2,2,2,5,5,5,5), each=10),
               perc_1=rep(c(0.1, 0.3, 0.5, 0.2, 0.1, 0.1, 0.2375) , each=10),
               perc_2=rep(c(0.9, 0.7, 0.5, 0.2, 0.1, 0.1, 0.2375), each=10), 
               perc_3=rep(c(NA, NA, NA, 0.2, 0.2, 0.1, 0.2375), each=10),
               perc_4=rep(c(NA,NA, NA, 0.2, 0.2, 0.1, 0.2375), each=10),
               perc_5=rep(c(NA,NA,NA, 0.2, 0.4,0.6, 0.05), each=10))
ma[is.na(ma)]<-0
meta<-fread('meta_splits.tsv', fill = T)
meta<-meta[, 1:8]
meta<-cbind(meta, round(meta[,4:8]/max(na.omit(meta[, 4:8])),2))

colnames(meta)<-c('dataset', 'i', 'nr_splits', paste0(1:5, '_split'), paste0('perc_', 1:5))
ma<-meta[ma, roll='nearest', on=c('i', 'nr_splits', 'perc_1', 'perc_2', 'perc_3', 'perc_4', 'perc_5')]
ma<-unique(ma[ ,c('perc_1', 'perc_2', 'perc_3', 'perc_4', 'perc_5')])



import_splits<-function(file_id_ev, file_id_data, meta){
  mm<-apply(meta, 1, function(x) paste0(trimws(na.omit(x[c(4,5,6,7,8)])), collapse = '.0_'))
eigenvectors<-list()
for(m in mm){
  eigenvector<-fread(file.path(paste0(file_id_ev, '_',m,'.0.tsv')))
  eigenvectors[[m]]<-eigenvector[2:nrow(eigenvector), 2:ncol(eigenvector)]
}


stops<-apply(meta[, c('1_split','2_split','3_split','4_split','5_split')], 1, function(x) as.integer(x[!is.na(x)]))
start=0
sd_list<-list()
i=1
for(l in 1:length(stops)){
  start=0
  sub_data<-list()
  for(end in stops[[l]]){
    name<-paste0(file_id_data, '_',start,'_',end, '.tsv')
    fi<-fread(name, header=T)
    fi<-fi[, 2:ncol(fi)]
    sub_data[[name]]<-fi
    start<-end
  }
  sd_list[[mm[[l]]]]<-sub_data
  i=i+1
}
return(list(sd_list, eigenvectors))
}


add_extra_legend<-function(myplot, label){
  # Get the legend
  g = ggplotGrob(myplot)
  leg = g$grobs[[which(g$layout$name == "guide-box")]]
  
  # Construct the label grob 
  xpos = 5
  textgrob = textGrob(x = unit(xpos, "points"), label, gp = gpar(fontsize =25, cex = .75), just = "left")
  width = unit(1, "grobwidth",textgrob) + unit(2*xpos, "points")  # twice the x position
  height = unit(1, "grobheight", textgrob)+ unit(2*xpos, "points")
  rectgrob = rectGrob(x = unit(0, "points"), just = "left", 
                      gp = gpar(colour = "black", fill = NA), height = height, width = width)
  labelGrob = gTree("labelGrob", children = gList(rectgrob, textgrob))
  
  # Add the label grob to a new row added to the legend
  pos = subset(leg$layout, grepl("guides", name), t:r)
  
  leg = gtable_add_rows(leg, height, pos = pos$t+1)
  leg = gtable_add_grob(leg, labelGrob, t = pos$t+2, l = pos$l)
  
  # Adjust the middle width of the legend to be the maximum of the original width 
  # or the width of the grob
  leg$widths[pos$l] = max(width, leg$widths[pos$l])
  
  # Add some space between the two parts of the legend
  leg$heights[pos$t+1] = unit(5, "pt")
  
  # Return the modified legend to the origial plot
  g$grobs[[which(g$layout$name == "guide-box")]] = leg
  
  # Adjust the width of the column containing the legend to be the maximum 
  # of the original width or the width of the label
  g$widths[g$layout[grepl("guide-box", g$layout$name), "l"]] = max(width, sum(leg$widths))
  
  # Draw the plot
  grid.newpage()
  #grid.draw(g)
  g
}


compare_local_proj_vs_global_proj<-function(eigenvector, data.complete, clu, sublist, das, subtitle, extra){
  PCs<-as.matrix(data.complete) %*% as.matrix(eigenvector)
  
  PCs<-as.data.table(PCs)
  PCs$dataset<-'projected.complete'
  PCs$color <-'p.complete'
  pca<-prcomp(data.complete, center = F)
  print('PCA COMPLETE')
  ev<-sapply(1:(ncol(PCs)-2), function(x) canonize(pca$rotation[,x], as.matrix(eigenvector)[,x]))
  pp<-as.matrix(data.complete) %*% as.matrix(ev)
  pp<-as.data.table(pp)
  pp$dataset<-'pca.complete'
  pp$color<-'global'
  print(colnames(PCs))
  colnames(pp)
  
  PCs<-rbind(pp,PCs, use.names=F)
  
  
  c0 <-ncol(PCs)-2
  i=1
  for(sub in sublist){
    PCL<-as.matrix(sub) %*% as.matrix(eigenvector)
    PCL<-as.data.table(PCL)
    PCL$dataset<-paste0('sub', i)
    colnames(PCL)<-c(paste0('V',as.numeric(gsub("V", "", colnames(PCL)[1:(ncol(PCL)-1)]))-1), 'dataset')
    print(paste0('sub', i))
    PCL$color<-'pooled'
    
    sub_local<-prcomp(sub, center = F)
    lim <- min(c0, ncol(sub_local$rotation))
    print(lim)
    ev<-sapply(1:lim, function(x) canonize(sub_local$rotation[,x], as.matrix(eigenvector)[,x]))
    pp<-as.matrix(sub) %*% as.matrix(ev)
    pp<-as.data.table(pp)
    pp$dataset<-paste0('sub.local', i)
    pp$color<-'local'
    i=i+1
    PCs<-rbind(PCs, PCL, pp, use.names=T, fill=T)
  }
  
 
  print(subtitle)
  PCs$dataset<-as.factor(PCs$dataset)
  colnames(PCs)<-gsub("V", "PC", colnames(PCL))
  ggs1<-ggplot(PCs[color %in% c('global', 'local', 'pooled')], aes(PC1,PC2,colour=as.factor(color)))+geom_point()+
    ggtitle('PC1 vs. PC2 (global, pooled and local) projections')+
    scale_color_manual("Model",values = brewer.pal(name = 'RdBu', n=11)[c(2,3,11)])+
    labs(subtitle = subtitle)+
    theme(plot.title = element_text(size = 35, hjust = 0.5), 
                                    legend.key.size = unit(1, "cm"),
                                    legend.title = element_text(size = 25),
                                    strip.text = element_text(size = 35), 
                                    axis.text = element_text(size=20),
                                    axis.title = element_text(size=35),
                                    plot.subtitle = element_text(size = 25, hjust = 0.5),
                                    legend.text=element_text(size=20), legend.position = 'bottom')+
    guides(color = guide_legend(override.aes = list(size = 5)))
  #add_extra_legend(ggs1, extra)
}







#### data improt
complete<-fread('data_new.tsv')

### Proxy covariance weighted
file_id_data<-'proxy'
file_id_ev<-'proxy_eigenvectors_weighted_2.0'
angles <-fread('proxy_angles_unequal_splits_weighted_2.0.tsv')
spl<-import_splits(file_id_ev = file_id_ev, file_id_data = file_id_data, meta)
sd_list<-spl[[1]]
eigenvector <-spl[[2]]
ps<-list()
i=1
k = 1
j = 2
for(l in 1:length(sd_list)){
#for(l in 1:1){
tf <- ma[l]!=0
extra = paste0('Angles to reference\neigenvectors:\nPC',k,': ', round(angles[l,k+1, with=F]), '\nPC', j, ': ',round(angles[l,j+1,with=F]))

subtitle = paste('Proxy covariance (weighted), splits:', paste0(ma[l, tf, with=F], collapse = ' '))
ps[[i]]<-compare_local_proj_vs_global_proj(eigenvector[[l]], complete, clu, sublist = sd_list[[l]], 
                                           names(sd_list)[l], subtitle = subtitle, extra = extra)
na<-file.path(plotdir, paste0('proxy_weighted_',l,'.eps'))
ggsave(ps[[i]], filename = na, width = 15, height = 10, dpi = 'print')
i=i+1
}
plot(ps[[1]])


#### proxy balcan
file_id_data<-'proxy'
file_id_ev<-'proxy_eigenvectors_balcan_2.0'
angles <-fread('proxy_angles_unequal_splits_balcan_2.0.tsv')
spl<-import_splits(file_id_ev = file_id_ev, file_id_data = file_id_data, meta)
sd_list<-spl[[1]]
eigenvector <-spl[[2]]
ps<-list()
i=1
k = 1
j = 2
for(l in 1:length(sd_list)){
  tf <- ma[l]!=0
  extra = paste0('Angles to reference\neigenvectors:\nPC',k,': ', round(angles[l,k+1, with=F]), '\nPC', j, ': ',round(angles[l,j+1,with=F]))
  subtitle = paste('Single round (stacked), splits:', paste0(ma[l, tf, with=F], collapse = ' '))
  na<-file.path(plotdir, paste0('proxy_balcan_',l,'.eps'))
  ps[[i]]<-compare_local_proj_vs_global_proj(eigenvector[[l]], complete, clu, sublist = sd_list[[l]], 
                                             names(sd_list)[l], subtitle = subtitle, extra=extra)
  ggsave(ps[[i]], filename = na, width = 15, height = 10, dpi = 350)
  i=i+1
}


### proxy data sub
file_id_data<-'power_datasub'
file_id_ev<-'power_subspace_eigenvectors'
angles <-fread('power_subspace_angles_unequal_splits.tsv')

spl<-import_splits(file_id_ev = file_id_ev, file_id_data = file_id_data, meta)
sd_list<-spl[[1]]
eigenvector <-spl[[2]]
ps<-list()
i=1
k = 1
j = 2
for(l in 1:length(sd_list)){
  #for(l in 1:1){
  tf <- ma[l]!=0
  subtitle = paste('Power iteration, splits:', paste0(ma[l, tf, with=F], collapse = ' '))
  extra = paste0('Angles to reference\neigenvectors:\nPC',k,': ', round(angles[l,k+1, with=F]), '\nPC', j, ': ',round(angles[l,j+1,with=F]))
  
  ps[[i]]<-compare_local_proj_vs_global_proj(eigenvector[[l]], complete, clu, 
                                             sublist = sd_list[[l]], names(sd_list)[l],
                                             subtitle = subtitle, extra = extra)
  na<-file.path(plotdir, paste0('subspace_iteration_',l,'.eps'))
  ggsave(ps[[i]], filename = na, width = 15, height = 10, dpi = 'print')
  i=i+1
}

### weigthed power iteration
file_id_data<-'power_datasub'
file_id_ev<-'power_subspace_weighted_eigenvectors'
angles <-fread('power_subspace_weighted_angles_unequal_splits.tsv')

spl<-import_splits(file_id_ev = file_id_ev, file_id_data = file_id_data, meta)
sd_list<-spl[[1]]
eigenvector <-spl[[2]]
ps<-list()
i=1
k = 1
j = 2
for(l in 1:length(sd_list)){
  #for(l in 1:1){
  tf <- ma[l]!=0
  subtitle = paste('Power iteration, splits:', paste0(ma[l, tf, with=F], collapse = ' '))
  extra = paste0('Angles to reference\neigenvectors:\nPC',k,': ', round(angles[l,k+1, with=F]), '\nPC', j, ': ',round(angles[l,j+1,with=F]))
  
  ps[[i]]<-compare_local_proj_vs_global_proj(eigenvector[[l]], complete, clu, 
                                             sublist = sd_list[[l]], names(sd_list)[l],
                                             subtitle = subtitle, extra = extra)
  na<-file.path(plotdir, paste0('subspace_iteration_weighted_',l,'.eps'))
  ggsave(ps[[i]], filename = na, width = 15, height = 10, dpi = 'print')
  i=i+1
}

