require(data.table)
require(ggplot2)
require(tidyr)
require(cowplot)
require(optparse)
require(dplyr)
require(R.utils)
library(ggthemes)   
require(scales)
require(R.utils)
library(ggforce)
require(stringr)

palette_div<-c('#062865', '#2f497d', '#4c6d96', '#6793af', '#81bac8', '#ffc4b3', '#f68888', '#d75161', '#ab203f', '#720022')
palette_seq<-c('#062865', '#203a72', '#324d80', '#43618d', '#52759b', '#618aa9', '#70a0b7', '#7eb6c5', '#8dcdd4', '#9ce4e2')

my_theme <-
  theme_classic() + theme(
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 12),
    plot.title = element_text(size = 20, hjust = 0.5),
    axis.text = element_text(size = 8),
    legend.text = element_text(size = 8),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    legend.key.size = unit(1.5, 'lines'))

read_files_from_dir <- function(directory, prefix, suffix, splits) {
  myfiles <-  list.files(directory, recursive = F)
  myangles <- list()
  mylastlines <- list()
  for (i in splits) {
    for (f in which(str_detect(myfiles, paste0("^", prefix, i, ".*.", suffix, '$')))) {
      filename <- myfiles[f]
      print(file.path(directory, filename))
      angles <-
        fread(
          file.path(directory, filename),
          sep = '\t',
          header = F,
          fill=T
        )
      if (is.null(angles)) {
        next
      }
      colnames(angles) <-c('iterations', sapply(1:(ncol(angles) - 1), function(x) paste0('', x)))
      angles$filename <- filename
      angles$dataset<-prefix
      angles$splits <- i
      angles <-
        as.data.table(pivot_longer(angles,-c(iterations, dataset, splits, filename)))
      myangles[[filename]]<-angles
    }
  }
  
  myangles <- rbindlist(myangles)
  myangles <- myangles[!is.na(value)]
  colnames(myangles)<-c('it','filename', 'dataset', 'splits', 'rank','angle')
  
  return (myangles)
}

read_files_from_dir_indiv <- function(directory, prefix, suffix, splits) {
  myfiles <-  list.files(directory)
  myangles <- list()
  mylastlines <- list()
  for (i in splits) {
    for (f in which(str_detect(myfiles, paste0("^", prefix, i, ".*.", suffix, '$')))) {
      filename <- myfiles[f]
      print(file.path(directory, filename))
      angles <-
        fread(
          file.path(directory, filename),
          sep = '\t',
          header = F,
          fill=T
        )
      if (is.null(angles)) {
        print('not found')
        next
      }
      angles<-angles[,1:3]
      colnames(angles) <- c('rank', 'it', 'angle')
      angles$filename <- filename
      angles$dataset <- prefix
      angles$splits <- i
      myangles[[filename]] <- angles
      
      
    }
  }
  if (length(myangles)!=0){
    myangles <- rbindlist(myangles)
    colnames(myangles)<-c('rank', 'it', 'angle', 'filename','dataset','splits')
    return (myangles)
  }
  else{
    return(NULL)
  }
}




make_plot_correlation<-function(sumstat, evrank, lim, splitl){
  angles.plot<-ggplot(sumstat[splits %in% splitl & rank==evrank], aes(it, mean_angle, col=dataset))+
    geom_line(aes(linetype=dataset) ,size=1)+
    facet_wrap(~splits, scales = 'free')+
    xlim(lim)+
    my_theme+ylab('Pearson correlation')+ 
    xlab('# iterations')+
    scale_color_manual('Algorithm', values = palette_div[c(1,3,8)])+
    theme(axis.line=element_line(),
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          legend.position = c(1, 0.75),
          legend.justification = c("right", "top"),
          legend.box.just = "right",
          legend.margin = margin(6, 6, 6, 6))+
    guides(linetype=FALSE, color=guide_legend(keyheight = 0.75, title = element_text('Algorithm', size = 10)))
  return(angles.plot)
}




combine_data<-function(directory, prefix, suffix, qr, splits){
  # read the files
  angles_qr <- read_files_from_dir(directory, paste0(prefix,  '_', qr) , suffix, splits)
  angles_guo <- read_files_from_dir_indiv(directory, paste0(prefix,  '_guo_', qr), suffix, splits)
 
  angles_guo$dataset <- 'Vector'
  angles_qr$dataset<-'Matrix'
  
  # combine the data sets
  a <- rbind(angles_guo, angles_qr)
  a$angle<-abs(a$angle)
  sumstat <- a %>% group_by(rank, it, dataset, splits) %>% summarise(mean_angle=mean(angle))
  sumstat <- as.data.table(sumstat)
  sumstat$rank<-as.numeric(sumstat$rank)
  
  return(sumstat)
}

create_summary<-function(directory, prefix, suffix, splits){
  # read gradient files
  directory_grad<-file.path(directory, 'gradient')
  sumstat<-combine_data(directory_grad, prefix,  suffix, 'central_qr_', splits)
  sumstat$qr <- 'central'
  sumstat2<-combine_data(directory_grad, prefix,  suffix, 'fed_qr_', splits )
  sumstat2$qr <- 'federated'
  sumstat<-rbind(sumstat, sumstat2)
  sumstat$gradient <- 'gradient'
  
  summary <- sumstat
  
  # read
  directory_pow<-file.path(directory, 'power')
  sumstat<-combine_data(directory_pow, prefix,  suffix, 'central_qr_', splits)
  sumstat$qr <- 'central'
  sumstat2<-combine_data(directory_pow, prefix,  suffix, 'fed_qr_', splits )
  sumstat2$qr <- 'federated'
  sumstat<-rbind(sumstat, sumstat2)
  sumstat$gradient <- 'power'
  summary<- rbind(summary, sumstat)
  return(summary)
}

make_plot<-function(sumstat, evrank, lim, splitl){
  angles.plot<-ggplot(sumstat[splits %in% splitl & rank==evrank], aes(it, mean_angle, col=dataset))+
    geom_line(aes(linetype=gradient))+
    facet_wrap(~splits+qr, scales = 'free')+
    xlim(lim)+
    my_theme+ylab('Mean angle [degree]')+ 
    xlab('# iterations')+
    scale_color_manual('Algorithm', values = palette_div[c(1,3,8)])+
    theme(axis.line=element_line(),
          strip.background = element_blank(),
          legend.position = c(1, 1),
          legend.justification = c("right", "top"),
          legend.box.just = "right",
          legend.margin = margin(0.25, 0.25, 0.25, 0.25),
          legend.title = element_text(size=10))+
    guides(linetype=guide_legend(keyheight = 0.5, title = element_text('QR version', size = 8)),
           color=guide_legend(keyheight = 0.5, title = element_text('Algorithm', size = 8)))
  return(angles.plot)
}

make_plot_all<-function(sumstat, lim, splitl){
  angles.plot<-ggplot(sumstat[splits %in% splitl], aes(it, mean_angle, col=dataset))+
    geom_line(aes(linetype=gradient))+
    facet_grid(c( 'rank', 'qr'), scales = 'free')+
    xlim(lim)+
    my_theme+ylab('Mean angle [degree]')+ 
    xlab('# iterations')+
    scale_color_manual('Algorithm', values = palette_div[c(1,3,8)])+
    theme(axis.line=element_line(),
          strip.background = element_blank(),
          legend.position = c(1, 1),
          legend.justification = c("right", "top"),
          legend.box.just = "right",
          legend.margin = margin(0.25, 0.25, 0.25, 0.25),
          legend.title = element_text(size=10))+
    guides(linetype=guide_legend(keyheight = 0.5, title = element_text('QR version', size = 8)),
           color=guide_legend(keyheight = 0.5, title = element_text('Algorithm', size = 8)))
  return(angles.plot)
}

#### GENERIC ####
dataset <- 'mnist_4'
gradient <- 'gradient'
directory <-file.path('/home/anne/Documents/featurecloud/pca/vertical-pca/results/', dataset, 'horizontal')
prefix = 'raw'
bound = c(0,120)
figure_folder = '/home/anne/Documents/featurecloud/pca/vertical-pca/figures/'
splits <- c(2)
suffix = 'angles.u'

left.angles<-create_summary(directory, prefix, suffix, splits)

suffix = 'angles.v'
right.angles<-create_summary(directory, prefix, suffix, splits)

p.left<-make_plot_all(left.angles[dataset=='Matrix'], bound, splits)
p.right<-make_plot_all(right.angles, bound, splits)
p.right
p.left
plist<-list()

for (i in 1:10){
  p<-make_plot(left.angles, i, bound, splits)
  #ggsave(p, file=paste0(figure_folder,
                        #dataset, '/', suffix, '_eigenvector_comparison', i, '.pdf'),height = 5, units = 'cm')
  plist[[i]]<-p
}


# transform data to wide format for tikz

summarised<-left.angles %>% pivot_wider(values_from = 'mean_angle', names_from = c('rank', 'dataset', 'splits', 'qr', 
                                                                       'gradient') , names_sep = '.')

