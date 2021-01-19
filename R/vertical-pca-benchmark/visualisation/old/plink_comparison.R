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

read_files_from_dir <- function(directory, prefix, suffix) {
  myfiles <-  list.files(directory, recursive = F)
  myangles <- list()
  mylastlines <- list()
  for (i in c(1,2, 3, 5, 10)) {
    for (f in which(str_detect(myfiles, paste0("^", prefix, i, ".*.", suffix)))) {
      filename <- myfiles[f]
      if(countLines(file.path(directory, filename))<=9){
        next
      }
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

directory <- paste0('/home/anne/Documents/featurecloud/pca/vertical-pca/results/', dataset)
dataset = '1000g/chr1'
suffix = 'angles'
i = 'all'
prefix = 'chr1_fed_qr'

angles_plink <- read_files_from_dir(directory, paste0(prefix,  '_') , 'angles_precomp')
sumstat <- angles_plink %>% group_by(rank, it, dataset, splits) %>% summarise(mean_angle=mean(angle))
sumstat <- as.data.table(sumstat)
sumstat$rank<-as.numeric(sumstat$rank)
sumstat$qr <- 'federated'

sumstat2<-sumstat

prefix <- 'chr1_central_qr'
angles_plink <- read_files_from_dir(directory, paste0(prefix,  '_') , 'angles_precomp')
sumstat <- angles_plink %>% group_by(rank, it, dataset, splits) %>% summarise(mean_angle=mean(angle))
sumstat <- as.data.table(sumstat)
sumstat$rank<-as.numeric(sumstat$rank)
sumstat$qr <- 'central'

sumstat<-rbind(sumstat, sumstat2)

angles.plot<-ggplot(sumstat[splits %in% c(10)], aes(it, mean_angle, col=as.factor(rank)))+
  geom_line(size = 1) +facet_wrap(~qr, scales = 'free')+
  my_theme+ylab('Mean angle [degree]')+ 
  xlab('# iterations')+
  scale_color_manual('Eigenvector\nrank', values = palette_div)+
  theme(axis.line=element_line(),
        legend.position = c(1, 0.95),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6))+
  guides(linetype=FALSE, color=guide_legend(keyheight = 0.75, title = element_text('Eigenvector\nrank', size = 10)))
angles.plot
ggsave(angles.plot, file=paste0('/home/anne/Documents/featurecloud/pca/vertical-pca/figures/',
                      dataset, '/', suffix, '_eigenvector_comparison_plink_', 10, '.pdf'),height = 10, units = 'cm')

angles.plot<-ggplot(sumstat[splits %in% c(2,3,5,10)], aes(it, mean_angle, col=as.factor(rank)))+
  geom_line(size = 1) +facet_grid( rows = vars(splits), cols = vars(qr))+
  my_theme+ylab('Mean angle [degree]')+ 
  xlab('# iterations')+
  scale_color_manual('Eigenvector\nrank', values = palette_div)+
  theme(axis.line=element_line(),
        legend.position = c(1, 0.95),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6))+
  guides(linetype=FALSE, color=guide_legend(keyheight = 0.75, title = element_text('Eigenvector\nrank', size = 10)))
angles.plot
ggsave(angles.plot, file=paste0('/home/anne/Documents/featurecloud/pca/vertical-pca/figures/',
                                dataset, '/', suffix, '_eigenvector_comparison_plink_', i, '.pdf'),height = 10, units = 'cm')
