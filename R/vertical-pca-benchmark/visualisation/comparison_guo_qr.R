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
  
  read_files_from_dir <- function(directory, prefix, suffix, splits) {
    myfiles <-  list.files(directory, recursive = F)
    myangles <- list()
    mylastlines <- list()
    for (i in splits) {
      for (f in which(str_detect(myfiles, paste0("^", prefix, i, ".*.", suffix, '$')))) {
        filename <- myfiles[f]
        if(countLines(file.path(directory, filename))<=9){
          next
        }
        print(file.path(directory, filename))
        angles <-
          fread(
            file.path(directory, filename),
            sep = '\t',
            skip = 9,
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
        if(countLines(file.path(directory, filename))<=9){
          next
        }
        
        angles <-
          fread(
            file.path(directory, filename),
            sep = '\t',
            skip = 9,
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
  
  my_theme <-
    theme_classic() + theme(
      axis.title = element_text(size = 12),
      legend.title = element_text(size = 12),
      plot.title = element_text(size = 20, hjust = 0.5),
      axis.text = element_text(size = 8),
      legend.text = element_text(size = 8),
      plot.subtitle = element_text(size = 10, hjust = 0.5),
      legend.key.size = unit(1.5, 'lines'))
  
  make_plot<-function(sumstat, evrank, lim, splitl){
    angles.plot<-ggplot(sumstat[splits %in% splitl & rank==evrank], aes(it, mean_angle, col=dataset))+
      geom_line(aes(linetype=qr))+
      facet_wrap(~splits, scales = 'free')+
      xlim(lim)+
      my_theme+ylab('Mean angle [degree]')+ 
      xlab('# iterations')+
      scale_color_manual('Algorithm', values = palette_div[c(1,3,8)])+
      theme(axis.line=element_line(),
            strip.background = element_blank(),
            strip.text.x = element_blank(),
            legend.position = c(1, 1),
            legend.justification = c("right", "top"),
            legend.box.just = "right",
            legend.margin = margin(0.25, 0.25, 0.25, 0.25),
            legend.title = element_text(size=10))+
      guides(linetype=guide_legend(keyheight = 0.5, title = element_text('QR version', size = 8)),
             color=guide_legend(keyheight = 0.5, title = element_text('Algorithm', size = 8)))
    return(angles.plot)
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
  angles_hybrid<-read_files_from_dir(directory, paste0(prefix,  '_hybrid_', qr), suffix, splits)
  angles_reit<- read_files_from_dir_indiv(directory, paste0(prefix,  '_hybrid_reit_', qr), suffix, splits)
  
  # find the last iteration per data set in the hybrid setting
  final_it <- angles_hybrid %>% group_by(filename, splits) %>% summarise(last_it = max(it))
  final_it<-as.data.table(final_it)
  
  # correct the number of iterations in the hybrid setting
  if(!is.null(angles_reit) && nrow(angles_reit)!=0){
    angles_reit$filename<-sapply(angles_reit$filename, function(x) gsub('_reit', '', x))
    angles_reit <- left_join(angles_reit, final_it, by = c('filename', 'splits'))
    angles_reit$it <- angles_reit$it+angles_reit$last_it
    angles_reit<- angles_reit[,1:6]
    angles_reit$dataset<-'raw_hybrid_'
    colnames(angles_reit)<-c('rank', 'it', 'angle', 'filename','dataset','splits')
    hybr <- rbind(angles_hybrid, angles_reit)
  }else{
    hybr <- angles_hybrid
  }

  angles_guo$dataset <- 'Vector'
  angles_qr$dataset<-'Matrix'
  hybr$dataset<-'Hybrid'

  # combine the data sets
  a <- rbind(angles_guo, angles_qr)
  a <- rbind(a, hybr)
  a$angle<-abs(a$angle)
  sumstat <- a %>% group_by(rank, it, dataset, splits) %>% summarise(mean_angle=mean(angle))
  sumstat <- as.data.table(sumstat)
  sumstat$rank<-as.numeric(sumstat$rank)
  
  return(sumstat)
}

#dataset <- 'mfeat'
#directory <- paste0('/home/anne/Documents/featurecloud/pca/vertical-pca/results/', dataset)
#suffix = 'angles'
#prefix = 'mfeat\\-zer'
#bound = c(0,100)




# dataset <- 'MMRF-COMMPASS'
# directory <- paste0('/home/anne/Documents/featurecloud/pca/vertical-pca/results/', dataset)
# suffix = 'angles'
# prefix = 'MMRF-COMMPASS'
# bound = c(0,300)
# figure_folder = '/home/anne/Documents/featurecloud/pca/vertical-pca/figures/'
# splits <- c(2,3,5,10)



#### GENERIC ####
dataset <- 'mnist'
directory <- paste0('/home/anne/Documents/featurecloud/pca/vertical-pca/results/', dataset)

prefix = 'raw'
bound = c(0,120)
figure_folder = '/home/anne/Documents/featurecloud/pca/vertical-pca/figures/'
splits <- c(2,5,10)
suffix = 'angles'
sumstat<-combine_data(directory, prefix,  suffix, 'central_qr_', splits)
sumstat$qr <- 'central'
sumstat2<-combine_data(directory, prefix,  suffix, 'fed_qr_', splits )
sumstat2$qr <- 'federated'
sumstat<-rbind(sumstat, sumstat2)
for (i in 1:10){
  p<-make_plot(sumstat, i, bound, splits)
  ggsave(p, file=paste0(figure_folder,
                              dataset, '/', suffix, '_eigenvector_comparison', i, '.pdf'),height = 5, units = 'cm')
}

######        MNIST ###
p1<-make_plot(sumstat, 1, c(0,20), splits)
ggsave(p1, file=paste0(figure_folder,
                       dataset, '/', suffix, '_eigenvector_comparison_1_zi.pdf'), 
       height = 5, units = 'cm')

p8<-make_plot(sumstat, 8, c(0,300), splits)
ggsave(p8, file=paste0(figure_folder,
                       dataset, '/', suffix, '_eigenvector_comparison_8_zo.pdf'), height = 5, units = 'cm')



### CHROMOSME 1####
#dataset <- '1000g/chr1'
#prefix = 'chr1'
dataset <- '1000g/chr2'
prefix = 'chr2'

directory <- paste0('/home/anne/Documents/featurecloud/pca/vertical-pca/results/', dataset)
suffix = 'angles'

bound = c(0,500)
figure_folder = '/home/anne/Documents/featurecloud/pca/vertical-pca/figures/'
splits <- c(2,5,10)

suffix = 'angles'
sumstat<-combine_data(directory, prefix,  suffix, 'central_qr_', splits)
sumstat$qr <- 'central'
sumstat2<-combine_data(directory, prefix,  suffix, 'fed_qr_', splits )
sumstat2$qr <- 'federated'
sumstat<-rbind(sumstat, sumstat2)
for (i in 1:5){
  bound = c(0,30)
  p<-make_plot(sumstat, i, bound, splits)
  ggsave(p, file=paste0(figure_folder,
                        dataset, '/', suffix, '_eigenvector_comparison', i, '.pdf'),height = 5, units = 'cm')
}
for (i in 6:7){
  bound = c(0,50)
  p<-make_plot(sumstat, i, bound, splits)
  ggsave(p, file=paste0(figure_folder,
                        dataset, '/', suffix, '_eigenvector_comparison', i, '.pdf'),height = 5, units = 'cm')
}
for (i in 7:10){
  bound = c(0,500)
  p<-make_plot(sumstat, i, bound, splits)
  ggsave(p, file=paste0(figure_folder,
                        dataset, '/', suffix, '_eigenvector_comparison', i, '.pdf'),height = 5, units = 'cm')
}



 

# suffix = 'cor'
# sumstat$qr <- 'central'
# sumstat2<-combine_data(directory, prefix,  suffix, 'fed_qr_', splits )
# sumstat2$qr <- 'federated'
# sumstat<-rbind(sumstat, sumstat2)
# 
# for (i in 1:10){
#   p<-make_plot_correlation(sumstat, i, bound,  splits)
#   ggsave(p, file=paste0(figure_folder,
#                              dataset, '/', suffix, '_eigenvector_comparison', i, '.pdf'),height = 5, units = 'cm')
# }




##UNEQUAL###
dataset <- 'MMRF-COMMPASS'
directory <- paste0('/home/anne/Documents/featurecloud/pca/vertical-pca/results/unequal/', dataset)
suffix = 'angles'
prefix = "MMRF-COMMPASS"
bound = c(0,120)
figure_folder = '/home/anne/Documents/featurecloud/pca/vertical-pca/figures/unequal/'
splits <- 1:5

suffix = 'angles'
sumstat<-combine_data(directory, prefix,  suffix, 'central_qr_', splits)
sumstat$qr <- 'central'
sumstat2<-combine_data(directory, prefix,  suffix, 'fed_qr_', splits )
sumstat2$qr <- 'federated'
sumstat<-rbind(sumstat, sumstat2)
sumstat$splits <- ((sumstat$splits+1)%/%20)+1
sumstat[,mean_angle:=mean(mean_angle), by = c('rank', 'it', 'dataset', 'splits', 'qr')]
for (i in 1:10){
  p<-make_plot(sumstat, i, bound, splits)
  ggsave(p, file=paste0(figure_folder,
                        dataset, '/', suffix, '_eigenvector_comparison', i, '.pdf'),height = 5, units = 'cm')
}
