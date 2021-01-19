require(data.table)
require(ggplot2)
require(tidyr)
require(cowplot)
require(optparse)
require(dplyr)

palette_div <-
  c(
    '#062865',
    '#2f497d',
    '#4c6d96',
    '#6793af',
    '#81bac8',
    '#ffc4b3',
    '#f68888',
    '#d75161',
    '#ab203f',
    '#720022'
  )
palette_seq <-
  c(
    '#062865',
    '#203a72',
    '#324d80',
    '#43618d',
    '#52759b',
    '#618aa9',
    '#70a0b7',
    '#7eb6c5',
    '#8dcdd4',
    '#9ce4e2'
  )

my_theme <-
  theme_classic() + theme(
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 12),
    plot.title = element_text(size = 20, hjust = 0.5),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 12),
    plot.subtitle = element_text(size = 12, hjust = 0.5)
  )

#
read_files_from_dir <- function(directory, suffix) {
  myfiles <-  list.files(directory)
  myangles <- list()
  mylastlines <- list()
  for (i in c(2, 3, 5, 10)) {
    for (f in grep(paste0("^raw_", i, ".*.", suffix), myfiles, perl = T)) {
      filename <- myfiles[f]
      
      print(file.path(dir, filename))
      angles <-
        fread(
          file.path(dir, filename),
          sep = '\t',
          skip = 9,
          header = F
        )
      if (is.null(angles)) {
        next
      }
      colnames(angles) <-
        c('iterations', sapply(1:(ncol(angles) - 1), function(x)
          paste0('', x)))
      angles$dataset <- filename
      angles$splits <- i
      
      # also save last line in a separate table
      mylastlines[[filename]] <- angles[nrow(angles)]
      angles <-
        as.data.table(pivot_longer(angles,-c(iterations, dataset, splits)))
      myangles[[filename]] <- angles

      
    }
  }
  
  mylastlines <- rbindlist(mylastlines)
  mylastlines <- mylastlines[, 2:ncol(mylastlines)]
  myangles <- rbindlist(myangles)
  myangles <- myangles[!is.na(value)]
  
  return (list(myangles, mylastlines))
}

summarise_angle<-function(myangles){
  myangles = as.data.table(
    myangles %>% group_by(iterations, name, splits) %>% summarise(
      median = median(abs(value)),
      q1 = quantile(abs(value), 0.25),
      q3 = quantile(abs(value), 0.75)
    )
  )
  lev <-sapply(1:(length(unique(myangles$name))), function(x) paste0('', x))
  myangles$name <- factor(myangles$name, levels = lev)
  return(myangles)
}
make_angle_vs_iteration_ribbon_plot<-function(myangles,ev, sp, dir, prefix){
p1<-ggplot(myangles[name %in% ev & splits %in% sp  & iterations<60], aes(x = iterations, y = median, color = as.factor(name))) +
  geom_line() +
  geom_ribbon(aes(
    x = iterations,
    ymin = q1,
    ymax = q3,
    color=as.factor(name),
    fill =as.factor(name)
  ) ,
  alpha = 0.1) +
  scale_fill_manual(name = "Eigenvector rank", values = palette_div[c(1,3, 5,7, 9)]) +
  scale_alpha_manual( values = palette_div[c(1,3, 5,7, 9)]) +
  scale_color_manual( values = palette_div[c(1,3, 5,7, 9)]) +
  my_theme +
  ylab('Angle [degree]') + xlab('#Iterations') + guides(color = guide_legend(title = "Eigenvector rank"))+
  theme(legend.position=c(.6,.95), legend.direction = 'horizontal')
  name = paste0(prefix, '_vs_iterations_ribbon_evs_',paste0(ev, collapse = '_'), '_splits_', sp, '.pdf')
  ggsave(p1, file=paste0(dir, name))
  p1
}





make_correlation_vs_iteration_ribbon_plot<-function(myangles,ev, sp, dir, prefix){
  p1<-ggplot(myangles[name %in% ev & splits %in% sp  & iterations<60], aes(x = iterations, y = median, color = as.factor(name))) +
    geom_line() +
    geom_ribbon(aes(
      x = iterations,
      ymin = q1,
      ymax = q3,
      color=as.factor(name),
      fill =as.factor(name)
    ) ,
    alpha = 0.1) +
    scale_fill_manual(name = "Eigenvector rank", values = palette_div[c(1,3, 5,7, 9)]) +
    scale_alpha_manual( values = palette_div[c(1,3, 5,7, 9)]) +
    scale_color_manual( values = palette_div[c(1,3, 5,7, 9)]) +
    my_theme +
    ylab('Pearson correlation') + xlab('#Iterations') + guides(color = guide_legend(title = "Eigenvector rank"))+
    theme(legend.position=c(.6,.15), legend.direction = 'horizontal')
  name = paste0(prefix, '_vs_iterations_ribbon_evs_',paste0(ev, collapse = '_'), '_splits_', sp, '.pdf')
  ggsave(p1, file=paste0(dir, name))
  p1
}
make_correlation_vs_iteration_plot<-function(myangles, ev, sp, dir, prefix){
  p1<-ggplot(myangles[name %in% ev & splits==sp & iterations<60], aes(x = iterations, y = median, color =as.factor(name))) +
    geom_line() + 
    my_theme  + 
    scale_color_manual(name = "Eigenvector rank",values = palette_div[c(1, 3,5,7, 9)])+
    ylab('Pearson correlation') + xlab('#Iterations') +theme(legend.position=c(.6,.15), legend.direction = 'horizontal')
  name = paste0(prefix,'_vs_iterations_evs_',paste0(ev, collapse = '_'), '_splits_', sp, '.pdf')
  ggsave(p1, file=paste0(dir, name))
}
angles_at_convergence<-function(mylastlines, dir){
  mylastlines <-as.data.table(pivot_longer(mylastlines,-c(dataset, splits)))
  mylastlines <- mylastlines[!is.na(value)]
  mylastlines <-mylastlines %>% group_by(name, splits) %>% summarise(median = median(value))
  lev <-sapply(1:(length(unique(mylastlines$name))), function(x) paste0('', x))
  mylastlines$name <- factor(mylastlines$name, levels = lev)
  
  p1<-ggplot(mylastlines, aes(x = name, y = median, col = as.factor(splits))) +
    geom_point() +
    scale_color_manual('Dataset partitions', values = palette_div[c(1, 3, 6, 9)]) +
    ylab('Median angle [degree]') + xlab('Eigenvector rank') +
    theme(legend.position='top', legend.direction = 'horizontal')
  name = paste0('angles_at_convergence',paste0(ev, collapse = '_'), '_splits.pdf')
  ggsave(p1, file=file.path(dir, name))
}

#Our transformation function
scaleFUN <- function(x) sprintf("%.3f", x)


correlation_at_convergence<-function(mylastlines, dir){
  mylastlines <-as.data.table(pivot_longer(mylastlines,-c(dataset, splits)))
  mylastlines <- mylastlines[!is.na(value)]
  mylastlines <-mylastlines %>% group_by(name, splits) %>% summarise(median = median(round(abs(value),2)))
  lev <-sapply(1:(length(unique(mylastlines$name))), function(x) paste0('', x))
  mylastlines$name <- factor(mylastlines$name, levels = lev)
  
  p1<-ggplot(mylastlines, aes(x = name, y = median, col = as.factor(splits))) +
    geom_point() +
    scale_color_manual('Dataset partitions', values = palette_div[c(1, 3, 6, 9)]) +
    ylab('Pearson correlation') + xlab('Eigenvector rank') +
    my_theme+
    theme(legend.position='top', legend.direction = 'horizontal')+
    scale_y_continuous(labels=scaleFUN)
  name <- paste0('correlation_at_convergence',paste0(ev, collapse = '_'), '_splits.pdf')
  ggsave(p1, file=file.path(dir, name))
  p1
}
make_angle_vs_iteration_plot<-function(myangles, ev, sp, dir, prefix){
  p1<-ggplot(myangles[name %in% ev & splits %in% sp & iterations<60], aes(x = iterations, y = median, color =as.factor(name))) +
    geom_line(size=1.5) + 
    my_theme  + 
    scale_color_manual(name = "Eigenvector rank",values = palette_div[c(1, 7)])+
    ylab('Median angle [degree]') + xlab('#Iterations') +
    theme(legend.position='top', legend.direction = 'horizontal')+
    facet_wrap(~splits)
  name = paste0(prefix,'_vs_iterations_evs_',paste0(ev, collapse = '_'), '_splits_', paste0(sp, collapse = '_'), '.pdf')
  ggsave(p1, file=paste0(dir, name))
  p1
}
  


dir <-'/home/anne/Documents/featurecloud/pca/vertical-pca/results/mnist'
outdir <- '/home/anne/Documents/featurecloud/pca/vertical-pca/figures/'
ev <-c('1', '10')

res <- read_files_from_dir(dir, 'angles')
myangles<-res[[1]]
myangles<-summarise_angle(myangles)
make_angle_vs_iteration_plot(myangles, ev, c(2,5,10), outdir, 'angles')

make_angle_vs_iteration_plot(myangles, ev, 5, outdir, 'angles')
make_angle_vs_iteration_ribbon_plot(myangles, ev, 5, outdir, 'angles')
mylastlines <- res[[2]]
angles_at_convergence(mylastlines = mylastlines , outdir)




res <- read_files_from_dir(dir, 'cor')
myangles<-res[[1]]
myangles<-summarise_angle(myangles)
make_correlation_vs_iteration_plot(myangles, ev, 2, outdir, 'correlation')
make_correlation_vs_iteration_plot(myangles, ev, 5, outdir, 'correlation')
make_correlation_vs_iteration_ribbon_plot(myangles, ev, 5, outdir, 'correlation')
mylastlines <- res[[2]]
correlation_at_convergence(mylastlines = mylastlines, outdir)




