require(data.table)
require(ggplot2)
require(dplyr)
require(tidyr)
require(tidyverse)
require(cowplot)
require(gridExtra)
library(ggpubr)
require(facetscales)
library(grid)
library(gtable)


process_angles<-function(mat_file, vec_file, name){
  data1<-fread(file.path(mat_file, 'summary.angles.u.tsv'))
  colnames(data1)<- c('iterations', paste0('Eigenvector ', 1:10), 'sites', 'algorithm', 'maxit', 'filename')
  data1<-data1[!is.na(data1$iterations)]
  data1$index<-1:nrow(data1)
  # get last iteration and index
  max_rows <- data1 %>% group_by( sites, algorithm, maxit, filename)  %>% slice_max(iterations)
  max_rows<-as.data.table(max_rows)
  iterations<-max_rows %>%ungroup() %>% select(iterations, sites, algorithm)
  
  
  globalmax<- data1[,max(iterations)]
  repeats<- abs(max_rows$iterations-globalmax)
  
  
  index<-unlist(sapply(max_rows$iterations[max_rows$iterations!=globalmax], function(x) seq(from=x+1, globalmax)))
  max_rows <- data1[rep(x = max_rows$index, repeats)]
  
  max_rows$iterations<-index
  max_rows<- as.data.table(max_rows)
  data1 <- rbind(data1, max_rows)
  
  summary1 <- data1  %>% group_by(iterations, sites,algorithm) %>% summarise(across(c(paste0('Eigenvector ', 1:10)),mean)) %>% pivot_longer(-c(iterations, sites,algorithm))
  summary1<-as.data.table(summary1)
  summary1$algorithm<-as.factor(summary1$algorithm)
  colnames(summary1)<-c('iterations', 'sites', 'algorithm', 'eigenvector', 'value')
  
  data<-fread(file.path(vec_file, 'summary.angles.u.tsv'))
  colnames(data)<-  c('eigenvector', 'iterations', 'value', 'sites','algorithm', 'dummy', 'filename')
  max_it<- data %>% group_by( sites, algorithm, filename)  %>% slice_max(iterations) %>% ungroup() %>%  select(iterations, sites, algorithm)
  data<- data  %>%group_by(eigenvector, filename) %>% mutate(counter = row_number(iterations))
  offset<-data %>% group_by(eigenvector, sites, algorithm, filename)%>% summarise(offset=max(iterations)) %>% summarise(mean_offset = floor(mean(offset)))
  offset<- as.data.table(offset)
  offset$eigenvector<-offset$eigenvector+1
  offset<-rbind(offset, list(eigenvector=1,sites=5, algorithm='guo', mean_offset = 0))
  #offset<-rbind(offset, list(eigenvector=1,sites=10, algorithm='guo', mean_offset = 0))
  data <- data %>% group_by(counter, eigenvector, sites, algorithm)%>% summarise(avg_angle = mean(value))
  data<-data %>% left_join(offset)
  
  
  data<-as.data.table(data)
  
  do<-list()
  counter <-1
  offset=offset[order(eigenvector)]
  for (i in 2:nrow(offset)){
    print(i)
    do[[counter]]<-data.table(c(1:as.numeric(offset[i,4])), offset[i,1], offset[i,2], 'guo', 90, '1')
    counter<-counter+1
  }
  
  
  do<-rbindlist(do)
  
  data$counter<-data$counter +data$mean_offset
  colnames(do)<-colnames(data)
  
  data<-rbind(data, do)
  colnames(data)<-c('iterations', 'eigenvector', 'sites', 'algorithm', 'value', 'maxit')
  data$eigenvector<-paste0('Eigenvector ', data$eigenvector)
  data <-data %>% select(iterations, eigenvector, sites, algorithm, value)
  summary1<-rbind(summary1, data)
  summary1$algorithm<-recode(summary1$algorithm, 'approximative-init'='AI-FULL', 
                             'random-init'='RI-FULL', 'randomized-no-approx'='RI-RAND', 
                             'randomized-approx-projected'='AI-RAND',
                             'guo'='GUO')
  summary1 <- summary1[algorithm %in% c('AI-FULL', 'RI-FULL', 'AI-RAND', 'RI-RAND', 'GUO')]
  summary1$eigenvector<-ordered(summary1$eigenvector,  levels = paste0('Eigenvector ', 1:10))
  
  summary1$name<-name
  
  iterations<-rbind(as.data.table(iterations), as.data.table(max_it))
  iterations$algorithm<-recode(iterations$algorithm, 'approximative-init'='AI-FULL', 'random-init'='RI-FULL', 'randomized-no-approx'='RI-RAND', 'randomized-approx-projected'='AI-RAND', 'guo'='GUO')
  iterations <- iterations[algorithm %in% c('AI-FULL', 'RI-FULL', 'AI-RAND', 'RI-RAND', 'GUO') ]
  iterations$name<-name

  return(list(summary1, iterations))
}

breaks_fun <- function(x) {
  if (max(x)>1001){
    seq(500, max(x), 500)
  }
  else if (max(x) > 500) {
    seq(250, max(x), 250)
  } else if (max(x)> 250) {
    seq(100, max(x), 150)
  }
  else{
    seq(50, max(x), 100)
  }
}

scales_x <- list(
  `Eigenvector 1` = scale_x_continuous(limits = c(0, 30)),
  `Eigenvector 5` = scale_x_continuous(limits = c(0, 300)),
  `Eigenvector 10` =scale_x_continuous(limits = c(0,2000))
  
)

extra_theme <-theme(axis.text.y = element_text( size=8),
                    axis.title = element_text(size=9),
                    plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),
                    legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
                    axis.text.x = element_blank(),
                    legend.text = element_text(size = 8))

#ggplot(summary1, aes(iterations, value,color=algorithm))+geom_line(aes(linetype=algorithm))+facet_grid(c( 'eigenvector','sites' ))
mat_file.100<-'~/Documents/featurecloud/pca/approximative-vertical/res-prelim/chr1.100000/vertical/matrix'
vec_file.100<-'~/Documents/featurecloud/pca/approximative-vertical/res-prelim/chr1.100000/vertical/vector'
res.100<-process_angles(mat_file = mat_file.100 , vec_file = vec_file.100, name = '100000')
convergence.100000<-res.100[[1]]

mat_file.500<-'~/Documents/featurecloud/pca/approximative-vertical/res-prelim/chr1.500000/vertical/matrix'
vec_file.500<-'~/Documents/featurecloud/pca/approximative-vertical/res-prelim/chr1.500000/vertical/vector'
res.500<-process_angles(mat_file = mat_file.500 , vec_file = vec_file.500, name = '500000')
convergence500000<-res.500[[1]]

mat_file.all<-'~/Documents/featurecloud/pca/approximative-vertical/res-prelim/chr1.all/vertical/matrix'
vec_file.all<-'~/Documents/featurecloud/pca/approximative-vertical/res-prelim/chr1.all/vertical/vector'
res.all<-process_angles(mat_file = mat_file.all , vec_file = vec_file.all, name = '1069419')
convergence.all<-res.all[[1]]

scales_x <- list(
  `Eigenvector 1` = scale_x_continuous(limits = c(0, 10)),
  `Eigenvector 5` = scale_x_continuous(limits = c(0, 300)),
  `Eigenvector 10` =scale_x_continuous(limits = c(0,2000))
  
)

breaks_fun <- function(x) {
  if (max(x)>1000){
    seq(500, max(x), 500)
  }
  else if (max(x) > 500) {
    seq(250, max(x), 250)
  } else if (max(x)> 250) {
    seq(100, max(x), 150)
  }
  else{
    seq(50, max(x), 100)
  }
}

data<-rbind(convergence.100000, convergence500000, convergence.all)

data$iterations<-as.integer(data$iterations)
colnames(data)<-c('iterations', 'sites', 'Algorithm', 'eigenvector', 'value', 'Features')
data$Features<- ordered(data$Features, levels=c('100000', '500000', '1069419'))
conv<-ggplot(data[sites ==5&eigenvector %in% c('Eigenvector 1', 'Eigenvector 5', 'Eigenvector 10')], 
             aes(iterations, value, group_=interaction(Algorithm, Algorithm), color=Algorithm))+
  geom_line(aes(linetype=Algorithm),size=0.5)+
  xlab('Iterations')+
  theme_classic()+
  ylab('Angles w.r.t reference')+
  facet_grid_sc(rows=vars(Features), cols = vars(eigenvector), scales = list(x = scales_x))+
  scale_x_continuous(breaks=breaks_fun, limits = c(0, NA))+
  theme(legend.position = 'None',
        strip.text = element_text(size=8), panel.border = element_rect(fill=NA), 
        panel.spacing = unit(0.05, 'cm'),strip.background = element_blank(), axis.title.x = element_blank() )+
  scale_color_viridis_d(option = 'A', end=0.7)+extra_theme

# conv<-ggplot(data[sites ==5&eigenvector %in% c('Eigenvector 1', 'Eigenvector 5', 'Eigenvector 10')],
#              aes(iterations, value, group_=interaction(Algorithm, Algorithm), color=Algorithm))+
#   geom_line(size=0.5)+
#   geom_point(data=data[sites ==5&eigenvector %in% c('Eigenvector 1', 'Eigenvector 5', 'Eigenvector 10')][iterations %% 100 ==0],
#              aes(x = iterations, y = value ,shape=Algorithm))+
#   xlab('Iterations')+
#   theme_classic()+
#   ylab('Angles w.r.t reference')+
#   facet_grid_sc(rows=vars(Features), cols = vars(eigenvector), scales = list(x = scales_x))+
#   scale_x_continuous(breaks=breaks_fun, limits = c(0, NA))+
#   theme(legend.position = 'None',
#         strip.text = element_text(size=8), panel.border = element_rect(fill=NA),
#         panel.spacing = unit(0.05, 'cm'),strip.background = element_blank(), axis.title.x = element_blank() )+
#   scale_color_viridis_d(option = 'A', end=0.7)+extra_theme
conv

# Convert the plot to a grob
text = 'Chromosome 1'
size=10
gt <- ggplotGrob(conv)
# Get the positions of the right strips in the layout: t = top, l = left, ...
strip <-c(subset(gt$layout, grepl("strip-r", gt$layout$name), select = t:r))
# Text grob
text.grob = textGrob(text, rot = -90,
                     gp = gpar(fontsize = size))
# New olumn to the right of current strip
# Adjusts its width to text size
width = unit(2, "grobwidth", text.grob) + unit(1, "lines")
gt <- gtable_add_cols(gt, width, max(strip$r))  
# Add text grob to new column
gt <- gtable_add_grob(gt, text.grob, 
                      t = min(strip$t), l = max(strip$r) + 1, b = max(strip$b))
# Draw it
conv<-ggplotify::as.ggplot(gt)

c<-conv+theme()
c
ggsave(conv,file='/home/anne/Documents/manuscripts/vertical-pattern-recognition/figures/chrom1_convergence.pdf', height=7, units = 'cm', width = 15 )


extra_theme <-theme(axis.text.y = element_text( size=8),
      axis.title.y = element_text(size=9),
      plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),
      axis.text.x = element_blank())


iterations<-rbind(res.100[[2]], res.500[[2]], res.all[[2]])
pl.iter<-ggplot(iterations, aes(algorithm, iterations, fill=interaction(algorithm, name),color=algorithm))+
  geom_boxplot()+xlab('')+
  theme_classic()+
  ylab('Iterations to convergence')+
  guides(fill=FALSE, color=F)+
  scale_color_viridis_d(option = 'A', end=0.7)+
  scale_fill_manual(values=c(viridisLite::magma(5, end=0.7, alpha = 0.1), viridisLite::magma(5, end=0.7, alpha = 0.4)
                             , viridisLite::magma(5, end=0.7, alpha = 0.7)))+
  extra_theme
pl.iter




time_reader<-function(mat_file, name){
  time <- fread(file.path(mat_file, '/../time.log'))
  time$V1<-recode(time$V1, 'approximative-init'='AI-FULL', 'random-init'='RI-FULL', 'randomized-no-approx'='RI-RAND', 'randomized-approx-projected'='AI-RAND', 'guo'='GUO')
  time<- time[V1 %in%  c('AI-FULL', 'RI-FULL', 'AI-RAND', 'RI-RAND', 'GUO')]
  time$name<-name
  return(time)
}
time.100<-time_reader(mat_file.100, '100000')
time.500 <-time_reader(mat_file.500 , '500000')
time.all<-time_reader(mat_file.all, '1069419')
time<-rbind(time.100, time.500, time.all)

pl.runtime<-ggplot(time, aes(V1, V4, fill=interaction(V1, name)))+geom_boxplot()+ylab('Runtime[s]')+xlab('')+
  theme_classic()+
  scale_fill_manual(values=c(viridisLite::magma(5, end=0.7, alpha = 0.1), viridisLite::magma(5, end=0.7, alpha = 0.4)
                             , viridisLite::magma(5, end=0.7, alpha = 0.7)))+guides(color=FALSE)+
  extra_theme

pl.runtime

eo_reader<-function(mat_file, vec_file, name){
  eo<- fread(file.path(mat_file, 'summary.eo.tsv'))
  eog<-fread(file.path(vec_file, 'summary.eo.tsv'))
  eo<-rbind(eo, eog)
  eo$V4<-recode(eo$V4, 'approximative-init'='AI-FULL', 'random-init'='RI-FULL', 'randomized-no-approx'='RI-RAND', 'randomized-approx-projected'='AI-RAND', "guo"="GUO")
  eo<- eo[V4 %in% c('AI-FULL', 'RI-FULL', 'AI-RAND', 'RI-RAND', 'GUO') ]
  eo <- eo %>% group_by(V3,V4, V5, V6) %>% summarise(sum=sum(V2))
  eo<-as.data.table(eo)
  eo$name<-name 
  return(eo)
}

eo.100<-eo_reader(mat_file.100, vec_file.100, '100000')
eo.500<-eo_reader(mat_file.500, vec_file.500, '500000')
eo.all<-eo_reader(mat_file.all, vec_file.all, '1069419')
eo<-rbind(eo.100, eo.500, eo.all)

pl.eo<-ggplot(eo, aes(V4, sum, fill=interaction(V4,name), color=V4))+geom_boxplot()+ylab('Runtime[s] (Matrix oper.)')+xlab('')+
  theme_classic()+theme(axis.text.x = element_text(angle = 90) )+guides(fill=FALSE, color=F)+scale_color_viridis_d(option = 'A', end=0.7)+
  scale_fill_manual(values=c(viridisLite::magma(5, end=0.7, alpha = 0.1), viridisLite::magma(5, end=0.7, alpha = 0.4)
                             , viridisLite::magma(5, end=0.7, alpha = 0.7)))+
  extra_theme
pl.eo





transmission_read<-function(mat_file, vec_file, name){
  transmission<-fread(file.path(mat_file, 'summary.transmission.tsv'))
  transmission_guo<-fread(file.path(vec_file, 'summary.transmission.tsv'))
  
  transmission<-rbind(transmission_guo, transmission)
  transmission<-transmission[V2 %in% c("H_local=CS", "H_global=SC")]
  colnames(transmission) <- c('index', 'key', 'iteration', 'site', 'size', 'kA', 'sites', 'algorithm', 'maxit', 'filename')
  transmission<-transmission[site==1]
  volume<-transmission %>% select(size, sites, algorithm, filename)%>% group_by(sites, algorithm, filename)%>% summarise(total_volume=(sum(size)*4)/1000000000)
  
  volume<-as.data.table(volume)
  volume$algorithm<-recode(volume$algorithm, 'approximative-init'='AI-FULL', 'random-init'='RI-FULL', 'randomized-no-approx'='RI-RAND', 'randomized-approx-projected'='AI-RAND', 'guo'='GUO')
  volume<- volume[algorithm %in% c( 'RI-FULL', 'AI-RAND', 'RI-RAND', 'AI-FULL', 'GUO')]
  volume$name <- name
  return(volume)
}

volume.100<-transmission_read(mat_file.100, vec_file.100, '100000')
volume.500<-transmission_read(mat_file.500, vec_file.500, '500000')
volume.all<-transmission_read(mat_file.all, vec_file.all, '1069419')
volume<-rbind(volume.100, volume.500, volume.all)
volume$name<- ordered(volume$name, levels=c('100000', '500000', '>1Mio'))
############FAKE ALPHA LEGEND
pl.volume.legend<-ggplot(volume, aes(algorithm, total_volume, fill= algorithm, alpha=name))+geom_boxplot()+
  theme_classic()+ylab('Total transmitted\n data [GB]')+xlab('')+
  scale_alpha_manual(values=c(0.2,0.4,0.7))+
  guides(fill=F,alpha=guide_legend('#Features', override.aes=list(fill=hcl(c(15,195),100,0,alpha=c(0.2,0.4,0.7)))))
pl.volume.legend

pl.volume<-ggplot(volume, aes(algorithm, total_volume, fill=interaction(algorithm, name), color=algorithm))+
  geom_boxplot()+
  ylab('Total transmitted data [GB]')+
  xlab('')+
  theme_classic()+
  guides(fill=FALSE, color=F)+
  scale_color_viridis_d(option = 'A', end=0.7)+
  scale_fill_manual(values=c(viridisLite::magma(5, end=0.7, alpha = 0.2), 
                             viridisLite::magma(5, end=0.7, alpha = 0.4)
                             , viridisLite::magma(5, end=0.7, alpha = 0.7)))+extra_theme
  
pl.volume

legend_b <- get_legend(
  pl.volume.legend +theme(legend.position = 'bottom',
                          legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0),
                          legend.text = element_text(size = 8),
                          legend.title = element_text(size = 8))
)
text = 'Chromosome 1'
size=12
text.grob = textGrob(text, rot = -90,
                     gp = gpar(fontsize = size))
gp<-ggarrange(pl.iter, pl.eo, pl.volume, ncol=3)
gp<-plot_grid(gp, text.grob, rel_widths = c(1,0.05))
#gp<-plot_grid(gp, legend_b, ncol=1, rel_heights = c(1,0.1))
# Convert the plot to a grob
gp
ggsave(gp,file='/home/anne/Documents/manuscripts/vertical-pattern-recognition/figures/chrom1_benchmark.pdf', height=5, units = 'cm', width = 15)



# mev<-fread(file.path(mat_file,'summary.mev.u.tsv'))
# colnames(mev)<-c('iterations', 'mev', 'sites', 'algorithm', 'maxit', 'filename')
# mev <- mev %>% group_by(iterations, sites, algorithm) %>% summarise(mean_mev = mean(mev))
# mev<-as.data.table(mev)
# mev$algorithm<-recode(mev$algorithm, 'approximative-init'='AI-FULL', 'random-init'='RI-FULL', 'randomized-no-approx'='RI-RAND', 'randomized-approx-projected'='AI-RAND', 'guo'='GUO')
# mev<- mev[algorithm %in% c( 'RI-FULL', 'AI-RAND', 'RI-RAND', 'AI-FULL') ]
# ggplot(mev, aes(iterations,mean_mev, color=algorithm))+geom_line()+facet_wrap(~sites)+
#   xlab('MEV')+theme_classic()+ylab('iterations')+xlim(0,100)+scale_color_viridis_d(option = 'A', end=0.7)
# 
# mev<-fread(file.path(vec_file, 'summary.mev.u.tsv'))
# colnames(mev)<-c('vector','iterations', 'mev', 'sites', 'algorithm', 'maxit', 'filename')
# mev <- mev %>% group_by(iterations, sites, algorithm) %>% summarise(mean_mev = mean(mev))
# mev<-as.data.table(mev)
# mev$algorithm<-recode(mev$algorithm, 'approximative-init'='AI-FULL', 'random-init'='RI-FULL', 'randomized-no-approx'='RI-RAND', 'randomized-approx-projected'='AI-RAND', 'guo'='GUO')
# ggplot(mev, aes(iterations,mean_mev, color=algorithm))+geom_line()+facet_wrap(~sites)+
#   xlab('MEV')+theme_classic()+ylab('iterations')+xlim(0,100)+scale_color_viridis_d(option = 'A', end=0.7)
# 
