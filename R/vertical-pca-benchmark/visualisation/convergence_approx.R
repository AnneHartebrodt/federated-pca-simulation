require(data.table)
require(ggplot2)
require(dplyr)
require(tidyr)
require(tidyverse)
require(cowplot)


data1<-fread('~/Documents/featurecloud/pca/approximative-vertical/results/matrix/summary.angles.u.tsv')
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

summary1 <- data1  %>% group_by(iterations, sites,algorithm) %>% summarise(across(c(paste0('Eigenvector ', 1:10)),mean)) %>% 
  pivot_longer(-c(iterations, sites,algorithm))
summary1<-as.data.table(summary1)
summary1$algorithm<-as.factor(summary1$algorithm)
colnames(summary1)<-c('iterations', 'sites', 'algorithm', 'eigenvector', 'value')
#ggplot(summary1, aes(iterations, value,color=algorithm))+geom_line(aes(linetype=algorithm))+facet_grid(c( 'eigenvector','sites' ))





data<-fread('~/Documents/featurecloud/pca/approximative-vertical/results/vector/summary.angles.u.tsv')
colnames(data)<-  c('eigenvector', 'iterations', 'value', 'sites','algorithm', 'dummy', 'filename')
max_it<- data %>% group_by( sites, algorithm, filename)  %>% slice_max(iterations) %>% ungroup() %>%  select(iterations, sites, algorithm)
data<- data  %>%group_by(eigenvector, filename) %>% mutate(counter = row_number(iterations))
offset<-data %>% group_by(eigenvector, sites, algorithm, filename)%>% summarise(offset=max(iterations)) %>% summarise(mean_offset = floor(mean(offset)))
offset<- as.data.table(offset)
offset$eigenvector<-offset$eigenvector+1
offset<-rbind(offset, list(eigenvector=1,sites=5, algorithm='guo', mean_offset = 0))
offset<-rbind(offset, list(eigenvector=1,sites=10, algorithm='guo', mean_offset = 0))
data <- data %>% group_by(counter, eigenvector, sites, algorithm)%>% summarise(avg_angle = mean(value))
data<-data %>% left_join(offset)


data<-as.data.table(data)

do<-list()
counter <-1
offset<-offset[order(eigenvector)]
for (i in 3:20){
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

scales_x <- list(
  `Eigenvector 1` = scale_x_continuous(limits = c(0, 30)),
  `Eigenvector 5` = scale_x_continuous(limits = c(0, 300)),
  `Eigenvector 10` =scale_x_continuous(limits = c(0,1100))
  
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
extra_theme <-theme(axis.text.y = element_text( size=8),
                    axis.title = element_text(size=9),
                    plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),
                    legend.margin = margin(t = -12, r = 0, b = 0, l = 0),
                    axis.text.x = element_text(angle = 90 , size = 8, hjust = 1, vjust = 0.5, margin = margin(t = 1, r = 0, b = 0, l = 0)))

summary1$sites<-recode(summary1$sites, '5'='5 Clients', '10'='10 Clients')
colnames(summary1)<- c('iterations', 'sites', 'Algorithm', 'eigenvector', 'value')
conv<-ggplot(summary1[eigenvector %in% c('Eigenvector 1', 'Eigenvector 5', 'Eigenvector 10')], 
             aes(iterations, value, group=interaction(Algorithm, Algorithm), color=Algorithm))+
  geom_line(aes(linetype=Algorithm))+
  xlab('Iterations')+
  theme_classic()+
  ylab('Angles w.r.t reference')+
  facet_grid_sc(rows=vars(sites), cols = vars(eigenvector), scales = list(x = scales_x))+
  scale_x_continuous(breaks=breaks_fun, limits = c(0, NA))+
  theme(legend.position = 'bottom',
        strip.text = element_text(size=8), panel.border = element_rect(fill=NA), 
        panel.spacing = unit(0.05, 'cm'),strip.background = element_blank() )+
  scale_color_viridis_d(option = 'A', end=0.7)+extra_theme
conv

# Convert the plot to a grob
text = 'MNIST'
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

ggsave(conv,file='/home/anne/Documents/manuscripts/vertical-pattern-recognition/figures/convergence.pdf', height=6, units = 'cm', width = 15 )




extra_theme <-theme(axis.text.y = element_text( size=8),
                    axis.title = element_text(size=9),
                    plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),
                    legend.margin = margin(t = -12, r = 0, b = 0, l = 0),
                    axis.text.x = element_text(angle = 90 , size = 8, hjust = 1, vjust = 0.5, margin = margin(t = 1, r = 0, b = -12, l = 0)))

time <- fread('~/Documents/featurecloud/pca/approximative-vertical/results/time.log')
time$V1<-recode(time$V1, 'approximative-init'='AI-FULL', 'random-init'='RI-FULL', 'randomized-no-approx'='RI-RAND', 'randomized-approx-projected'='AI-RAND', 'guo'='GUO')
time<- time[V1 %in%  c('AI-FULL', 'RI-FULL', 'AI-RAND', 'RI-RAND', 'GUO')]
pl.runtime<-ggplot(time, aes(V1, V4, color=V1))+geom_boxplot()+ylab('Runtime[s]')+xlab('')+
  theme_classic()+scale_color_viridis_d(option = 'A', end=0.7)+
  extra_theme+
  guides(color=FALSE)

pl.runtime

eo<- fread('~/Documents/featurecloud/pca/approximative-vertical/results/matrix/summary.eo.tsv')
eog<-fread('~/Documents/featurecloud/pca/approximative-vertical/results/vector/summary.eo.tsv')
eo<-rbind(eo, eog)
eo$V4<-recode(eo$V4, 'approximative-init'='AI-FULL', 'random-init'='RI-FULL', 'randomized-no-approx'='RI-RAND', 'randomized-approx-projected'='AI-RAND', "guo"="GUO")
eo<- eo[V4 %in% c('AI-FULL', 'RI-FULL', 'AI-RAND', 'RI-RAND', 'GUO') ]
eo <- eo %>% group_by(V3,V4, V5, V6) %>% summarise(sum=sum(V2))
pl.eo<-ggplot(eo, aes(V4, sum, fill=interaction(V4,V3), color=V4))+geom_boxplot()+ylab('Runtime[s] (Matrix oper.)')+xlab('')+
  theme_classic()+
  extra_theme+
  guides(fill=FALSE, color=F)+
  scale_color_viridis_d(option = 'A', end=0.7)+
  scale_fill_manual(values=c(viridisLite::magma(5, end=0.7, alpha = 0.6), viridisLite::magma(5, end=0.7, alpha = 0.8)))
  
  
pl.eo

iterations<-rbind(as.data.table(iterations), as.data.table(max_it))
iterations$algorithm<-recode(iterations$algorithm, 'approximative-init'='AI-FULL', 'random-init'='RI-FULL', 'randomized-no-approx'='RI-RAND', 'randomized-approx-projected'='AI-RAND', 'guo'='GUO')
iterations <- iterations[algorithm %in% c('AI-FULL', 'RI-FULL', 'AI-RAND', 'RI-RAND', 'GUO') ]
pl.iter<-ggplot(iterations, aes(algorithm, iterations, fill=interaction(algorithm, sites),color=algorithm))+
  geom_boxplot()+
  xlab('')+
  theme_classic()+
  ylab('Iterations to nconvergence')+
  extra_theme +
  guides(fill=FALSE, color=F)+scale_color_viridis_d(option = 'A', end=0.7)+
  scale_fill_manual(values=c(viridisLite::magma(5, end=0.7, alpha = 0.6), viridisLite::magma(5, end=0.7, alpha = 0.8)))
pl.iter




transmission<-fread('~/Documents/featurecloud/pca/approximative-vertical/results/matrix/summary.transmission.tsv')
transmission_guo<-fread('~/Documents/featurecloud/pca/approximative-vertical/results/vector/summary.transmission.tsv')
transmission<-rbind(transmission_guo, transmission)
transmission<-transmission[V2 %in% c("H_local=CS", "H_global=SC")]
colnames(transmission) <- c('index', 'key', 'iteration', 'site', 'size', 'kA', 'sites', 'algorithm', 'maxit', 'filename')
transmission<-transmission[site==1]
volume<-transmission %>% select(size, sites, algorithm, filename)%>% group_by(sites, algorithm, filename)%>% summarise(total_volume=(sum(size)*4)/1000000000)

volume<-as.data.table(volume)
volume$algorithm<-recode(volume$algorithm, 'approximative-init'='AI-FULL', 'random-init'='RI-FULL', 'randomized-no-approx'='RI-RAND', 'randomized-approx-projected'='AI-RAND', 'guo'='GUO')
volume<- volume[algorithm %in% c( 'RI-FULL', 'AI-RAND', 'RI-RAND', 'AI-FULL', 'GUO')]
pl.volume<-ggplot(volume, aes(algorithm, total_volume, fill=interaction(algorithm, sites),color=algorithm))+geom_boxplot()+
  ylab('Total transmitted\n data [GB]')+
  xlab('')+
  theme_classic()+
  extra_theme+
  guides(fill=FALSE, color=F)+scale_color_viridis_d(option = 'A', end=0.7)+
  scale_fill_manual(values=c(viridisLite::magma(5, end=0.7, alpha = 0.6), viridisLite::magma(5, end=0.7, alpha = 0.8)))
pl.volume
############FAKE ALPHA LEGEND
pl.volume.legend<-ggplot(volume, aes(algorithm, total_volume, fill=algorithm, alpha=as.factor(sites)))+geom_boxplot()+
  theme_classic()+ylab('Total transmitted data [GB]')+xlab('')+
  scale_alpha_manual(values=c(0.2,0.7))+
  guides(fill=F,alpha=guide_legend('#Sites', override.aes=list(fill=hcl(c(15,195),100,0,alpha=c(0.2,0.7)))))
pl.volume.legend


legend_b <- get_legend(
  pl.volume.legend +theme(legend.position = 'bottom',
                          legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0),
                          legend.text = element_text(size = 8),
                          legend.title = element_text(size = 8))
)
text = 'MNIST'
size=10
text.grob = textGrob(text, rot = -90,
                     gp = gpar(fontsize = size))
gp<-ggarrange(pl.iter, pl.eo, pl.volume, ncol=3)
gp<-plot_grid(gp, text.grob, rel_widths = c(1,0.05))
gp<-plot_grid(gp, legend_b, ncol=1, rel_heights = c(1,0.1))
# Convert the plot to a grob
gp
ggsave(gp,file='/home/anne/Documents/manuscripts/vertical-pattern-recognition/figures/benchmark.pdf', height=6, units = 'cm', width = 15)

# 
# mev<-fread('~/Documents/featurecloud/pca/approximative-vertical/results/matrix//summary.mev.u.tsv')
# colnames(mev)<-c('iterations', 'mev', 'sites', 'algorithm', 'maxit', 'filename')
# mev <- mev %>% group_by(iterations, sites, algorithm) %>% summarise(mean_mev = mean(mev))
# mev<-as.data.table(mev)
# mev$algorithm<-recode(mev$algorithm, 'approximative-init'='AI-FULL', 'random-init'='RI-FULL', 'randomized-no-approx'='RI-RAND', 'randomized-approx-projected'='AI-RAND', 'guo'='GUO')
# mev<- mev[algorithm %in% c( 'RI-FULL', 'AI-RAND', 'RI-RAND', 'AI-FULL') ]
# ggplot(mev, aes(iterations,mean_mev, color=algorithm))+geom_line()+facet_wrap(~sites)+
#   xlab('MEV')+theme_classic()+ylab('iterations')+xlim(0,100)+scale_color_viridis_d(option = 'A', end=0.7)
# 
# mev<-fread('~/Documents/featurecloud/pca/approximative-vertical/results/vector/summary.mev.u.tsv')
# colnames(mev)<-c('vector','iterations', 'mev', 'sites', 'algorithm', 'maxit', 'filename')
# mev <- mev %>% group_by(iterations, sites, algorithm) %>% summarise(mean_mev = mean(mev))
# mev<-as.data.table(mev)
# mev$algorithm<-recode(mev$algorithm, 'approximative-init'='AI-FULL', 'random-init'='RI-FULL', 'randomized-no-approx'='RI-RAND', 'randomized-approx-projected'='AI-RAND', 'guo'='GUO')
# ggplot(mev, aes(iterations,mean_mev, color=algorithm))+geom_line()+facet_wrap(~sites)+
#   xlab('MEV')+theme_classic()+ylab('iterations')+xlim(0,100)+scale_color_viridis_d(option = 'A', end=0.7)
# # 