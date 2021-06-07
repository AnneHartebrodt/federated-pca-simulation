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
  require(viridis)
  require(gridExtra)
  library(ggpubr)
  
  # make custom palettes
  palette_div<-c('#062865', '#2f497d', '#4c6d96', '#6793af', '#81bac8', '#ffc4b3', '#f68888', '#d75161', '#ab203f', '#720022')
  palette_seq<-c('#062865', '#203a72', '#324d80', '#43618d', '#52759b', '#618aa9', '#70a0b7', '#7eb6c5', '#8dcdd4', '#9ce4e2')
  
  # create a minimal theme
  my_theme <-
    theme_classic() + theme(
      axis.title = element_text(size = 12),
      legend.title = element_text(size = 12),
      plot.title = element_text(size = 20, hjust = 0.5),
      axis.text = element_text(size = 8),
      legend.text = element_text(size = 8),
      plot.subtitle = element_text(size = 10, hjust = 0.5),
      legend.key.size = unit(1.5, 'lines'))
  
  
  # read the data and transform into long
  chr1<-fread('/home/anne/Documents/featurecloud/pca/vertical-pca/results-new-tests/1000g/chr1.tsv')
  chr1$test<- 'Chromosome 1'
  chr2<- fread('/home/anne/Documents/featurecloud/pca/vertical-pca/results-new-tests/1000g/chr2.tsv')
  chr2$test <- 'Chromosome 2'
  mnist<- fread('/home/anne/Documents/featurecloud/pca/vertical-pca/results-new-tests/mnist/long.dummy.angles.u.summary.tsv')
  mnist$test<- 'MNIST'
  data<- rbind(chr1, chr2, mnist)
  #d <- data %>% pivot_longer(-c(iterations))
  #d<-as.data.table(d)
  #d<-d[!is.na(value)]
  # select the correct configuration
  #selection<-c("matrix_5_power_central_qr_10", "matrix_5_power_federated_qr_10","vector_5_gradient_central_qr_10")
  #selection3<-c("matrix_3_power_central_qr", "matrix_3_power_federated_qr", "vector_3_power_central_qr", "vector_3_power_federated_qr")
  
  
  #selection<-c("vector_2_gradient_central_qr")
  #selection3<-c("matrix_3_power_central_qr", "matrix_3_power_federated_qr", "vector_3_power_central_qr", "vector_3_power_federated_qr")
  #selection<-c("matrix_5_power_central_qr_1", "matrix_5_power_federated_qr_1","vector_5_gradient_central_qr_1")
  # make the plot
  data$rank <-as.ordered(data$rank)
  data[, facet_title:=ordered(paste('Eigenvector', rank), levels = paste('Eigenvector', unique(rank)))]
  #data$facet_title <- ordered(data$facet_title, levels=levels(data$rank))
  data[, name_qr := paste0(matrix, eigenvector_update, qr_method, orthonormalisation_skip)]
  data$rank <-as.factor(data$rank)
  data$name_qr<-as.factor(data$name_qr)
  data$name_qr<- recode(data$name_qr, matrixpowerfederated_qr1='FED-GS', matrixpowerfederated_qr100='NO-GS', vectorgradientcentral_qr1='GUO')
  data <- as.data.table(data)
  
  data[, breaks:=(data$iterations - (data$iterations %% 50))]
  
  
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
  
  scales_x <- list(
    `Eigenvector 1` = scale_x_continuous(limits = c(0, 15)),
    `Eigenvector 2` = scale_x_continuous(limits = c(0, 30)),
    `Eigenvector 5` = scale_x_continuous(limits = c(0, 150)),
    `Eigenvector 7` = scale_x_continuous(limits = c(0, 300)),
    `Eigenvector 10` =scale_x_continuous(limits = c(0,1400))
    
  )
  chr1.plot<-ggplot(data[facet_title %in% paste('Eigenvector', c(1,2,5,7,10)) & test =='Chromosome 1'], aes(iterations, mean_value, col=as.factor(name_qr)))+
    geom_line(size=0.75)+
    theme_bw()+ylab('')+   
    xlab('')+
    scale_color_manual(values = viridis(4)[1:3])+
    theme(axis.line=element_line(),
          strip.background = element_blank(),
          legend.position = 'none',
          legend.justification = c("right", "top"),
          legend.box.just = "right",
          
          legend.title = element_blank(),
          panel.grid = element_blank(),
          axis.text = element_text(size=7),
          plot.margin=unit(c(0.1,0.1,-0.2,0.1), "cm"))+
    #facet_wrap('facet_title', nrow=1, scales ='free_x')+
    facet_grid_sc(rows=vars(test), cols = vars(facet_title), scales = list(x = scales_x))+
    scale_x_continuous(breaks=breaks_fun, limits = c(0, NA))+
    guides(color=guide_legend(keyheight = 0.6, title = element_text('Configuration', size = 8)))
 chr1.plot
  
  #ggsave(angles.plot, filename = '/home/anne/Documents/featurecloud/pca/vertical-pca/figures/angles_chr1.pdf', width=20, height=10, units = 'cm')

  
  chr2.plot<-ggplot(data[facet_title %in% paste('Eigenvector', c(1,2,5,7,10))& test =='Chromosome 2'], aes(iterations, mean_value, col=as.factor(name_qr)))+
    geom_line(size=0.75)+
    theme_bw()+ylab('angle w.r.t reference')+   
    xlab('')+
    scale_color_manual(values = viridis(4)[1:3])+
    theme(axis.line=element_line(),
          legend.position = 'none',
          strip.background = element_blank(),
          panel.grid = element_blank(),
          axis.text = element_text(size=7),
          plot.margin=unit(c(-0.1,0.1,-0.1,0.1), "cm"))+
    #facet_wrap('facet_title', nrow=1, scales ='free_x')+
    facet_grid_sc(rows=vars(test), cols = vars(facet_title), scales = list(x = scales_x))+
    scale_x_continuous(breaks=breaks_fun, limits = c(0, NA))+
    guides(color=guide_legend(keyheight = 0.6, title = element_text('Configuration', size = 8)))
  chr2.plot
  
  #ggsave(angles.plot, filename = '/home/anne/Documents/featurecloud/pca/vertical-pca/figures/angles_chr2.pdf', width=20, height=10, units = 'cm')
  scales_x <- list(
    `Eigenvector 1` = scale_x_continuous(limits = c(0, 30)),
    `Eigenvector 2` = scale_x_continuous(limits = c(0, 100)),
    `Eigenvector 5` = scale_x_continuous(limits = c(0, 400)),
    `Eigenvector 7` = scale_x_continuous(limits = c(0, 500)),
    `Eigenvector 10` =scale_x_continuous(limits = c(0,900))
    
  )
  
  sub <- data[facet_title %in% paste('Eigenvector', c(1,2,5,7,10))& test =='MNIST']
  mnist.plot<-ggplot(sub, aes(iterations, mean_value, col=as.factor(name_qr)))+
    geom_line(size=0.75)+
    theme_bw()+ylab('')+   
    xlab('iterations')+
    scale_color_manual(values = viridis(4)[1:3])+
    theme(axis.line=element_line(),
          strip.background = element_blank(),
          legend.position = 'bottom',
          legend.box.just = "right",
          legend.title = element_blank(),
          panel.grid = element_blank(),
          axis.text = element_text(size=7),
          plot.margin=unit(c(-0.2,0.1,0.1,0.1), "cm"))+
    facet_grid_sc(rows=vars(test), cols = vars(facet_title), scales = list(x = scales_x))+
    guides(color=guide_legend(keyheight = 0.6, title = element_text('Configuration', size = 8)))
  mnist.plot
  
  full.plot<- ggarrange(chr1.plot, chr2.plot, mnist.plot, nrow=3, common.legend=TRUE)
  full.plot
  ggsave(full.plot, filename = '/home/anne/Documents/featurecloud/pca/vertical-pca/figures/angles_all.pdf', width=20, height=15, units = 'cm')
  