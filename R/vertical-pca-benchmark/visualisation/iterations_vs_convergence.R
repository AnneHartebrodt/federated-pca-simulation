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
  data<-fread('/home/anne/Documents/featurecloud/pca/vertical-pca/results-new-tests/mnist/long.dummy.angles.u.summary.tsv')
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
  data$rank <-as.factor(data$rank)
  data[, facet_title:=paste('Eigenvector', rank)]
  data$facet_title <- as.ordered(data$facet_title)
  levels(data$facet_title)<- paste('Eigenvector', unique(sort(data$rank)))
  data[, name_qr := paste0(matrix, eigenvector_update, qr_method, orthonormalisation_skip)]
  data$rank <-as.factor(data$rank)
  data$name_qr<-as.factor(data$name_qr)
  data$name_qr<- recode(data$name_qr, matrixpowerfederated_qr1='FED-GS', matrixpowerfederated_qr100='NO-GS', vectorgradientcentral_qr1='GUO')
  data <- as.data.table(data)
  angles.plot<-ggplot(data, aes(iterations, mean_value, col=as.factor(name_qr)))+
    geom_line(size=0.75)+
    theme_bw()+ylab('angle w.r.t reference')+   
    xlab('iterations')+
    scale_color_manual(values = viridis(4)[1:3])+
    theme(axis.line=element_line(),
          strip.background = element_blank(),
          legend.justification = c("right", "top"),
          legend.box.just = "right",
          legend.title = element_blank(),
          panel.grid = element_blank())+
    facet_wrap('facet_title', nrow=2, scales = 'free_x')+
    guides(color=guide_legend(keyheight = 0.6, title = element_text('Configuration', size = 8)))
  angles.plot
  
