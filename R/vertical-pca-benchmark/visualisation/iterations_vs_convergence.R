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
  data<-fread('/home/anne/Documents/featurecloud/pca/vertical-pca/results/results_for_david/chr2.tsv')
  d <- data %>% pivot_longer(-c(iterations))
  d<-as.data.table(d)
  d<-d[!is.na(value)]
  # select the correct configuration
  selection<-c("matrix_5_power_central_qr_10", "matrix_5_power_federated_qr_10","vector_5_gradient_central_qr_10")
  #selection3<-c("matrix_3_power_central_qr", "matrix_3_power_federated_qr", "vector_3_power_central_qr", "vector_3_power_federated_qr")
  
  
  #selection<-c("vector_2_gradient_central_qr")
  #selection3<-c("matrix_3_power_central_qr", "matrix_3_power_federated_qr", "vector_3_power_central_qr", "vector_3_power_federated_qr")
  
  # make the plot
  angles.plot<-ggplot(d, aes(iterations, value, col=as.factor(name)))+
    geom_line()+
    my_theme+ylab('Mean angle [degree]')+   
    xlab('#Iterations')+
   # scale_color_manual('Configuration', values = palette_div)+
    theme(axis.line=element_line(),
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          legend.position = c(1, 1),
          legend.justification = c("right", "top"),
          legend.box.just = "right",
          legend.margin = margin(0.25, 0.25, 0.25, 0.25),
          legend.title = element_text(size=10))#+
    #guides(color=guide_legend(keyheight = 0.5, title = element_text('Configuration', size = 8)))
  
  angles.plot
  
