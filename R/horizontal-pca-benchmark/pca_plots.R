require(data.table)
require(ggplot2)
require(tidyr)
require(stringr)
source('/home/anne/Documents/featurecloud/pca/federated_dp_pca/R/horizontal-pca-benchmark/library.R')

type = "Thyroid_gland"

#read_and_pca<-function(type){
  data<-fread(paste0('/home/anne/Documents/featurecloud/data/tcga/cancer_type/',type,'/sites/data_all/',type,'.tsv'))
  names<-fread(paste0('/home/anne/Documents/featurecloud/data/tcga/cancer_type/',type,'/sites/data_all/',type,'.names.tsv'),
               header=F, sep='\t')
  
  names[,c('type', 'tss', 'id') := tstrsplit(V1, '-')]
  pca <- prcomp(data, scale. = F, center = F)
  
  x <- as.data.table(pca$x)[, 1:10]
  x$id<-names$id
  x$tss<-names$tss
  x$name<-names$V1
 # return(x)
#}
#https://www.r-bloggers.com/detecting-outlier-samples-in-pca/
robust.outlier<-function(x, m) which( (abs(x - median(x)) / mad(x)) > m)

ou<-lapply(3:6, function (m) data.table(ol=unique(as.numeric(unlist(apply(x[,1:10], 2, function(a) robust.outlier(a, m))))),
                                           factor=m))

ou<-rbindlist(ou)
ou[, count:=.N, by=.(ol)]

ou.names<-x[ou]$name
ou<-data.table(name = ou.names, method = 'mad')
ol$method<-'angle'
ol1<-right_join(ol, ou, by= c('algorithm'='name'))
ol1$col<- 'both'
ol1$col<- ifelse(!is.na(ol1$method.x) & is.na(ol1$method.y) , 'angle', 'both')
ol1$col<- ifelse(is.na(ol1$method.x) & !is.na(ol1$method.y) , 'mad', 'both')

g<-ggplot(x, aes(PC1, PC2, col=as.factor(tss), shape=tss))+
  geom_point(aes(shape=tss))+
  theme_bw()+
  scale_color_manual(name='TSS',values = c(palette_div, palette_div,palette_div, palette_div))+
  scale_shape_manual(name = 'TSS',values = c(1:40))+
  guides(color=guide_legend(ncol=3))+theme(legend.key.size = unit(5, 'mm'), 
                                           legend.title = element_text(size = 10))+
  geom_label_repel(data = subset(x, name %in% ol1$algorithm),
                   aes(label=id, fill=ol1$col), size=2)+
  scale_fill_manual("Outlier\nidentification", values = c('#fde600','#fffac5','#fcf27e' ))


g

pal<-c('#d4b8aa','#d6c4a8', '#dca8a3' )