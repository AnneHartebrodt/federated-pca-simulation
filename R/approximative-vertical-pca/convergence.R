require(data.table)
require(ggplot2)
require(tidyr)
require(optparse)
require(dplyr)
require(R.utils)

require(scales)
library(ggforce)
require(stringr)
require(viridis)
require(gridExtra)
library(ggpubr)
require(facetscales)

option_list = list(
  make_option(c("-f", "--infiles"), action="store", default=NA, type='character',
              help="infiles, comma separated list"),
  make_option(c("-o", "--outfile"), action="store", default=NA, type='character',
              help="outfile"),
  make_option(c("-n", "--names"), action="store", default=NA, type='character',
              help="test names, comma separated list")
  
  
)
opt = parse_args(OptionParser(option_list=option_list))

infiles<-opt$infiles
infiles<-'/home/anne/Documents/featurecloud/pca/vertical-pca/results-new-tests/1000g/chr2.tsv,/home/anne/Documents/featurecloud/pca/vertical-pca/results-new-tests/mnist/long.dummy.angles.u.summary.tsv'
#infiles<-'/home/anne/Documents/featurecloud/pca/approximative-vertical/results-new-tests/mnist/long.dummy.angles.u.summary.tsv'
outfile<-opt$outfile
outfile <- '/home/anne/Documents/manuscripts/approximative-vertical-pca/angles_all.pdf'
names<-opt$names
names<-'MNIST'


infiles<- str_split(infiles, ',')[[1]]
names <- str_split(names, ',')[[1]]
# read the data and transform into long
assertthat::are_equal(length(infiles), length(names))

files <-list()
for(i in 1:length(infiles)){
  chr1<-fread(infiles[[i]])
  chr1$test<- names[[i]]
  files[[i]]<- chr1
}



data<- rbindlist(files)

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



sub <- data[facet_title %in% paste('Eigenvector', c(1,5,10))& test =='MNIST']
mnist.plot<-ggplot(sub, aes(iterations, mean_value, col=as.factor(name_qr)))+
  geom_line(size=0.75)+
  theme_bw()+ylab('')+   
  xlab('')+
  scale_color_manual(values = viridis(4)[1:3])+
  theme(axis.line=element_line(),
        strip.background = element_blank(),
        axis.title = element_blank(),
        legend.position = 'none',
        legend.box.just = "right",
        legend.title = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(size=6),
        plot.margin=unit(c(-0.1,0.1,0.1,0.1), "cm"),
        legend.text = element_text(size = 6))+
  facet_grid_sc(rows=vars(test), cols = vars(facet_title), scales = list(x = scales_x))+
  guides(color=guide_legend(keyheight = 0.6, title = element_text('Configuration', size = 8)))
mnist.plot

#sub.plot<- ggarrange(chr2.plot, mnist.plot, nrow=1)
full.plot<- ggarrange(chr1.plot, chr2.plot,mnist.plot, nrow=3, common.legend=TRUE)


# Annotate the figure by adding a common labels
a<- annotate_figure(full.plot,
                    bottom = text_grob("iterations", size = 10),
                    left = text_grob("angle w.r.t reference", rot = 90, size = 10)
)
#full.plot
ggsave(a, filename = outfile, width=10, height=10, units = 'cm')



scales_x <- list(
  `Eigenvector 1` = scale_x_continuous(limits = c(0, 15)),
  #`Eigenvector 2` = scale_x_continuous(limits = c(0, 30)),
  `Eigenvector 5` = scale_x_continuous(limits = c(0, 150)),
  #`Eigenvector 7` = scale_x_continuous(limits = c(0, 300)),
  `Eigenvector 10` =scale_x_continuous(limits = c(0,1400))
  
)
chr1.plot<-ggplot(data[facet_title %in% paste('Eigenvector', c(1,5,10)) & test =='Chrom. 1'], aes(iterations, mean_value, col=as.factor(name_qr)))+
  geom_line(size=0.75)+
  theme_bw()+ylab('')+   
  xlab('')+
  scale_color_manual(values = viridis(4)[1:3])+
  theme(axis.line=element_line(),
        strip.background = element_blank(),
        legend.position = 'none',
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        axis.title = element_blank(),
        legend.title = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(size=7),
        plot.margin=unit(c(0.1,0.1,0.1,0.1), "cm"))+
  #facet_wrap('facet_title', nrow=1, scales ='free_x')+
  facet_grid_sc(rows=vars(test), cols = vars(facet_title), scales = list(x = scales_x))+
  scale_x_continuous(breaks=breaks_fun, limits = c(0, NA))+
  guides(color=guide_legend(keyheight = 0.6, title = element_text('Configuration', size = 8)))
chr1.plot

#ggsave(angles.plot, filename = '/home/anne/Documents/featurecloud/pca/vertical-pca/figures/angles_chr1.pdf', width=20, height=10, units = 'cm')
scales_x <- list(
  `Eigenvector 1` = scale_x_continuous(limits = c(0, 30)),
  `Eigenvector 5` = scale_x_continuous(limits = c(0, 400)),
  `Eigenvector 10` =scale_x_continuous(limits = c(0,900))
  
)

chr2.plot<-ggplot(data[facet_title %in% paste('Eigenvector', c(1,5,10))& test =='Chrom. 2'], aes(iterations, mean_value, col=as.factor(name_qr)))+
  geom_line(size=0.75)+
  theme_bw()+ylab('')+   
  xlab('')+
  scale_color_manual(values = viridis(4)[1:3])+
  theme(axis.line=element_line(),
        legend.position = 'none',
        axis.title = element_blank(),
        strip.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(size=7),
        plot.margin=unit(c(-0.1,0.1,0.1,0.1), "cm"))+
  #facet_wrap('facet_title', nrow=1, scales ='free_x')+
  facet_grid_sc(rows=vars(test), cols = vars(facet_title), scales = list(x = scales_x))+
  scale_x_continuous(breaks=breaks_fun, limits = c(0, NA))+
  guides(color=guide_legend(keyheight = 0.6, title = element_text('Configuration', size = 8)))
chr2.plot

#ggsave(angles.plot, filename = '/home/anne/Documents/featurecloud/pca/vertical-pca/figures/angles_chr2.pdf', width=20, height=10, units = 'cm')
scales_x <- list(
  `Eigenvector 1` = scale_x_continuous(limits = c(0, 30)),
  `Eigenvector 5` = scale_x_continuous(limits = c(0, 400)),
  `Eigenvector 10` =scale_x_continuous(limits = c(0,900))
  
)