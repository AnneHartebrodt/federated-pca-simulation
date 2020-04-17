require(data.table)
require(ggplot2)
require(cowplot)
require(optparse)

option_list = list(
  make_option(c("-r", "--resultfolder"), action="store", default=NA, type='character',
              help="annoation file in gtf format"),
  make_option(c("-s", "--study.sizes"), action="store", default=NA, type='character',
              help="annoation file in gtf format"),
  make_option(c("-p", "--plotdir"), action="store", default=NA, type='character',
              help="The directory for plot")
)
opt = parse_args(OptionParser(option_list=option_list))

plotdir <-opt$plotdir

eigenvector.path <-file.path(opt$plotdir, 'leave1out', 'eigenvectors')
dir.create(eigenvector.path, recursive = T)

eigenvalue.path <-file.path(opt$plotdir, 'leave1out', 'eigenvalues')
dir.create(eigenvalue.path, recursive = T)

make_angle_plot<-function(angles, dataset){
  p.angles<-ggplot(angles, aes(x = as.factor(eigenvector.rank), y = angle, fill=nr.dropped))+
    geom_boxplot()+
    xlab('Eigenvectors rank')+
    ylab('Angle (degree)')+
    ggtitle('Angles between leading eigenvectors of PCA between two adjacent datasets')+
    theme(plot.title = element_text(size = 25), 
          legend.key.size = unit(1.5, "cm"),
          legend.title = element_text(size = 20),
          strip.text = element_text(size = 25), axis.text = element_text(size=20),
          axis.title = element_text(size=25))+
    ggtitle(dataset)
  p.angles
}


make_eigen_boxplot<-function(filename, nr_dropped){
  dirs<-list.dirs('/home/anne/Documents/featurecloud/results/leave1out', recursive = F)
  summaries<-list()
  max.inds<-list()
  eigs<-list()
  namesd<-c('outlier_free_1', 'outlier_free_5', 'outlier_free_10', 'cheat_10')
  indxd <-c('1', '5', '10', '10')
  id <- c('1', '5', '10', 'O')
  for(d in dirs){
    for(i in 1:4){
      # load data and create individial plot
      name<-file.path(d, namesd[i], paste0('eigenvalues.tsv'))
      if(file.exists(name)){
        if(file.info(name)$size != 0){
    eig<-fread(name, fill= T)
    eig$rn <-seq(1, nrow(eig))
    eig$nr.dropped<-id[i]
    eigs[[name]]<-eig
        }
      }
    }
  }
  
  eig<-rbindlist(eigs, fill = T)
  eig<-melt(eig, id.vars = c('V1', "V2", "rn", "nr.dropped"), variable.name = 'eigenvalue.rank', value.name = 'eigenvalue')
  eig<-eig[, .(V1, rn, nr.dropped, eigenvalue.rank, eigenvalue)]
  eig$eigenvalue.rank<-sapply(eig$eigenvalue.rank, function(x) gsub('V', '', x))
  eig$eigenvalue.rank<-as.factor(as.numeric(eig$eigenvalue.rank)-2)
  colnames(eig)<-c('study.id' , 'j','nr.dropped', 'eigenvalue.rank', 'eigenvalue')
  eig$nr.dropped<-as.factor(eig$nr.dropped)
  eig<-eig[!is.na(eig$eigenvalue)]
  for(dataset in unique(eig$study.id)){
  p.eigenvalues<-ggplot(eig[study.id==dataset], aes(x = as.factor(eigenvalue.rank), y = eigenvalue, fill=nr.dropped))+
    geom_boxplot()+
    xlab('Eigenvalue rank')+
    ylab('Eigenvalue') + 
    theme(plot.title = element_text(size = 25), 
                              legend.key.size = unit(1.5, "cm"),
                              legend.title = element_text(size = 20),
                              strip.text = element_text(size = 25), axis.text = element_text(size=20),
                              axis.title = element_text(size=25))+ggtitle(dataset)
    ggsave(p.eigenvalues, filename = file.path(eigenvalue.path,  paste0(dataset, '.pdf')), width = 10,
         height = 10, dpi = 'print')
  }
}
  
pp1<-make_eigen_boxplot(filenames, drops)

study_sizes<-fread('/home/anne/Documents/featurecloud/results/usability_study/ordered_decr.tsv', header = F)


dirs<-list.dirs(opt$resultfolder, recursive = F)
summaries<-list()
max.inds<-list()

namesd<-c('outlier_free_1', 'outlier_free_5', 'outlier_free_10', 'cheat_10')
indxd <-c('1', '5', '10', '10')

id <- c('1', '5', '10', 'O')
for(d in dirs){
  if(!(basename(d) %in% study_sizes$V1)){
  next
  }
  alist<-list()
  for(i in 1:4){
  # load data and create individial plot
  name<-file.path(d, namesd[i], paste0('angles_dropout', indxd[i], '.tsv'))
  if(file.exists(name)){
  if(file.info(name)$size != 0){
  angles<-fread(name)
  colnames(angles)<-c('i', 'j', 'eigenvector.rank', 'angle')
  angles$nr.dropped<-id[i]
  alist[[i]]<-angles
  
  summary<-angles[, list(min=min(angle), median = median(angle), max = max(angle), max.ind = which.max(angle)) , by = eigenvector.rank]
  summary$nr.dropped<-id[i]
  summary$study.id<-basename(d)
  summaries[[name]]<-summary
  }
  }
  }
  if(length(alist)>0){
  alist<-rbindlist(alist, fill=T)
  p.eigenvalues<-make_angle_plot(alist,basename(d))
  ggsave(p.eigenvalues, filename = file.path(eigenvector.path,  paste0(basename(d), '.pdf')), width = 20,
         height = 20, dpi = 'print')
  }
}

summaries<-rbindlist(summaries)

summaries<-melt(summaries, id.vars = c('eigenvector.rank', 'nr.dropped', 'study.id'), variable.factor = T)
summaries$nr.dropped<-as.factor(summaries$nr.dropped)
summaries$nr.dropped<-ordered(summaries$nr.dropped, levels=c('1', '5', '10', 'O'))

study.sizes<-fread(opt$study.sizes)
colnames(study.sizes)<-c('study.id', 'nr.samples', 'database')
study.sizes[study.id=='BEATAML1.0-COHORT', 1]<-'BEATAML1'
summaries<-merge(summaries, study.sizes, by.x = 'study.id', by.y = 'study.id')

summaries$eigenvector.rank<-summaries$eigenvector.rank+1

gp<-ggplot(summaries[eigenvector.rank %in% c(1,3,6,9) & variable!='max.ind'], aes(nr.dropped, value,fill = variable))+
 geom_boxplot(outlier.shape = NA)+
  scale_fill_brewer(palette = 'Blues')+
  geom_point(position=position_jitterdodge(dodge.width=0.8), size=1)+
  facet_wrap(~eigenvector.rank, ncol = 2)+
  xlab('Number of records dropped')+ylab('Angle (degree)')+
  theme(plot.title = element_text(size = 35, hjust = 0.5), 
        legend.key.size = unit(1, "cm"),
        legend.title = element_text(size = 25),
        strip.text = element_text(size = 35), 
        axis.text = element_text(size=25),
        axis.title = element_text(size=35),
        plot.subtitle = element_text(size = 25, hjust = 0.5),
        legend.text=element_text(size=20))+
  labs(fill='Statistic', alpha='#Samples')
gp  
#ggsave(gp, file='/home/anne/Documents/paper_fed_PCA/figures/angles_dropout.pdf', width = 20, height = 15, dpi = 'print')
ggsave(gp, file=file.path(plotdir, 'angles_dropout_summary.eps'), width = 15, height = 10)
ggsave(gp, file=file.path(plotdir, 'angles_dropout_summary.png'), width = 15, height = 10, dpi = 350)
