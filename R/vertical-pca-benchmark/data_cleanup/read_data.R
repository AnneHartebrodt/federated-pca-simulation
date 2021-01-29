suppressMessages(require(data.table))
suppressMessages(require(ggplot2))
suppressMessages(require(tidyr))
suppressMessages(require(cowplot))
suppressMessages(require(optparse))
suppressMessages(require(dplyr))
suppressMessages(require(stringr))
suppressMessages(require(optparse))

option_list = list(
  make_option(c("-b", "--basedir"), action="store", default=NA, type='character',
              help="base directory"),
  make_option(c("-s", "--suffix"), action="store", default=NA, type='character',
              help="file ending of files to be combined"),
  make_option(c("-c", "--colname"), action="store", default=NA, type='character',
              help="column name"),
  make_option(c("-o", "--outfile"), action="store", default=NA, type='character',
              help="output file"),
  make_option(c("-d", "--scriptdir"), action="store", default=NA, type='character',
              help="output file"),
  make_option(c("-a", "--scalability"), action="store_true", default=FALSE, type='logical',
              help="output file"),
  make_option(c("-m", "--matrix_only"), action="store_true", default=FALSE, type='logical',
              help="output file")

)
opt = parse_args(OptionParser(option_list=option_list))

scriptdir<-opt$scriptdir
source(file.path(scriptdir, "/R/vertical-pca-benchmark/data_cleanup/library.R"))

basedir<-opt$basedir
suffix <- opt$suffix
column_name <- opt$colname
outfile<-opt$outfile
scalabililty<-opt$scalability
print(scalabililty)

print(basedir)
#source("/home/anne/Documents/featurecloud/pca/federated_dp_pca/R/vertical-pca-benchmark/data_cleanup/library.R")
#basedir<-'/home/anne/Documents/featurecloud/pca/vertical-pca/results/mnist'
#suffix<-'angles.u'
#column_name<-'angle'

if(scalabililty){
  data_orientation<-c('0.1','0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9')
}else{
data_orientation<-c('vertical', 'horizontal')
}
vector_or_subspace<-c('matrix')
nr_sites<- c(2,3,5,10)
eigenvector_update<-c('power', 'gradient')
qr_method<-c('central_qr', 'federated_qr')



mf_m<-read_all_files(basedir, suffix)
mf_m<-rbindlist(mf_m)
# rename
colnames(mf_m)<-c( 'iterations', 1:10, 'empty', 
                  'orientation', 'matrix', 'sites', 'eigenvector_update',
                  'qr_method', 'filename')

# remove the one empty column and transform into long format
# use characters as column names!
mf_m <- mf_m %>% select('iterations', paste0(1:10, ''), 
                        'orientation', 'matrix', 'sites', 'eigenvector_update',
                        'qr_method', 'filename')%>%  
                pivot_longer(-c('iterations',
                                   'orientation', 'matrix', 'sites', 'eigenvector_update',
                                   'qr_method', 'filename'), values_to = column_name, names_to = 'rank')

if (!opt$matrix_only){
# all data frames need to have the same shape
vector_or_subspace<-c('vector')
mf_v<-read_all_files(basedir, suffix)
mf_v<-rbindlist(mf_v)

colnames(mf_v)
colnames(mf_v)<-c('rank', 'iterations', column_name, 'empty', 
                  'orientation', 'matrix', 'sites', 'eigenvector_update',
                  'qr_method', 'filename')
# remove the empty column.
mf_v<- mf_v %>% select('rank', 'iterations',
                       'orientation', 'matrix', 'sites', 'eigenvector_update',
                       'qr_method', 'filename', column_name)

#combine data
data<-rbind(mf_m, mf_v)
}else{
data<-mf_m
}
fwrite(data, file = file.path(basedir, outfile), sep='\t', quote = F)

