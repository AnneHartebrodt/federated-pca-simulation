require(data.table)
require(ggplot2)
require(tidyr)
require(cowplot)
require(optparse)
require(dplyr)
require(stringr)
require(optparse)

option_list = list(
  make_option(c("-b", "--basedir"), action="store", default=NA, type='character',
              help="base directory"),
  make_option(c("-s", "--suffix"), action="store", default=NA, type='character',
              help="file ending of files to be combined"),
  make_option(c("-c", "--colname"), action="store", default=NA, type='character',
              help="column name"),
  make_option(c("-o", "--outfile"), action="store", default=NA, type='character',
              help="output file")
)
opt = parse_args(OptionParser(option_list=option_list))

read_all_files<-function(base_directory, suffix){
  data_list<-list()
  counter<-1
  for(o in data_orientation){
    for(v in vector_or_subspace)
      for(s in nr_sites){
        for(e in eigenvector_update){
          for(q in qr_method){
            
            dir.name <- file.path(base_directory, o, v, s, e, q)
            if(!file.exists(dir.name)){
              next
            }
            myfiles <-  list.files(dir.name, recursive = F)
            for (f in which(str_detect(myfiles, paste0(suffix, '$')))) {
              #print(file.path(dir.name, myfiles[f]))
              data<- fread(file.path(dir.name, myfiles[f]),
                           sep = '\t',
                           header = F)
              data$orientation<-o
              data$matrix<-v
              data$sites<-s
              data$eigenvector_update<-e
              data$qr_method<-q
              data$filename<-file.path(dir.name, myfiles[f])
              head(data)
              data_list[[file.path(dir.name, myfiles[f])]]<-data
            }
          }
        }
      }
  }
  return(data_list)
}

basedir<-opt$b
suffix <- opt$s
column_name <- opt$c
outfile<-opt$o

data_orientation<-c('vertical', 'horizontal')
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
                        'qr_method', 'filename')  %>%  
                pivot_longer(-c('iterations',
                                   'orientation', 'matrix', 'sites', 'eigenvector_update',
                                   'qr_method', 'filename'), values_to = column_name, names_to = 'rank')


# all data frames need to have the same shape
vector_or_subspace<-c('vector')
mf_v<-read_all_files(basedir, 'u')
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
fwrite(data, file = file.path(basedir, outfile), sep='\t')

