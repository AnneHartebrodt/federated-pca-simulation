require(data.table)
require(ggplot2)
require(tidyr)
require(cowplot)
require(optparse)
require(dplyr)
require(stringr)
require(optparse)

option_list = list(
  make_option(c("-f", "--infile"), action="store", default=NA, type='character',
              help="infile"),
  make_option(c("-o", "--outfile"), action="store", default=NA, type='character',
              help="outfile"),
  make_option(c("-c", "--column"), action="store", default=NA, type='character',
              help="columnname")
  
  
)
opt = parse_args(OptionParser(option_list=option_list))

infile<-opt$infile
outfile<-opt$outfile
column<-opt$column

infile<-'/home/anne/Documents/featurecloud/pca/vertical-pca/results/mnist_january/angles.u.tsv'
column<-'angle'

data<-fread(infile)
summary<-data %>% group_by(iterations,orientation, matrix, sites, eigenvector_update, qr_method, rank) %>%
        summarise(mean_value = mean(get(column)))
summary<-as.data.table(summary)

fwrite(summary, outfile, sep='\t')

## tranform to wide format
## split into vertical and horizontal
wide.vertical<-summary[orientation=='vertical'] %>% pivot_wider(id_cols = c(iterations, rank), names_from = c(matrix, sites, eigenvector_update, qr_method), values_from = mean_value)
wide.vertical<-as.data.table(wide.vertical)
fwrite(wide.vertical, paste0('wide.vertical.', outfile), sep='\t')


wide.hori<-summary[orientation=='horizontal'] %>% pivot_wider(id_cols = c(iterations, rank), names_from = c(matrix, sites, eigenvector_update, qr_method), values_from = mean_value)
wide.hori<-as.data.table(wide.hori)
fwrite(wide.hori, paste0('wide.horizontal.', outfile), sep='\t')