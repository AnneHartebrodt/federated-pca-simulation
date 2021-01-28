suppressMessages(require(data.table))
suppressMessages(require(ggplot2))
suppressMessages(require(tidyr))
suppressMessages(require(cowplot))
suppressMessages(require(optparse))
suppressMessages(require(dplyr))
suppressMessages(require(stringr))
suppressMessages(require(optparse))

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
print(infile)
#infile<-'/home/anne/Documents/featurecloud/pca/vertical-pca/results/MMRF-COMMPASS/eigenval.tsv'
#column<-'eigenval'


data<-fread(infile)
data[, 9]<-as.numeric(unlist(data[, 9]))

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