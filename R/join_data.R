require(data.table)
require(optparse)

option_list = list(
  make_option(c("-d", "--dir"), action="store", default=NA, type='character',
              help="coding genes output file"),
  make_option(c("-o", "--out"), action="store", default=NA, type='character',
              help="outdir")
)
opt = parse_args(OptionParser(option_list=option_list))


#
opt$d<-'/home/anne/Documents/featurecloud/data/tcga/data_clean'
opt$o<-'/home/anne/Documents/featurecloud/data/tcga/data_fake'
dirs<- list.dirs(opt$d, recursive = F)
for(d in 1: length(dirs)){
  nr_lines<-nrow(fread(file.path(dirs[d],'coding_only.tsv'), header = T))
  for(dd in (d+1):length(dirs)){
    nr_lines.2<-nrow(fread(file.path(dirs[dd],'coding_only.tsv'), header = T))
    print(paste0("cat ", dirs[d], '/coding_only.tsv', ' > ', paste0(opt$o,'/', basename(dirs[d]), '_', basename(dirs[dd]), '.tsv')))
    print(paste0("tail -n ", nr_lines.2,' ', dirs[dd], '/coding_only.tsv', ' >> ', paste0(opt$o, '/', basename(dirs[d]), '_', basename(dirs[dd]), '.tsv')))
    print(paste0(basename(dirs[d]), '_', basename(dirs[dd]), 'clusters.tsv'))
    da<-data.table(seq(1:(nr_lines+nr_lines.2)), c(rep(1, nr_lines), rep(2, nr_lines.2)))
   # fwrite(da, file = paste0(basename(dirs[d]), '_', basename(dirs[dd]), 'clusters.tsv'), sep='\t')
  }
}


 