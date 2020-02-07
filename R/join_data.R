require(data.table)
require(optparse)

option_list = list(
  make_option(c("-d", "--dir"), action="store", default=NA, type='character',
              help="coding genes output file")
)
opt = parse_args(OptionParser(option_list=option_list))

dirs<- list.dirs(opt$d, recursive = F)

for(d in 1: length(dirs)){
  nr_lines<-nrow(fread(file.path(dirs[d],'coding_only.tsv')))
  for(dd in d:length(dirs)){
    nr_lines.2<-nrow(fread(file.path(dirs[dd],'coding_only.tsv')))
    print(paste0("cat ", dirs[d], '/coding_only.tsv', ' > ', basename(dirs[d]), '_', basename(dirs[dd]), '.tsv'))
    print(paste0("cat ", dirs[dd], '/coding_only.tsv', ' >> ', basename(dirs[d]), '_', basename(dirs[dd]), '.tsv'))
    dd<-data.table(seq(1:(nr_lines+nr_lines.2)), c(rep(1, nr_lines), rep(2, nr_lines.2)))
    fwrite(dd, file = paste0(basename(dirs[d]), '_', basename(dirs[dd]), 'clusters.tsv'), sep='\t')
  }
}


 