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
dd <-list()
for(d in dirs){
if(file.exists(file.path(d, 'coding_only.tsv'))){
  dd[d]<-d
}
}
dirs <- as.character(unlist(dd))
for (k in c(2,5)){

  for(i in 1:3){
    fil <-sample(dirs, k, replace = FALSE)
    dal <- list()
    cl = 1
    for(f in fil){
      nr_lines<-nrow(fread(input = file.path(f,'coding_only.tsv')))
      if(cl ==1){
      system(paste0("cat ", f, '/coding_only.tsv', ' > ', paste0(opt$o,'/', k, '_', i, '.tsv')))
      }else{
        system(paste0("tail -n ", nr_lines,' ', f, '/coding_only.tsv', ' >> ', paste0(opt$o,'/', k, '_', i, '.tsv')))
      }
      dal[[cl]]<-data.table(seq(1:nr_lines), rep(cl, nr_lines), basename(f))
      cl = cl+1
    }
    dal <-rbindlist(dal)
    fwrite(dal, file = paste0(opt$o, '/clusters/', k, '_', i, '_clusters.tsv'), sep='\t')
  }
  
}
