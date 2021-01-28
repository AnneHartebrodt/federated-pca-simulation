suppressMessages(require(data.table))
suppressMessages(require(ggplot2))
suppressMessages(require(tidyr))
suppressMessages(require(cowplot))
suppressMessages(require(optparse))
suppressMessages(require(dplyr))
suppressMessages(require(stringr))
suppressMessages(require(optparse))

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