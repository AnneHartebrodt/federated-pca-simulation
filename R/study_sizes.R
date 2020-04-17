require(data.table)
require(optparse)
require(jsonlite)

option_list = list(
  make_option(c("-f", "--inputfile"), action="store", default=NA, type='character',
              help="annoation file in gtf format"),
  make_option(c("-o", "--outputfile"), action="store", default=NA, type='character',
              help="The directory for plot")
)
opt = parse_args(OptionParser(option_list=option_list))
# The cancer genome atlas
#Data has been dowloaded from the cancer genome atlas website (https://portal.gdc.cancer.gov/) using
#this query 
#'files.access in ["open"] and files.analysis.workflow_type in ["HTSeq - FPKM"] and files.data_type in 
#'["Gene Expression Quantification"] and files.experimental_strategy in ["RNA-Seq"]' using the json download.

print (opt$inputfile)
f<-as.data.table(fromJSON(opt$inputfile, simplifyDataFrame = T, flatten = T))
f<-as.data.table(tidyr::unnest(f))
tcga<-f[,.N, by=project.project_id]
tcga$database<-'TCGA'
fwrite(tcga, file = opt$outputfile, sep = '\t')
