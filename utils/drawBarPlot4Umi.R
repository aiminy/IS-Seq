#!/usr/bin/env Rscript

initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)
other.name <- file.path(script.basename, "ipak.R")
print(paste("Sourcing",other.name,"from",script.name))
source(other.name)

# CRAN PACKAGES
cranpackages <- c("curl","stringr","ggplot2","ggpubr","cowplot","usethis", "covr", "httr", "roxygen2", "rversions","devtools","plyr","optparse")
ipak(cranpackages, repository='CRAN')

# # BIOCONDUCTOR
# #biocpackages <- c("AnnotationDbi", "baySeq", "Biobase", "BiocGenerics", 
#                   "BiocParallel", "DEDS", "DESeq2", "EBSeq", "edgeR", "IHW", "iCOBRA",
#                   "limma", "Linnorm", "MAST", "monocle", "NOISeq", "qvalue", "ROTS", "RUVSeq", 
#                   "S4Vectors", "scater", "scDD", "scde", "scone", "scran", "SCnorm", 
#                   "SingleCellExperiment", "SummarizedExperiment", "zinbwave")
# # ipak(biocpackages, repository='Bioconductor')
# 
# # GITHUB
# #githubpackages <- c('nghiavtr/BPSC', 'VCCRI/cidr', 'cz-ye/DECENT', 
#                     'mohuangx/SAVER', 'statOmics/zingeR')
# #ipak(githubpackages, repository = 'github')

packages <- c("stringr","ggplot2","ggpubr","plyr","optparse")
zzz<-lapply(packages, function(xxx) suppressMessages(require(xxx, character.only = TRUE,quietly=TRUE,warn.conflicts = FALSE)))

option_list = list(
  make_option(c("-i", "--input_file_dir"), type="character", default=NULL,
              help="input file dir", metavar="character"),
  make_option(c("-o", "--output_file_dir"), type="character", default="./",
              help="output file dir name [default= %default]", metavar="character")
);

example.use <- "Example:
      Rscript /home/ayan/Aimin/ispipe/utils/drawBarPlot4Umi.R -i /home/ayan/Aimin/ispipe/ISseqOutput/Mar06/UmiBased/UmiCluster -o /home/ayan/Aimin/ispipe/ISseqOutput/Mar06/UmiBased/UmiCluster/UmiBarPlots"

opt_parser = OptionParser(option_list=option_list,epilogue=example.use);
opt = parse_args(opt_parser);

if (is.null(opt$input_file_dir)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file dir)", call.=FALSE)
}

drawBarPlot4Umi <- function(input.file,output.dir) {

  if(!dir.exists(output.dir)){dir.create(output.dir,recursive = TRUE)}     
  
  input.file.pattern= "*umi_cluster$"
  
  file.name.4 <- list.files(input.file,pattern=input.file.pattern,all.files = TRUE,full.names = TRUE,recursive = FALSE,include.dirs = TRUE)     
  
  inds <- file.size(file.name.4) != 0 
  X <- file.name.4[inds]
  
  re.out<-lapply(1:length(X),function(u,X){
    
   input <- read.table(X[u],header = F)
   fileName <-  basename(X[u])
   y <- str_sub(fileName,str_locate(fileName,"Barcode")[,2]+2,str_locate(fileName,".fq")[,1]-1)
   z <- data.frame(input[,-3],fName=rep(y,dim(input)[1]))
   z
  
  },X)
  
  z.FileName <-  do.call(rbind.data.frame,re.out)
  
  fName <- unique(as.character(z.FileName$fName))
  
  null <- lapply(fName, function(u,z.FileName){
    
  df <- z.FileName[z.FileName$fName==u,]
  df <- df[1:100,]
  Category <- as.character(df$V1)
  color.codes <- rainbow(length(Category))
  df$V1 <- factor(df$V1, levels = df$V1[order(df$V2,decreasing = TRUE)])
  
  sp <- ggplot(df, aes(x = V1, y = V2, fill = color.codes)) +
    geom_bar(stat="identity",position = "dodge")  +
    theme(plot.title=element_text(hjust=0.5)) +  
    ylab("Count")  +  theme(legend.position="none")+
    theme(axis.text.x = element_text(angle = 90, hjust = 1,size = 4))+xlab("UMI")
  ggtitle(u)
  
  spp <- list(sp=sp)
  multi.page <- ggarrange(plotlist=spp,nrow = 1, ncol = 1)
  
  if(!dir.exists(output.dir)){dir.create(output.dir,recursive = TRUE)}
  ggexport(multi.page, filename = file.path(output.dir,paste0(u,"_UMI.pdf")))
  
  },z.FileName)

}

input.file <- opt$input_file_dir
output.dir <- opt$output_file_dir

null <- drawBarPlot4Umi(input.file,output.dir)

save.image(file=file.path(output.dir,"umi.RData"))
