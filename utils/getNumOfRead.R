#!/usr/bin/env Rscript

packages <- c("optparse","stringr","ggplot2")
zzz<-lapply(packages, function(xxx) suppressMessages(require(xxx, character.only = TRUE,quietly=TRUE,warn.conflicts = FALSE)))

option_list = list(
  make_option(c("-i", "--input_file"), type="character", default=NULL,
              help="input file", metavar="character"),
  make_option(c("-p", "--input_pattern"), type="character", default=NULL,
              help="input pattern", metavar="character"),
  make_option(c("-o", "--out_file"), type="character", default=NULL,
              help="output file", metavar="character")
);

example.use <- "Example: Rscript $HOME/Aimin/ispipe/utils/getNumOfRead.R -i /home/ayan/Aimin/Seagate/ISseqOutput/Sept2/CutAdapt/align -p sam -o /home/ayan/Aimin/Seagate/ISseqOutput/Sept2/TotalReads/totalReads.rds\n"

opt_parser = OptionParser(option_list=option_list,epilogue=example.use);
opt = parse_args(opt_parser);

if (is.null(opt$input_file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input count file)", call.=FALSE)
}

input.file <- opt$input_file
output.file <- opt$out_file
input.pattern <- opt$input_pattern

#input.file <- '/home/ayan/Aimin/Seagate/ISseqOutput/Sept2/CutAdapt/align'
#output.file <- '/home/ayan/Aimin/Seagate/ISseqOutput/Sept2/TotalReads/totalReads.rds'

getReadCount <- function(input.file,input.pattern,output.file) {

  #input.file <- '/home/ayan/Aimin/Seagate/ISseqOutput/Sept2/CutAdapt/missingIS' 
  #input.pattern <-'*mem_sort_inMask_qual.bam'
  
  cat(input.file,'\n')
  cat(input.pattern,'\n')

    
  fastq.files <- list.files(path = input.file,pattern=paste0(input.pattern,"$"),full.names = TRUE, recursive = F)
  read.count.stat <- lapply(1:length(fastq.files), function(i){
    
  u <- fastq.files[[i]]
  
  cmd = paste0('samtools flagstat ',u)
  
  cmd.out <- system(cmd, intern = TRUE)
  cmd.out.1 <- strsplit(cmd.out,split=' ')
  cmd.out.2 <- lapply(cmd.out.1,function(u){
    z <- data.frame(QCpassedReads=u[1],QCfailedReads=u[3])
    z
  })
  cmd.out.3 <- do.call(rbind,cmd.out.2)
  row.names(cmd.out.3) <- c('total','secondary','supplementary','duplicates','mapped','PairedInSequencing','read1','read2','properlyPaired','with_itself_and_mate_mapped','singletons','with_mate_mapped_to_a_different_chr','with_mate_mapped_to_a_different_chr_mapQ_above_or_equal_5')
  cmd.out.3
  
})
  
  names(read.count.stat) <- basename(fastq.files)
  
  x <- names(read.count.stat)
  
  y <- str_sub(x,unlist(lapply(str_locate_all(x,'_'),"[",c(3)))+1,unlist(lapply(str_locate_all(x,"_"),"[",c(5)))-1)
  
  names(read.count.stat) <- y
  
  T <- lapply(read.count.stat, function(u){
    
    t <- as.integer(as.character(u[1,1]))
    
  })
  
  total.read.data <- do.call(rbind,T)
  
  # output.dir <- dirname(output.file)
  # if(!dir.exists(output.dir)){dir.create(output.dir,recursive = TRUE)}
  # saveRDS(total.read.data,file = output.file)
  
  total.read.data
  
}

getMapped <- function(input.file,input.pattern,output.file) {
  
  #input.file <- '/home/ayan/Aimin/Seagate/ISseqOutput/Sept2/CutAdapt/missingIS' 
  #input.pattern <-'*mem_sort_inMask_qual.bam'
  
  cat(input.file,'\n')
  cat(input.pattern,'\n')
  
  
  fastq.files <- list.files(path = input.file,pattern=paste0(input.pattern,"$"),full.names = TRUE, recursive = F)
  read.count.stat <- lapply(1:length(fastq.files), function(i){
    
    u <- fastq.files[[i]]
    
    cmd = paste0('samtools flagstat ',u)
    
    cmd.out <- system(cmd, intern = TRUE)
    cmd.out.1 <- strsplit(cmd.out,split=' ')
    cmd.out.2 <- lapply(cmd.out.1,function(u){
      z <- data.frame(QCpassedReads=u[1],QCfailedReads=u[3])
      z
    })
    cmd.out.3 <- do.call(rbind,cmd.out.2)
    row.names(cmd.out.3) <- c('total','secondary','supplementary','duplicates','mapped','PairedInSequencing','read1','read2','properlyPaired','with_itself_and_mate_mapped','singletons','with_mate_mapped_to_a_different_chr','with_mate_mapped_to_a_different_chr_mapQ_above_or_equal_5')
    cmd.out.3
    
  })
  
  names(read.count.stat) <- basename(fastq.files)
  
  x <- names(read.count.stat)
  
  y <- str_sub(x,unlist(lapply(str_locate_all(x,'_'),"[",c(3)))+1,unlist(lapply(str_locate_all(x,"_"),"[",c(5)))-1)
  
  names(read.count.stat) <- y
  
  T <- lapply(read.count.stat, function(u){
    
    t <- as.integer(as.character(u[5,1]))
    
  })
  
  total.read.data <- do.call(rbind,T)
  
  # output.dir <- dirname(output.file)
  # if(!dir.exists(output.dir)){dir.create(output.dir,recursive = TRUE)}
  # saveRDS(total.read.data,file = output.file)
  
  total.read.data
  
}

T <- getReadCount(input.file,input.pattern,output.file)
M <- getMapped(input.file,input.pattern,output.file)

T.M <- merge(T,M,by=0)

m <- T.M[,3]/T.M[,2]

T.M.1 <- data.frame(T.M,Percentage=paste(round(100*m, 2), "%", sep=""))

colnames(T.M.1) <- c('barcode','Total','Mapped','Mapped/Total')

T.M.1

output.dir <- dirname(output.file)
if(!dir.exists(output.dir)){dir.create(output.dir,recursive = TRUE)}
saveRDS(T.M.1,file = output.file)


# association.file <- "/media/aws_share/D32_Platform_Development/test/ISAtest/MiSeqTest/Association_pool_ISA_AVRO_TEST4_add_hg19.csv"
# 
# y.lab <- "Number of total reads"
# angle <- 90
# size <- 10
# plot.title <- "Number of total reads"
# 
# getTotalReadsPlot <- function(total.read.data,association.file,y.lab,angle,size,plot.title) {
# 
#   data.association <- read.csv(association.file)
#   
#   index <-paste0(data.association$Fusion.Primer.LTR..ID,'_',data.association$Fusion.Primer.LC..ID)
#   sample.name <- unique(paste0(data.association$Sample.Type,'_',data.association$TimePoint.R.d..C.m.))
#   
#   z <- data.frame(sampleName=sample.name)
#   row.names(z) <- index
#   
#   res <- merge(total.read.data,z,by=0)
#   
#   df <- data.frame(Time=res$sampleName,IScount=res$V1)
#   
#   Category <- as.character(df$Time)
#   color.codes <- rainbow(length(Category))
#   
#   df$Time <- factor(df$Time,levels = Category)
#   
#   sp <- ggplot(df, aes(x = Time, y = IScount, fill = Time)) +
#     geom_bar(stat="identity",position = "dodge") + geom_text(aes(label=IScount)) + theme(plot.title=element_text(hjust=0.5)) +
#     ylab(y.lab)+xlab("trt")+ggtitle(plot.title) + theme_bw() +theme(axis.text.x = element_text(angle = angle, hjust = 1,size = size),legend.position = "bottom", legend.justification = 0)
#   
# }
