#!/usr/bin/env Rscript

packages <- c("optparse","stringr","rtracklayer","curl")

zzz<-lapply(packages, function(xxx) suppressMessages(require(xxx, character.only = TRUE,quietly=TRUE,warn.conflicts = FALSE)))

option_list = list(
  make_option(c("-i", "--input_file"), type="character", default=NULL,
              help="input file", metavar="character"),
  make_option(c("-g", "--gene_annotation"), type="character", default=NULL,
          help="input gene annotation", metavar="character"),
  make_option(c("-o", "--out_file"), type="character", default=NULL,
              help="output file", metavar="character")
);

example.use <- "Example: Rscript $HOME/Aimin/ispipe/R/makeREFIndex.R -i ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh38.primary_assembly.genome.fa.gz -g ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.annotation.gtf.gz -o /home/ayan/Aimin/ispipe/utilsRefData/hg38/GRCh38.primary_assembly.genome.fa\n"

opt_parser = OptionParser(option_list=option_list,epilogue=example.use);
opt = parse_args(opt_parser);

if (is.null(opt$input_file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input count file)", call.=FALSE)
}

input.file <- opt$input_file
input.gene.annotation <- opt$gene_annotation
output.file <- opt$out_file

out.dir.name <- dirname(output.file)
if (!dir.exists(out.dir.name)){dir.create(out.dir.name, recursive = TRUE)}

t.input <- file.info(input.file)$mtime
t.output <- file.info(output.file)$mtime

x <- dirname(output.file)
y <- basename(x)
z <- file.path(x,paste0(y,"ChrOnly.fa"))

if(is.na(t.output)){
  
  cat("Download fasta file\n")
  cmd0 <- paste0('wget ',input.file,' -P ',x)
  system(cmd0)
  
  cat("gunzip file\n")
  cmd1 <- paste0('gunzip ',paste0(output.file,'.gz'))
  system(cmd1)
  
  cat("samtools faidx fasta file\n")
  cmd2 <- paste0('samtools faidx ',output.file)
  system(cmd2)
  
  
  cat("get chr1-22 and X, Y,M fasta file only\n")
  chr<-c(paste("chr",1:22,sep=""),"chrX","chrY","chrM")
  
  # chr
  # output.file <- '/home/ayan/Aimin/ispipe/utilsRefData/hg38/GRCh38.primary_assembly.genome.fa'
  # x <- dirname(output.file)
  # y <- basename(x)
  # z <- file.path(x,paste0(y,"ChrOnly.fa"))
  
  z.output <- file.info(z)$mtime
  
  if(is.na(z.output)){
     
  null <- lapply(chr, function(u){
    
    cmd3 <- paste0('samtools faidx ',output.file,' ',u,' >> ',z)
    system(cmd3)
    
  })
    
  }
  
  cat("make bwa index\n")
  cmd4= paste0('bwa index -a bwtsw ',z)
  system(cmd4)
  
}

w <- file.path(x,paste0('repeatMasker',y,'BED'))

t.w <- file.info(w)$mtime

if(is.na(t.w)){
  
cat("get repeatMaskerBED\n")

mySession = browserSession("UCSC")
genome(mySession) <- y

query <- ucscTableQuery(mySession, "rmsk")

rmsk.table <- track(query)
con <- file(w)
export(rmsk.table,con,format = "bed",trackLine = FALSE)
}

gene.file <- file.path(x,paste0(y,'_genesKNOWN_sorted.bed'))
gene.file.time <- file.info(gene.file)$mtime

if(is.na(gene.file.time)){
  
  cat("get genesKNOWN_sorted.bed\n")
  
  url = input.gene.annotation
  tmp <- tempfile()
  curl_download(url, tmp)
  
  zz=gzfile(tmp,'rt')
  
  dat=read.table(zz,header=F,sep = "\t")
  
  dat.gene <- dat[dat$V3=='gene',]
  
  gene.name <- str_sub(unlist(lapply(base::strsplit(as.character(dat.gene$V9),';'),function(u)u[3])),start = 12L)
  
  gene.hg38.table <- data.frame(chr=dat.gene$V1,start=dat.gene$V4,end=dat.gene$V5,strand=dat.gene$V7,geneName=gene.name,type=rep('KNOWN',length(gene.name)))
  
  write.table(gene.hg38.table,file = gene.file, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = FALSE,col.names = FALSE)

}

outputFile1 <- file.path(x,"hg38.genome.sorted.txt")
outputFile2 <- file.path(x,"hg38.genome.num.txt")

hg38.genome.num.time <- file.info(outputFile2)$mtime

if(is.na(hg38.genome.num.time)){
  
cat("get hg38.genome.num file\n")
  
hg38.genome <- as.data.frame(SeqinfoForUCSCGenome(y))
hg38.genome.size <- data.frame(chr=row.names(hg38.genome),length=hg38.genome$seqlengths)
hg38.genome.num <- data.frame(chr.index=seq(1,length(row.names(hg38.genome))),chr=row.names(hg38.genome))

write.table(hg38.genome.size,file = outputFile1, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = FALSE,col.names = FALSE)

write.table(hg38.genome.num,file = outputFile2, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = FALSE,col.names = FALSE)

}
