#!/usr/bin/env Rscript

packages <- c("optparse","stringr","rtracklayer","curl")

zzz<-lapply(packages, function(xxx) suppressMessages(require(xxx, character.only = TRUE,quietly=TRUE,warn.conflicts = FALSE)))

option_list = list(
  make_option(c("-i", "--input_genome"), type="character", default=NULL,
              help="input genome", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="output", metavar="character")
);

example.use <- "Example: Rscript $HOME/ispipe/R/makeREFIndex4INSPIIRED.R -i hg18 -o /home/ubuntu/SHARE/ISseqOutput/INSPIIREDData/utilsRefData\n"

opt_parser = OptionParser(option_list=option_list,epilogue=example.use);
opt = parse_args(opt_parser);

if (is.null(opt$input_genome)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input count file)", call.=FALSE)
}

input.genome <- opt$input_genome

out.dir.name <- file.path(opt$output,input.genome)

if (!dir.exists(out.dir.name)){dir.create(out.dir.name, recursive = TRUE)}

z <- file.path(out.dir.name,paste0(input.genome,"ChrOnly.fa"))

print(input.genome)

if(input.genome=='hg18'){

library(BSgenome.Hsapiens.UCSC.hg18)
mySession = browserSession("UCSC")
genome(mySession) <- input.genome

library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg18.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg18.knownGene

}

hg18.fa <- getSeq(BSgenome.Hsapiens.UCSC.hg18)

cat("get chr1-22 and X, Y,M fasta file only\n")

chr<-c(paste("chr",1:22,sep=""),"chrX","chrY","chrM")

hg18.fa.chr.only <- hg18.fa[which(names(hg18.fa) %in% chr)]


  z.output <- file.info(z)$mtime
  
  if(is.na(z.output)){

   library(ShortRead)
   writeFasta(hg18.fa.chr.only,file=z)
  
   cat("make bwa index\n")
   cmd4= paste0('bwa index -a bwtsw ',z)
   system(cmd4)


   fa_2bit_out <- file.path(out.dir.name,paste0(input.genome,"ChrOnly.2bit"))

   print(fa_2bit_out)

   cmd=paste0('faToTwoBit ',z,' ',fa_2bit_out)
   system(cmd)
  
  }

  
w <- file.path(out.dir.name,paste0('repeatMasker',input.genome,'BED'))

t.w <- file.info(w)$mtime

if(is.na(t.w)){
  
cat("get repeatMaskerBED\n")

query <- ucscTableQuery(genome(mySession), table = "rmsk")

table <- getTable(query)

con <- file(w)

write.table(table[,c(6:8,11,2,10)],file=con,sep='\t', row.names=FALSE,col.names=FALSE,quote=FALSE)

}

gene.file <- file.path(out.dir.name,paste0(input.genome,'_genesKNOWN_sorted.bed'))
gene.file.time <- file.info(gene.file)$mtime


GetAnnotation <- function(txdb,chr,output.file) {

  gn <- sort(genes(txdb))

  myGeneSymbols <- select(org.Hs.eg.db,keys = mcols(gn)$gene_id,columns = c("SYMBOL","ENTREZID"),keytype = "ENTREZID")
  
  gn.1 <- data.frame(as.data.frame(gn),SYMBOL=myGeneSymbols$SYMBOL)
  
  gn.1[which(is.na(gn.1$SYMBOL)),]$SYMBOL <- as.character(gn.1[which(is.na(gn.1$SYMBOL)),]$gene_id)
  
  chr<-c(paste("chr",1:22,sep=""),"chrX","chrY","chrM")
  
  gn.2 <- gn.1[which(gn.1$seqnames %in% chr),]
  
  gn.3 <- data.frame(chr=gn.2$seqnames,start=gn.2$start,end=gn.2$end,strand=gn.2$strand,geneName=gn.2$SYMBOL,KNOWN=rep('KNOWN',dim(gn.2)[1]))
  
  write.table(gn.3,file=output.file, sep='\t', row.names=FALSE,col.names=FALSE,quote=FALSE)
}

if(is.na(gene.file.time)){
  
  cat("get genesKNOWN_sorted.bed\n")
  
  GetAnnotation(txdb,chr,gene.file)

}

outputFile1 <- file.path(out.dir.name,paste0(input.genome,".genome.sorted.txt"))
outputFile2 <- file.path(out.dir.name,paste0(input.genome,".genome.num.txt"))

hg18.genome.num.time <- file.info(outputFile2)$mtime

if(is.na(hg18.genome.num.time)){
  
cat("get hg18.genome.num file\n")
  
hg18.genome <- as.data.frame(SeqinfoForUCSCGenome(input.genome))
hg18.genome.size <- data.frame(chr=row.names(hg18.genome),length=hg18.genome$seqlengths)
hg18.genome.num <- data.frame(chr.index=seq(1,length(row.names(hg18.genome))),chr=row.names(hg18.genome))

write.table(hg18.genome.size,file = outputFile1, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = FALSE,col.names = FALSE)

write.table(hg18.genome.num,file = outputFile2, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = FALSE,col.names = FALSE)

}

