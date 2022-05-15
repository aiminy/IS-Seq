#!/usr/bin/env Rscript

#packages <- c("optparse","stringr","hiReadsProcessor")

packages <- c("optparse","stringr","GenomicRanges")

zzz<-lapply(packages, function(xxx) suppressMessages(require(xxx, character.only = TRUE,quietly=TRUE,warn.conflicts = FALSE)))

option_list = list(
  make_option(c("-i", "--input_file"), type="character", default=NULL,
              help="input file", metavar="character"),
  make_option(c("-r", "--ref_fa_file"), type="character", default=NULL,
              help="input ref fasta file", metavar="character"),
  make_option(c("-q", "--identity_threshold"), type="double", default=NULL,
              help="identity threshold value", metavar="double"),
  make_option(c("-o", "--out_file"), type="character", default=NULL,
              help="output file", metavar="character")
);

example.use <- "Example: Rscript $HOME/Aimin/ispipe/R/getMissIS.R -i /home/ayan/Aimin/Seagate/ISseqOutput/June24/CutAdapt/BAMSorted/R1_R2_Barcode_FB-P5-Rd1-LTR.1_FB-P7-Rd2-LC.1_aligned_mem_sort_inMask.bam -r ~/Aimin/ispipe/utilsRefData/mm10.ebwt/mm10ChrOnly.fa -q 95 -o /home/ayan/Aimin/Seagate/ISseqOutput/June24/CutAdapt/missingIS/R1_R2_Barcode_FB-P5-Rd1-LTR.1_FB-P7-Rd2-LC.1_aligned_mem_sort_inMask_missingIS.txt\n"

opt_parser = OptionParser(option_list=option_list,epilogue=example.use);
opt = parse_args(opt_parser);

if (is.null(opt$input_file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input count file)", call.=FALSE)
}

input.file <- opt$input_file
output.file <- opt$out_file

ref.fa <- opt$ref_fa_file
identity.threshold <- opt$identity_threshold

out.dir.name <- dirname(output.file)
if (!dir.exists(out.dir.name)){dir.create(out.dir.name, recursive = TRUE)}

x_name <- tools::file_path_sans_ext(basename(input.file))

input_file_le_30_bam <- file.path(out.dir.name,paste0(x_name,"_le_30.bam"))
input_file_qual_bam <- file.path(out.dir.name,paste0(x_name,"_qual.bam"))
input_file_le_30_bam_fa <- paste0(tools::file_path_sans_ext(input_file_le_30_bam),'.fa')
output_all.psl <- paste0(tools::file_path_sans_ext(input_file_le_30_bam),'.psl')
output_all_start_0_1_psl <- paste0(tools::file_path_sans_ext(input_file_le_30_bam),'_start_0_1.psl')

t.input <- file.info(input.file)$mtime
t.output <- file.info(output_all_start_0_1_psl)$mtime

if(is.na(t.output)|(t.input>t.output)){
  
cmd0 <- paste0('samtools view -bq 30 ',input.file,' -U ',input_file_le_30_bam,' -o ',input_file_qual_bam)
system(cmd0)

cmd1 <- paste0('samtools bam2fq ',input_file_le_30_bam,' | seqtk seq -A - > ',input_file_le_30_bam_fa)
system(cmd1)

cmd2 <- paste0('blat -stepSize=5 -repMatch=2253 -minScore=20 -minIdentity=0 ', ref.fa,' ',input_file_le_30_bam_fa,' ',output_all.psl)
cat(cmd2,"\n")
system(cmd2)

cmd3 <- paste0("awk '$12==0' ",output_all.psl,' | grep /1 > ',output_all_start_0_1_psl)
cat(cmd3,"\n")
system(cmd3)

}

#cols <- pslCols()

cols <- c("matches","misMatches","repMatches","nCount","qNumInsert","qBaseInsert","tNumInsert","tBaseInsert","strand" ,"qName","qSize","qStart","qEnd","tName","tSize","tStart","tEnd","blockCount","blockSizes","qStarts","tStarts")


#output_all_start_0_1_psl <- '~/Aimin/NoMachineLinux/home/ayan/Aimin/Seagate/ISseqOutput/June24/CutAdapt/missingIS/R1_R2_Barcode_FB-P5-Rd1-LTR.1_FB-P7-Rd2-LC.1_aligned_mem_sort_inMask_le_30_start_0_1.psl'

#output_all_start_0_1_psl <- '~/Aimin/NoMachineLinux/home/ayan/Aimin/Seagate/ISseqOutput/June24/CutAdapt/missingIS/R1_R2_Barcode_FB-P5-Rd1-LTR.1_FB-P7-Rd2-LC.4_aligned_mem_sort_inMask_le_30_start_0_1.psl'

#output_all_start_0_1_psl <- '~/Aimin/NoMachineLinux/home/ayan/Aimin/Seagate/ISseqOutput/June24/CutAdapt/missingIS/R1_R2_Barcode_FB-P5-Rd1-LTR.1_FB-P7-Rd2-LC.3_aligned_mem_sort_inMask_le_30_start_0_1.psl'
 

s <- file.info(output_all_start_0_1_psl)$size

if(s!=0) {
	
psl.df <- read.table(output_all_start_0_1_psl,sep="\t",header = FALSE,stringsAsFactors = FALSE)

colnames(psl.df) <- cols

calculateIdentityScore <- function(psl.df.IS.1) {
  y=rep(1,dim(psl.df.IS.1)[1])
  
  mp=psl.df.IS.1$misMatches
  gq=psl.df.IS.1$qNumInsert
  
  bq=psl.df.IS.1$qStart
  eq=psl.df.IS.1$qEnd
  
  bt=psl.df.IS.1$tStart
  et=psl.df.IS.1$tEnd
  
  zq <- eq-bq
  
  zt <- et-bt
  
  z = zq-zt
  
  z[which(z<0)] <- 0
  
  m=psl.df.IS.1$matches
  
  r=psl.df.IS.1$repMatches
  
  t = y*(m+r+mp)
  
  
  d = 1000*(y*mp+gq+round(3*log(1+z)))/t
  
  p=100-0.1*d
}

psl.df$Identity <- calculateIdentityScore(psl.df)
identity.threshold <- 95
psl.df <- psl.df[which(psl.df$Identity>=identity.threshold),]

psl.df$chr_pos <- paste0(psl.df$tName,"_",psl.df$tStart)

pos.count <- as.data.frame(table(psl.df$chr_pos))

colnames(pos.count) <- c('chr_pos','count')

pos.count.sorted <- pos.count[order(pos.count$count,decreasing = T),]

head(pos.count.sorted,20)

psl.df.count <- merge(psl.df,pos.count,by='chr_pos')

psl.df.count.1 <- psl.df.count[,c(15,17,10,24)]

psl.df.count.2 <- unique(psl.df.count.1)

chrOrder<-c(paste("chr",1:22,sep=""),"chrX","chrY")
  
psl.df.count.2$tName <- factor(psl.df.count.2$tName, levels=chrOrder)
  
psl.df.count.3 <- psl.df.count.2[order(psl.df.count.2$tName,psl.df.count.2$tStart),]

psl.df.count.3$tName_strand <- paste0(psl.df.count.3$tName,"_",psl.df.count.3$strand)

chr.strand <- unique(psl.df.count.3$tName_strand)

temp <- psl.df.count.3[,c(1,2,2,3,4)]
colnames(temp) <- c('chr','start','end','strand','count')

gr <- GenomicRanges::makeGRangesFromDataFrame(temp, keep.extra.columns=TRUE)

gr1 <- GenomicRanges::reduce(gr,min.gapwidth=7L,with.revmap=T)


if(length(which(width(gr1)!=1))==0){
  
  res  <- as.data.frame(gr)
  
  res$seqnames <- factor(res$seqnames, levels=chrOrder)
  
  res <- res[order(res$seqnames,res$start),]
  
  results  <- data.frame(chr=res$seqnames,start=res$start,chr.index=as.integer(res$seqnames),strand=res$strand,count=res$count)
}else{
gr2 <- gr1[which(width(gr1)!=1)]

gr3 <- gr1[-which(width(gr1)!=1)]

index <- gr2$revmap
Y <- lapply(index, function(u){
  
  y <- gr[unlist(u)]
  
})

mergeIS <- function(W,dist.cutoff) {
  
  if(dist.cutoff==0){
    res <- W
  }else{
    n <- dim(W)[1]
    
  if(n==2){
    
    d <- W[2,]$start-W[1,]$start
    
    if(d>dist.cutoff)
    {
      res <- W
    }else{
      res <- W[which.max(W$count),]
      res$count <- sum(W$count)
    }
    
  }else
  {
    pivot <- W[1,]$start
    
    d <- W$start-pivot
    
    index <- which(d<=dist.cutoff)
    
    if(length(index)!=1)
    {
      W.pivot <- W[index,]
      temp <- W.pivot[which.max(W.pivot$count),]
      temp$count <- sum(W.pivot$count)
    
      if(dim(W[-index,])[1]!=0){
        res <- rbind(temp,W[-index,])
        res <- mergeIS(res,dist.cutoff)
      }else{
        res <- temp
        }
      
    }else
    {
      temp <- W[-1,]
      res <- rbind(W[1,],mergeIS(temp,dist.cutoff))
    }
    
  }
    
  }
  
  res
  
}

YY <- lapply(1:length(Y), function(u){
  
  w.1 <- mergeIS(as.data.frame(Y[[u]]),dist.cutoff = 7)
  
})

YYY <- do.call(rbind,YY)

res <- rbind(as.data.frame(gr[unlist(gr3$revmap)]),YYY)

res$seqnames <- factor(res$seqnames, levels=chrOrder)

res <- res[order(res$seqnames,res$start),]

results  <- data.frame(chr=res$seqnames,start=res$start,chr.index=as.integer(res$seqnames),strand=res$strand,count=res$count) 
}

#x <- paste0("POOL-ISA-AVRO-TEST2_",substr(basename(output_all_start_0_1_psl),15,nchar(basename(output_all_start_0_1_psl))-44),"_final_parse_filterNo_grouped_IS.txt")

#output.file <- file.path('~/Aimin/NoMachineLinux/home/ayan/Aimin/Seagate/ISseqOutput/June24/CutAdapt/missingIS',x)
 
write.table(results,file=output.file ,sep = "\t",row.names = F,quote = F,col.names = F)
}

#t1 <- file.info('~/Aimin/NoMachineLinux/home/ayan/Aimin/Seagate/ISseqOutput/June24/CutAdapt/BAMSorted/R1_R2_Barcode_FB-P5-Rd1-LTR.1_FB-P7-Rd2-LC.1_aligned_mem_sort_inMask.bam')$mtime


#t2 <- file.info('~/Aimin/NoMachineLinux/home/ayan/Aimin/Seagate/ISseqOutput/June24/CutAdapt/missingIS/R1_R2_Barcode_FB-P5-Rd1-LTR.1_FB-P7-Rd2-LC.3_aligned_mem_sort_inMask_le_30_start_0_1.psl')$mtime



