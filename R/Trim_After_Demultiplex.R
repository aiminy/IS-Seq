# Rscript /home/ubuntu/intsitecaller/testCases/intSiteValidation/completeMetadata.RData /home/ubuntu/intsitecaller/MyTest/demultiplexedReps/clone1-1_R1.fastq.gz /home/ubuntu/intsitecaller/MyTest/demultiplexedReps/clone1-1_R2.fastq.gz /home/ubuntu/intsitecaller/testCases/intSiteValidation/p746vector.fasta /home/ubuntu/intsitecaller/testCases/intSiteValidation/hg18.2bit ~/SHARE/ISseqOutput/INSPIIRED_test_run

args = commandArgs(trailingOnly=TRUE)
if (length(args)<6) {
  stop("6 arguments are needed", call.=FALSE)
} else{
  completeMetadata=args[1]
  r1=args[2] 
  r2=args[3]
  vectorSeq=args[4]
  ref.bit=args[5]
  output.dir==args[6]
}

if(!dir.exists(output.dir)){dir.create(output.dir,recursive = TRUE)}

library(ShortRead)

#completeMetadata <- get(load('/home/ubuntu/intsitecaller/testCases/intSiteValidation/completeMetadata.RData'))

completeMetadata <- get(load(completeMetadata))

alias <- completeMetadata[1,]$alias

stats.bore <- data.frame(sample=alias)

#read1 <- '/home/ubuntu/intsitecaller/MyTest/demultiplexedReps/clone1-1_R1.fastq.gz'
#read2 <- '/home/ubuntu/intsitecaller/MyTest/demultiplexedReps/clone1-1_R2.fastq.gz'

read1 <- r1
read2 <- r2

vectorSeq <- vectorSeq
ref.bit <-  ref.bit

#'/home/ubuntu/intsitecaller/testCases/intSiteValidation/hg18.2bit'




reads <- lapply(list(read1, read2), sapply, readFastq)

stats.bore$barcoded <- sum(sapply(reads[[1]], length))

r <- lapply(reads, function(x){
  seqs <- x[[1]]
  if(length(seqs) > 0){
    #remove anything after 5 bases under Q30 in 10bp window
    ##trimmed <- trimTailw(seqs, badQuality, qualityThreshold,
    ##                     round(qualityWindow/2))
    ## this step is not necessary at  all
    ## trim if 5 bases are below '0'(fred score 15) in a window of 10 bases
    ## trimmed <- trimTailw(seqs, 5, '+', 5)
    ## trimmed <- trimTailw(seqs, 5, '#', 5)
    ## this step is necessary because many shortreads functions work on ACGT only
    ##trimmed <- trimmed[width(trimmed) > 65]
    trimmed <- seqs
    trimmed <- trimmed[!grepl('N', sread(trimmed))]
    if(length(trimmed) > 0){
      trimmedSeqs <- sread(trimmed)
      trimmedqSeqs <- quality(quality(trimmed))
      names(trimmedSeqs) <- names(trimmedqSeqs) <- 
        sapply(sub("(.+) .+","\\1",ShortRead::id(trimmed)),
               function(z){paste0(alias, "%", strsplit(z, "-")[[1]][2])})
    }
  }
  list(trimmedSeqs, trimmedqSeqs)
})


reads <- sapply(r, "[[", 1)
qualities <- sapply(r, "[[", 2)
#this is needed for primerID quality scores later on
R1Quality <- qualities[[1]]

print(t(stats.bore), quote=FALSE)

primer <- completeMetadata[1,]$primer
ltrbit <- completeMetadata[1,]$ltrBit

PairwiseAlignmentsSingleSubject2DF <- function(PA, shift=0) {
  stopifnot("PairwiseAlignmentsSingleSubject"  %in% class(PA))
  
  return(data.frame(
    width=width(pattern(PA)),
    score=score(PA),
    mismatch=width(pattern(PA))-score(PA),
    start=start(pattern(PA))+shift,
    end=end(pattern(PA))+shift
  ))
}


trim_Ltr_side_reads <- function(reads.p, primer, ltrbit, maxMisMatch=0) {
  
  stopifnot(class(reads.p) %in% "DNAStringSet")
  stopifnot(!any(duplicated(names(reads.p))))
  stopifnot(length(primer)==1)
  stopifnot(length(ltrbit)==1)
  
  ## allows gap, and del/ins count as 1 mismatch
  submat1 <- nucleotideSubstitutionMatrix(match=1,
                                          mismatch=0,
                                          baseOnly=TRUE)
  
  ## p for primer
  ## search for primer from the beginning
  aln.p <- pairwiseAlignment(pattern=subseq(reads.p, 1, 1+nchar(primer)),
                             subject=primer,
                             substitutionMatrix=submat1,
                             gapOpening = 0,
                             gapExtension = 1,
                             type="overlap")
  aln.p.df <- PairwiseAlignmentsSingleSubject2DF(aln.p)
  
  ## l for ltrbit
  ## search for ltrbit fellowing primer
  ## note, for SCID trial, there are GGG between primer and ltr bit and hence 5
  ## for extra bases
  aln.l <- pairwiseAlignment(pattern=subseq(reads.p, nchar(primer)+1, nchar(primer)+nchar(ltrbit)+1),
                             subject=ltrbit,
                             substitutionMatrix=submat1,
                             gapOpening = 0,
                             gapExtension = 1,
                             type="overlap")
  aln.l.df <- PairwiseAlignmentsSingleSubject2DF(aln.l, shift=nchar(primer))
  
  goodIdx <- (aln.p.df$score >= nchar(primer)-maxMisMatch &
                aln.l.df$score >= nchar(ltrbit)-maxMisMatch)
  
  reads.p <- subseq(reads.p[goodIdx], aln.l.df$end[goodIdx]+1)
  
  return(reads.p)
}


reads.p <- trim_Ltr_side_reads(reads[[2]], primer, ltrbit)

linker <- completeMetadata[1,]$linkerSequence

trim_primerIDlinker_side_reads <- function(reads.l, linker, maxMisMatch=3) {
  
  stopifnot(class(reads.l) %in% "DNAStringSet")
  stopifnot(!any(duplicated(names(reads.l))))
  stopifnot(length(linker)==1)
  
  pos.N <- unlist(gregexpr("N", linker))
  len.N <- length(pos.N)
  link1 <- substr(linker, 1, min(pos.N)-1)
  link2 <- substr(linker, max(pos.N)+1, nchar(linker))
  
  ## allows gap, and del/ins count as 1 mismatch
  submat1 <- nucleotideSubstitutionMatrix(match=1,
                                          mismatch=0,
                                          baseOnly=TRUE)
  
  ## search at the beginning for 1st part of linker
  aln.1 <- pairwiseAlignment(pattern=subseq(reads.l, 1, 2+nchar(link1)),
                             subject=link1,
                             substitutionMatrix=submat1,
                             gapOpening = 0,
                             gapExtension = 1,
                             type="overlap")
  aln.1.df <- PairwiseAlignmentsSingleSubject2DF(aln.1)
  
  ## search after 1st part of linker for the 2nd part of linker
  aln.2 <- pairwiseAlignment(pattern=subseq(reads.l, max(pos.N)-1, nchar(linker)+1),
                             subject=link2,
                             substitutionMatrix=submat1,
                             gapOpening = 0,
                             gapExtension = 1,
                             type="overlap")
  aln.2.df <- PairwiseAlignmentsSingleSubject2DF(aln.2, max(pos.N)-2)
  
  goodIdx <- (aln.1.df$score >= nchar(link1)-maxMisMatch &
                aln.2.df$score >= nchar(link2)-maxMisMatch)
  
  primerID <- subseq(reads.l[goodIdx],
                     aln.1.df$end[goodIdx]+1,
                     aln.2.df$start[goodIdx]-1)
  
  reads.l <- subseq(reads.l[goodIdx], aln.2.df$end[goodIdx]+1)
  
  stopifnot(all(names(primerID)==names(reads.l)))
  
  return(list("reads.l"=reads.l,
              "primerID"=primerID))
}

readslprimer <- trim_primerIDlinker_side_reads(reads[[1]], linker)

reads.l <- readslprimer$reads.l
primerIDs <- readslprimer$primerID
stats.bore$linkered <- length(reads.l)

ltrlinkeredQname <- intersect(names(reads.p), names(reads.l))

reads.p <- reads.p[ltrlinkeredQname]
reads.l <- reads.l[ltrlinkeredQname]
stats.bore$ltredlinkered <- length(reads.l)

print(t(stats.bore), quote=FALSE) 

trim_overreading <- function(reads, marker, maxMisMatch=3) {
  
  stopifnot(class(reads) %in% "DNAStringSet")
  stopifnot(!any(duplicated(names(reads))))
  stopifnot(length(marker)==1)
  
  
  submat1 <- nucleotideSubstitutionMatrix(match=1,
                                          mismatch=0,
                                          baseOnly=TRUE)
  
  ## allows gap, and del/ins count as 1 mismatch
  tmp <- pairwiseAlignment(pattern=reads,
                           subject=marker,
                           substitutionMatrix=submat1,
                           gapOpening = 0,
                           gapExtension = 1,
                           type="overlap")
  
  odf <- PairwiseAlignmentsSingleSubject2DF(tmp)
  
  odf$isgood <- FALSE
  ## overlap in the middle or at right
  odf$isgood <- with(odf, ifelse(mismatch<=maxMisMatch &
                                   start>1,
                                 TRUE, isgood))
  
  ## overlap at left
  odf$isgood <- with(odf, ifelse(mismatch<=maxMisMatch &
                                   start==1 &
                                   width>=nchar(marker)-1,
                                 TRUE, isgood))
  
  ## note with ovelrap alignmment, it only align with a minimum of 1/2 of the shorter one
  odf$cut <- with(odf, ifelse(isgood, odf$start-1, nchar(reads)))
  if( any(odf$cut < nchar(reads)) ) {
    odf$cut <- nchar(reads)-nchar(marker)/2
    odf$cut <- with(odf, ifelse(isgood, odf$start-1, cut))
  }
  
  reads <- subseq(reads, 1, odf$cut)
}

linker_common <- completeMetadata[1,]$linkerCommon

reads.p.o <- trim_overreading(reads.p, linker_common, 3)

largeLTRFrag <- completeMetadata[1,]$largeLTRFrag

reads.l.o <- trim_overreading(reads.l, substr(largeLTRFrag, 1, 20), 3)

mingDNA <- completeMetadata[1,]$mingDNA

# select sequence length > 30 for R1 and R2 
reads.pp <- reads.p.o[which(width(reads.p.o) > mingDNA)]

reads.ll <- reads.l.o[which(width(reads.l.o) > mingDNA)]

# select the reads have common names in R1 and R2

reads.p <- reads.pp
reads.l <- reads.ll

ltrlinkeredQname <- intersect(names(reads.p), names(reads.l))
reads.p <- reads.p[ltrlinkeredQname]
reads.l <- reads.l[ltrlinkeredQname]

stats.bore$lenTrim <- length(reads.p)

readpsl <- function(pslFile, toNull=NULL) {
  cols <- c("matches", "misMatches", "repMatches", "nCount", "qNumInsert",
            "qBaseInsert", "tNumInsert", "tBaseInsert", "strand", "qName",
            "qSize", "qStart", "qEnd", "tName", "tSize", "tStart", "tEnd",
            "blockCount", "blockSizes", "qStarts", "tStarts")
  cols.class <- c(rep("numeric",8), rep("character",2), rep("numeric",3),
                  "character", rep("numeric",4), rep("character",3))
  
  psl <- lapply(pslFile, function(f) {
    message("Reading ",f)
    data.table::fread( paste("zcat", f), sep="\t" )
  })
  psl <- data.table::rbindlist(psl)
  colnames(psl) <- cols
  
  if(length(toNull)>0) psl[, toNull] <- NULL
  
  return(as.data.frame(psl))
}


findVectorReads <- function(vectorSeq, primerLTR="GAAAATCTCTAGCA",
                            reads.p, reads.l,
                            debug=FALSE) {
  
  library(hiReadsProcessor)
  library(dplyr)
  
  Vector <- readDNAStringSet(vectorSeq)
  
  primerInVector <- matchPattern(pattern=primerLTR,
                                 subject=DNAString(as.character(Vector)),
                                 algorithm="auto",
                                 max.mismatch=4,
                                 with.indels=TRUE,
                                 fixed=TRUE)
  
  #print(primerInVector)
  if( length(primerInVector)<1 ) message("--- Cannot locate primer and ltrBit in vector ---")
  
  
  globalIdentity <- 0.75
  blatParameters <- c(minIdentity=70, minScore=15, stepSize=3, 
                      tileSize=8, repMatch=112312, dots=1000, 
                      q="dna", t="dna", out="psl")
  
  
  hits.v.p <- try(readpsl(blatSeqs(query=reads.p, subject=Vector,     
                                   blatParameters=blatParameters, parallel=F)))
  if( class(hits.v.p) == "try-error" ) hits.v.p <- data.frame()
  if ( debug ) save(hits.v.p, file="hits.v.p.RData")    
  
  hits.v.l <- try(readpsl(blatSeqs(query=reads.l, subject=Vector, 
                                   blatParameters=blatParameters, parallel=F)))
  if( class(hits.v.l) == "try-error" ) hits.v.l <- data.frame()
  if ( debug ) save(hits.v.l, file="hits.v.l.RData")    
  
  ## Sometimes the vector files received from collaborators are different from the 
  ## vector put in human host. So, it is not feasible to put a lot of constrains.
  ## Filtering on globalIdentity identities for both R1 and R2 seems to work well. 
  hits.v.p <- dplyr::filter(hits.v.p, ##tStart  > ltrpos &
                            ##tStart  < ltrpos+nchar(primerLTR)+10 &
                            ##strand == "+" &
                            matches > globalIdentity*qSize &
                              qStart  <= 5 ) 
  hits.v.l <- dplyr::filter(hits.v.l, matches>globalIdentity*qSize )
  ##strand=="-") 
  hits.v <- try(merge(hits.v.p[, c("qName", "tStart")],
                      hits.v.l[, c("qName", "tStart")],
                      by="qName")
                ,silent = TRUE)
  if( class(hits.v) == "try-error" ) hits.v <- data.frame()
  
  ##hits.v <- dplyr::filter(hits.v, tStart.y >= tStart.x &
  ##                                tStart.y <= tStart.x+2000)
  
  if ( debug ) {
    save(reads.p, file="reads.p.RData")
    save(reads.l, file="reads.l.RData")
  }
  
  vqName <- unique(hits.v$qName)
  
  message(paste0("Vector sequences found ", length(vqName)))
  
  return(vqName)
}

#vectorSeq <- '/home/ubuntu/intsitecaller/testCases/intSiteValidation/p746vector.fasta'

primerLTR <- paste0(primer, ltrbit)

vqName <- findVectorReads(vectorSeq,paste0(primer, ltrbit),reads.p, reads.l,debug=TRUE)


print(vqName)

toload <- names(reads.p)[!names(reads.p) %in% vqName]
reads.p <- reads.p[toload]
reads.l <- reads.l[toload]
stats.bore$vTrimed <- length(reads.p)


print(stats.bore)


reads.p.u <- unique(reads.p)
reads.l.u <- unique(reads.l)

print(reads.p)
print(reads.p.u)

print(reads.l)
print(reads.l.u)


reads.p30.u <- unique(subseq(reads.p,1,mingDNA))

stats.bore$uniqL <- length(reads.l.u)  
stats.bore$uniqP <- length(reads.p.u)  
stats.bore$uniqP30 <- length(reads.p30.u)  

names(reads.p.u) <- seq_along(reads.p.u)
names(reads.l.u) <- seq_along(reads.l.u)

print(names(reads.p.u))

print(stats.bore)

keys <- data.frame("R2"=match(reads.p, reads.p.u),
                   "R1"=match(reads.l, reads.l.u),
                   "names"=toload)

keys$readPairKey <- paste0(keys$R1, "_", keys$R2)

print(keys)

saveRDS(keys, file="keys.rds")

chunkSize <- 30000

stats.bore$lLen <- as.integer(mean(width(reads.l)))  
stats.bore$pLen <- as.integer(mean(width(reads.p)))


stats <- data.frame()

stats <- rbind(stats, stats.bore)

save(stats, file="stats.RData")

if(length(toload) > 0){
  
  chunks.p <- split(seq_along(reads.p.u), ceiling(seq_along(reads.p.u)/chunkSize))
  for(i in c(1:length(chunks.p))){
    writeXStringSet(reads.p.u[chunks.p[[i]]],
                    file=paste0("R2-", i, ".fa"),
                    append=FALSE)
  }
  
  chunks.l <- split(seq_along(reads.l.u), ceiling(seq_along(reads.l.u)/chunkSize))
  for(i in c(1:length(chunks.l))){    
    writeXStringSet(reads.l.u[chunks.l[[i]]],
                    file=paste0("R1-", i, ".fa"),
                    append=FALSE)
  }
  
  save(stats, file="stats.RData")
  alias #return 'value' which ultimately gets saved as trimStatus.RData
}else{
  stop("error - no curated reads")
}

print(stats)

#ref.bit <- '/home/ubuntu/intsitecaller/testCases/intSiteValidation/hg18.2bit'

cmd <- paste0('blat ', ref.bit,' ','R1-1.fa',' ','R1-1.fa.psl',' ','-tileSize=11 -stepSize=9 -minIdentity=85 -maxIntron=5 -minScore=27 -dots=1000 -out=psl -noHead')

system(cmd)

cmd <- paste0('blat ', ref.bit,' ','R2-1.fa',' ','R2-1.fa.psl',' ','-tileSize=11 -stepSize=9 -minIdentity=85 -maxIntron=5 -minScore=27 -dots=1000 -out=psl -noHead')

system(cmd)

    
#Rscript ~/intsitecaller/PslToIs.R R1-1.fa.psl R2-1.fa.psl keys.rds ~/intsitecaller/testCases/intSiteValidation/completeMetadata.RData /home/ubuntu/intsitecaller/clone1-1










