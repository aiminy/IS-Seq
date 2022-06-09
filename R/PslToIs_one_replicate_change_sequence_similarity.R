# This script use blat-aligned R1 and R2 psl file to call IS

# Rscript /home/ubuntu/intsitecaller/PslToIs_one_replicate_change_sequence_similarity.R /home/ubuntu/SHARE/ISseqOutput/Feb9G222/IsaByINSPIIRED/fa/HL60cl60HL60Poly100/R2_fastq_trim12nt_qcTrimmed_MatchBlastLtrLc_Barcode_FB-P7-Rd2-LC.20.fq_trimwithCutAdapt_HL60cl60HL60Poly100_ReadyToAlignSort.fa.psl /home/ubuntu/SHARE/ISseqOutput/Feb9G222/IsaByINSPIIRED/fa/HL60cl60HL60Poly100/R1_fastq_trim12nt_qcTrimmed_MatchBlastLtrLc_Barcode_FB-P5-Rd1-LTR.16.fq_trimwithCutAdapt_HL60cl60HL60Poly100_ReadyToAlignSort.fa.psl /home/ubuntu/SHARE/ISseqOutput/Feb9G222/IsaByINSPIIRED/fa/HL60cl60HL60Poly100/keys.rds ~/intsitecaller/testCases/intSiteValidation/completeMetadata.RData /home/ubuntu/SHARE/ISseqOutput/Feb9G222/IsaByINSPIIRED/fa/HL60cl60HL60Poly100/rev0 hg38 1 0

# Rscript ~/IS-Seq/R/PslToIs_one_replicate_change_sequence_similarity.R /home/ubuntu/SHARE/Aimin/INSPIIRED_test_output/clone1-1/R1-1.fa.psl /home/ubuntu/SHARE/Aimin/INSPIIRED_test_output/clone1-1/R2-1.fa.psl /home/ubuntu/SHARE/Aimin/INSPIIRED_test_output/clone1-1/keys.rds ~/IS-Seq/utilsRefData/INSPIIRED/completeMetadata.RData ~/SHARE/Aimin/INSPIIRED_test_output/clone1-1/IS_95 hg18 1 95

initial.options <- commandArgs(trailingOnly = FALSE)
print(initial.options)

file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.dirname <- dirname(script.name)
print(script.dirname)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)>=1) {
  # default output file
  input_R1_psl=args[1]
  input_R2_psl=args[2]
  input_keys=args[3]
  input_completeMetadata=args[4]
  output.dir=args[5]
  ref_genome=args[6]
  index=args[7]
  minPctIdent=args[8]
}

print(input_R1_psl)
print(input_R2_psl)
print(input_keys)
print(input_completeMetadata)
print(output.dir)

if(!dir.exists(output.dir)){dir.create(output.dir,recursive = TRUE)}

#aws_root <- '/home/ubuntu/SHARE/ISseqOutput/Oct2721/IsaByINSPIIRED'
#output.dir <- '/home/ubuntu/SHARE/ISseqOutput/Oct2721/IsaByINSPIIRED/output'

#if(!dir.exists(output.dir)){dir.create(output.dir,recursive = TRUE)}

#input.file <- '~/intsitecaller/testCases/intSiteValidation/completeMetadata.RData'

completeMetadata <- get(load(input_completeMetadata))

#psl.R2 <- list.files(file.path(aws_root,'align'), pattern="R2.*.fa.psl",recursive = T,full.names=T)

psl.R2 <- input_R2_psl 

#psl.R1 <- list.files(file.path(aws_root,'align'), pattern="R1.*.fa.psl",recursive = T,full.names=T)

psl.R1 <- input_R1_psl 

print(psl.R2)

print(psl.R1)

#source('~/ispipe/R/functions.R')

#' read psl gz files, assuming psl gz files don't have column header
#' @param pslFile character vector of file name(s)
#' @param toNull  character vector of column names to get rid of
#' @return data.frame, data.table of the psl table
#' @example 
readpsl <- function(pslFile, toNull=NULL) {
  stopifnot(require("data.table"))
  cols <- c("matches", "misMatches", "repMatches", "nCount", "qNumInsert",
            "qBaseInsert", "tNumInsert", "tBaseInsert", "strand", "qName",
            "qSize", "qStart", "qEnd", "tName", "tSize", "tStart", "tEnd",
            "blockCount", "blockSizes", "qStarts", "tStarts")
  cols.class <- c(rep("numeric",8), rep("character",2), rep("numeric",3),
                  "character", rep("numeric",4), rep("character",3))
  
  psl <- lapply(pslFile, function(f) {
    message("Reading ",f)
    data.table::fread(f, sep="\t" )
  })
  psl <- data.table::rbindlist(psl)
  colnames(psl) <- cols
  
  if(length(toNull)>0) psl[, toNull] <- NULL
  
  return(as.data.frame(psl))
}

hits.R2 <- readpsl(psl.R2)

hits.R1 <- readpsl(psl.R1)

keys <- readRDS(input_keys)

readsAligning <- length(which(keys$R2 %in% hits.R2$qName & keys$R1 %in% hits.R1$qName))

stats <- data.frame(readsAligning=readsAligning)

print(stats)

source(file.path(script.dirname,'quality_filter.R'))
source(file.path(script.dirname,'processBLATData.R'))

u=index

completeMetadata$refGenome <- ref_genome

if(completeMetadata[u,]$refGenome=='hg38'){
  library(BSgenome.Hsapiens.UCSC.hg38)
}else
{
  library(BSgenome.Hsapiens.UCSC.hg18)
}

print(hits.R2)

print(hits.R1)

save.image(file.path(output.dir,'temp0.RData'))

if(u==1){completeMetadata[u,]$minPctIdent <- minPctIdent}

hits.R2 <- qualityFilter(hits.R2, completeMetadata[u,]$maxAlignStart, completeMetadata[u,]$minPctIdent)

print(hits.R2)

hits.R2 <- processBLATData(hits.R2, "R2", completeMetadata[u,]$refGenome)  

#print(hits.R2)

#print(hits.R1)

save.image(file.path(output.dir,'temp.RData'))

#if(u==1){completeMetadata[u,]$minPctIdent <- 0}

hits.R1 <- qualityFilter(hits.R1, completeMetadata[u,]$maxAlignStart, completeMetadata[u,]$minPctIdent)

print(hits.R1)

hits.R1 <- processBLATData(hits.R1, "R1", completeMetadata[u,]$refGenome)
  
stopifnot(all(strand(hits.R1) == "+" | strand(hits.R1) == "-"))
stopifnot(all(strand(hits.R2) == "+" | strand(hits.R2) == "-"))
  
readsWithGoodAlgnmts <- length(which(keys$R2 %in% hits.R2$qName & keys$R1 %in% hits.R1$qName))

stats <- cbind(stats, readsWithGoodAlgnmts)
print(stats)
unique_key_pairs <- unique(keys[,c("R1", "R2", "readPairKey")])

library(GenomicRanges)
  
  red.hits.R1 <- GenomicRanges::reduce(GenomicRanges::flank(hits.R1, -1, start = TRUE), min.gapwidth = 0L, with.revmap = TRUE)
  
  red.hits.R2 <- GenomicRanges::reduce(
    GenomicRanges::flank(hits.R2, -1, start = TRUE), min.gapwidth = 0L, with.revmap = TRUE)
  
  pairs <- GenomicRanges::findOverlaps(
    red.hits.R1,
    red.hits.R2,
    maxgap = completeMetadata[u,]$maxFragLength,
    ignore.strand = TRUE
  )
  
  R1.loci <- red.hits.R1[queryHits(pairs)]
  R2.loci <- red.hits.R2[subjectHits(pairs)]
  
  #' Check isDownstream and isOppositeStrand
  R1.loci.starts <- BiocGenerics::start(R1.loci)
  R2.loci.starts <- BiocGenerics::start(R2.loci)
  
  R1.loci.strand <- BiocGenerics::strand(R1.loci)
  R2.loci.strand <- BiocGenerics::strand(R2.loci)

 keep.loci <- ifelse(
    R2.loci.strand == "+", 
    as.vector(R1.loci.starts > R2.loci.starts & 
                R1.loci.strand != R2.loci.strand), 
    as.vector(R1.loci.starts < R2.loci.starts & 
                R1.loci.strand != R2.loci.strand))
  
  keep.loci <- as.vector(
    keep.loci & R2.loci.strand != "*" & R1.loci.strand != "*")
  
  R1.loci <- R1.loci[keep.loci]
  R2.loci <- R2.loci[keep.loci]

library(IRanges)

  loci.key <- data.frame(
    "R1.loci" = queryHits(pairs)[keep.loci],
    "R2.loci" = subjectHits(pairs)[keep.loci])
  loci.key$lociPairKey <- paste0(loci.key$R1.loci, ":", loci.key$R2.loci)
  
  R1.mapping.index <- lapply(R1.loci$revmap, function(x){
    as.integer(hits.R1$qName[x])
  })
  
  R1.mapping.index.1 <- IRanges::IntegerList(R1.mapping.index)
  
  loci.key$R1.qNames <- IRanges::IntegerList(lapply(R1.loci$revmap, function(x){
    as.integer(hits.R1$qName[x])
  }))
  
  loci.key$R2.qNames <- IRanges::IntegerList(lapply(R2.loci$revmap, function(x){
    as.integer(hits.R2$qName[x])
  }))
  
  loci.key$R1.readPairs <- IRanges::IntegerList(lapply(
    loci.key$R1.qNames, function(x){
      which(unique_key_pairs$R1 %in% x)
    }))
  
  loci.key$R2.readPairs <- IRanges::IntegerList(lapply(
    loci.key$R2.qNames, function(x){
      which(unique_key_pairs$R2 %in% x)
    }))
  
  paired.loci <- GRanges(
    seqnames = seqnames(R2.loci), 
    ranges = IRanges(
      start = ifelse(strand(R2.loci) == "+", start(R2.loci), start(R1.loci)),
      end = ifelse(strand(R2.loci) == "+", end(R1.loci), end(R2.loci))),
    strand = strand(R2.loci),
    lociPairKey = loci.key$lociPairKey)
  
  paired.loci$readPairKeys <- CharacterList(lapply(
    1:length(paired.loci), 
    function(i){
      unique_key_pairs[intersect(
        loci.key$R1.readPairs[[i]], 
        loci.key$R2.readPairs[[i]]),
        "readPairKey"]
    }))
  paired.loci <- paired.loci[sapply(paired.loci$readPairKeys, length) > 0]
  
  read.loci.mat <- data.frame(
    "lociPairKey" = Rle(
      values = paired.loci$lociPairKey,
      lengths = sapply(paired.loci$readPairKeys, length)),
    "readPairKey" = unlist(paired.loci$readPairKeys)
  )
  
  numProperlyPairedAlignments <- nrow(
    keys[keys$readPairKey %in% read.loci.mat$readPairKey,])
  
  stats <- cbind(stats, numProperlyPairedAlignments)

  readPairCounts <- table(read.loci.mat$readPairKey)
  uniq.readPairs <- names(readPairCounts[readPairCounts == 1])
  multihit.readPairs <- names(readPairCounts[readPairCounts > 1])


  saveRDS(read.loci.mat,file.path(output.dir,'read.loci.mat.rds'))
  print(readPairCounts)


failedReads <- keys[!keys$readPairKey %in% read.loci.mat$readPairKey,]

chimera.reads <- failedReads[failedReads$R1 %in% hits.R1$qName & failedReads$R2 %in% hits.R2$qName,]
  
  if(dim(chimera.reads)[1]!=0){
    chimera.alignments <- GRangesList(lapply(1:dim(chimera.reads)[1], function(i){
      R1 <- hits.R1[hits.R1$qName == chimera.reads[i, "R1"]]
      R2 <- hits.R2[hits.R2$qName == chimera.reads[i, "R2"]]
      names(R1) <- rep(chimera.reads[i, "names"], length(R1))
      names(R2) <- rep(chimera.reads[i, "names"], length(R2))
      c(R2, R1)
    }))
    chimeraData <- list("read_info" = chimera.reads, "alignments" = chimera.alignments)

    save(chimeraData, file = file.path(output.dir,"chimeraData.RData"))

  #' Record chimera metrics
  chimeras <- length(unique(chimera.reads$names))
  stats <- cbind(stats, chimeras)
  }
  
 
#chimeras <- length(unique(chimera.reads$names))

#stats <- cbind(stats, chimeras)
 
  uniq.read.loci.mat <- read.loci.mat[
    read.loci.mat$readPairKey %in% uniq.readPairs,]
  
  uniq.templates <- paired.loci[
    match(uniq.read.loci.mat$lociPairKey, paired.loci$lociPairKey)]
  uniq.templates$readPairKeys <- NULL
  uniq.templates$readPairKey <- uniq.read.loci.mat$readPairKey
  
  uniq.keys <- keys[keys$readPairKey %in% uniq.readPairs,]
  uniq.reads <- uniq.templates[
    match(uniq.keys$readPairKey, uniq.templates$readPairKey)]
  names(uniq.reads) <- as.character(uniq.keys$names)
  uniq.reads$sampleName <- sapply(
    strsplit(as.character(uniq.keys$names), "%"), "[[", 1)
  uniq.reads$ID <- sapply(strsplit(as.character(uniq.keys$names), "%"), "[[", 2)
  
  allSites <- uniq.reads

  saveRDS(allSites, file=file.path(output.dir,"allSites.rds"))


dereplicateSites <- function(sites){
  #Reduce sites which have the same starts, but loose range info
  #(no need to add a gapwidth as sites are standardized)
  sites.reduced <- flank(sites, -1, start=TRUE)
  sites.reduced <- reduce(sites.reduced, min.gapwidth = 0L, with.revmap=TRUE)
  
  sites.reduced$counts <- sapply(sites.reduced$revmap, length)
  
  #Order original sites by revmap  
  dereplicatedSites <- sites[unlist(sites.reduced$revmap)]
  
  #Skip this step and provide similar output if length(sites) = 0
  if(length(sites) > 0){
    dereplicatedSites <- split(dereplicatedSites, Rle(values = seq(length(sites.reduced)), lengths = sites.reduced$counts))
  }  
  
  #Dereplicate reads with same standardized starts and provide the longeset width
  dereplicatedSites <- unlist(reduce(dereplicatedSites, min.gapwidth = 0L))
  mcols(dereplicatedSites) <- mcols(sites.reduced)
  
  dereplicatedSites
  
}


 sites.final <- dereplicateSites(allSites)
  if(length(sites.final)>0){
    sites.final$sampleName <- allSites[1]$sampleName
    sites.final$posid <- paste0(as.character(seqnames(sites.final)),
                                as.character(strand(sites.final)),
                                start(flank(sites.final, width=-1, start=TRUE)))
  }
  saveRDS(sites.final, file=file.path(output.dir,"sites.final.rds"))


library(igraph)

numAllSingleReads <- length(allSites)
  stats <- cbind(stats,numAllSingleReads)
  
  numAllSingleSonicLengths <- length(unique(granges(allSites)))
  stats <- cbind(stats, numAllSingleSonicLengths)

  numUniqueSites <- length(sites.final)
  stats <- cbind(stats, numUniqueSites)
 
  save.image(file=file.path(output.dir,'IS_call_before_clean.RData'))   

  #' Clean up environment for expansion and clustering of multihits
  rm(uniq.read.loci.mat, uniq.templates, uniq.keys, 
     uniq.reads, allSites, sites.final)
  gc()


   #' ########## IDENTIFY MULTIPLY-PAIRED READS (multihits) ##########
  #' Multihits are reads that align to multiple locations in the reference
  #' genome. There are bound to always be a certain proportion of reads aligning
  #' to repeated sequence due to the high level degree of repeated DNA elements
  #' within genomes. The final object generated, "multihitData", is a list of
  #' three objects. "unclusteredMultihits" is a GRanges object where every
  #' alignment for every multihit read is present in rows.
  #' "clusteredMultihitPositions" returns all the possible integration site
  #' positions for the multihit. Lastly, "clusteredMultihitLengths" contains the
  #' length of the templates mapping to the multihit clusters, used for
  #' abundance calculations.
  unclusteredMultihits <- GRanges()
  clusteredMultihitPositions <- GRangesList()
  clusteredMultihitLengths <- list()

  if(length(multihit.readPairs) > 0){
    #' Only consider readPairKeys that aligned to multiple genomic loci
    multi.read.loci.mat <- read.loci.mat[
      read.loci.mat$readPairKey %in% multihit.readPairs,]

    multihit.templates <- paired.loci[
      paired.loci$lociPairKey %in% multi.read.loci.mat$lociPairKey]
    multihit.expansion.map <- multihit.templates$readPairKeys
    multihit.templates$readPairKeys <- NULL
    multihit.templates <- multihit.templates[Rle(
      values = 1:length(multihit.templates),
      lengths = sapply(multihit.expansion.map, length)
    )]
    multihit.templates$readPairKey <- unlist(multihit.expansion.map)

    #' As the loci are expanded from the paired.loci object, unique templates
    #' and readPairKeys are present in the readPairKeys unlisted from the
    #' paired.loci object.
    multihit.templates <- multihit.templates[
      multihit.templates$readPairKey %in% multi.read.loci.mat$readPairKey]

    multihit.keys <- keys[keys$readPairKey %in% multihit.readPairs,]
    multihit.keys$sampleName <- sapply(strsplit(
      as.character(multihit.keys$names), "%"), "[[", 1)
    multihit.keys$ID <- sapply(strsplit(
      as.character(multihit.keys$names), "%"), "[[", 2)

    #' Medians are based on all the potential sites for a given read, which will
    #' be identical for all reads associated with a readPairKey.
    multihit.medians <- round(
      median(width(split(multihit.templates, multihit.templates$readPairKey))))
    multihit.keys$medians <- multihit.medians[multihit.keys$readPairKey]

    multihits.pos <- flank(multihit.templates, -1, start = TRUE)
    multihits.red <- reduce(multihits.pos, min.gapwidth = 5L, with.revmap = TRUE)  #! Should make 5L a option
    revmap <- multihits.red$revmap

    axil_nodes <- as.character(Rle(
      values = multihit.templates$readPairKey[sapply(revmap, "[", 1)],
      lengths = sapply(revmap, length)
    ))
    nodes <- multihit.templates$readPairKey[unlist(revmap)]
    edgelist <- unique(matrix( c(axil_nodes, nodes), ncol = 2 ))
    multihits.clusterData <- igraph::clusters(
      igraph::graph.edgelist(edgelist, directed=F))
    clus.key <- data.frame(
      row.names = unique(as.character(t(edgelist))),
      "clusID" = multihits.clusterData$membership)

    multihits.pos$clusID <- clus.key[multihits.pos$readPairKey, "clusID"]
    clusteredMultihitPositions <- split(multihits.pos, multihits.pos$clusID)
    clusteredMultihitNames <- lapply(
      clusteredMultihitPositions, function(x) unique(x$readPairKey))
    clusteredMultihitPositions <- GRangesList(lapply(
      clusteredMultihitPositions,
      function(x){
        unname(unique(granges(x)))
      }))

    clusteredMultihitLengths <- lapply(clusteredMultihitNames, function(x){
      readIDs <- unique(multihit.keys[multihit.keys$readPairKey %in% x,]$ID)
      data.frame(table(multihit.keys[multihit.keys$ID %in% readIDs,]$medians))
    })

    #' Expand the multihit.templates object from readPairKey specific to read
    #' specific.
    multihit.keys <- multihit.keys[order(multihit.keys$readPairKey),]
    multihit.readPair.read.exp <- IntegerList(lapply(
      unique(multihit.keys$readPairKey),
      function(x){which(multihit.keys$readPairKey == x)}))
    names(multihit.readPair.read.exp) <- unique(multihit.keys$readPairKey)
    unclusteredMultihits <- multihit.templates
    multihit.readPair.read.exp <- as(multihit.readPair.read.exp[
      as.character(unclusteredMultihits$readPairKey)], "SimpleList")
    unclusteredMultihits <- unclusteredMultihits[Rle(
      values = 1:length(unclusteredMultihits),
      lengths = sapply(multihit.readPair.read.exp, length)
    )]
    names(unclusteredMultihits) <- multihit.keys$names[
      unlist(multihit.readPair.read.exp)]
    unclusteredMultihits$ID <- multihit.keys$ID[
      unlist(multihit.readPair.read.exp)]
    unclusteredMultihits$sampleName <- multihit.keys$sampleName[
      unlist(multihit.readPair.read.exp)]
  }

  stopifnot(length(clusteredMultihitPositions)==length(clusteredMultihitLengths))
  multihitData <- list(unclusteredMultihits, clusteredMultihitPositions, clusteredMultihitLengths)
  names(multihitData) <- c("unclusteredMultihits", "clusteredMultihitPositions", "clusteredMultihitLengths")

  saveRDS(multihitData, file=file.path(output.dir,"multihitData.rds"))



  #' Record multihit metrics (reads, clusters, sonicLengths)
  multihitReads <- nrow(keys[keys$readPairKey %in% multihit.readPairs,])
  stats <- cbind(stats, multihitReads)

  multihitSonicLengths <- 0
  if( length(multihitData$clusteredMultihitLengths) > 0 ) {
    multihitSonicLengths <- sum(
      sapply(multihitData$clusteredMultihitLengths, nrow))
  }
  stats <- cbind(stats, multihitSonicLengths)

  multihitClusters <- length(multihitData$clusteredMultihitPositions)
  stats <- cbind(stats, multihitClusters)

  #' Finalize metrics by combining unique alignments and multihit clusters
  totalSonicLengths <- numAllSingleSonicLengths + multihitSonicLengths
  stats <- cbind(stats, totalSonicLengths)

  totalEvents <- numUniqueSites + multihitClusters
  stats <- cbind(stats, totalEvents)

  print(stats)

  saveRDS(stats, file=file.path(output.dir,"stats.rds"))


save.image(file=file.path(output.dir,'IS_call.RData'))
