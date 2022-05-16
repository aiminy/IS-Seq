get_reference_genome <- function(reference_genome) {
  pattern <- paste0("\\.", reference_genome, "$")
  match_index <- which(grepl(pattern, installed.genomes()))
  if (length(match_index) != 1) {
    write("Installed genomes are:", stderr())
    write(installed.genomes(), stderr())
    stop(paste("Cannot find unique genome for", reference_genome))
  }
  BS_genome_full_name <- installed.genomes()[match_index]
  library(BS_genome_full_name, character.only=T)
  get(BS_genome_full_name)
}

processBLATData <- function(algns, from, refGenome){
  stopifnot(from == "R1" | from == "R2")
  algns$from <- from
  algns$qtStart <- ifelse(
    algns$strand == "+",
    (algns$tStart - (algns$qStart)),
    (algns$tStart - (algns$qSize - algns$qEnd - 1)))
  algns$qtEnd <- ifelse(
    algns$strand == "+",
    (algns$tEnd + (algns$qSize - algns$qEnd - 1)),
    (algns$tEnd + (algns$qStart)))    
  
  algns.gr <- GRanges(seqnames=Rle(algns$tName),
                      ranges = IRanges(
                        start = (algns$qtStart + 1), 
                        end = (algns$qtEnd)), #Convert to 1-base
                      strand=Rle(algns$strand),
                      seqinfo=seqinfo(get_reference_genome(refGenome)))
  
  mcols(algns.gr) <- algns[,c("from", "qName", "matches", "repMatches", 
                              "misMatches", "qStart", "qEnd", "qSize", 
                              "tBaseInsert")]
  rm(algns)
  algns.gr
}
