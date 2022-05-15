# Inputs:
# /home/ubuntu/intsitecaller/testCases/intSiteValidation/Data/Undetermined_S0_L001_I1_001.fastq.gz

# /home/ubuntu/intsitecaller/testCases/intSiteValidation/completeMetadata.RData

# /home/ubuntu/intsitecaller/testCases/intSiteValidation/Data/Undetermined_S0_L001_R1_001.fastq.gz

# /home/ubuntu/intsitecaller/testCases/intSiteValidation/Data/Undetermined_S0_L001_R2_001.fastq.gz

# Outputs:

# trimmed and corrected index files
# MyTest/correctedI1-1.fasta  
# MyTest/trimmedI1-1.fasta

# demultiplexed reads
# MyTest/demultiplexedReps/*fastq.gz

library("ShortRead")

#  codeDir <- get(load("codeDir.RData"))
#  completeMetadata <- get(load("completeMetadata.RData"))
#  jobID <- get(load("jobID.RData"))

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
	  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)>=1) {
	  # default output file
      output.dir=args[1]
}

  if(!dir.exists(output.dir)){dir.create(output.dir,recursive = TRUE)}

  I1 <- readFastq("/home/ubuntu/intsitecaller/testCases/intSiteValidation/Data/Undetermined_S0_L001_I1_001.fastq.gz")

  print(I1)


  I1 <- trimTailw(I1, 2, "0", 12)
  I1 <- I1[width(I1)==max(width(I1))]
  I1 <- split(I1, ceiling(seq_along(I1)/500000))



  for(chunk in names(I1)){

      output_fasta=file.path(output.dir,paste0("trimmedI1-", chunk, ".fasta"))
      writeFasta(I1[[chunk]], file=output_fasta)
      cmd <- paste0('export PYTHONPATH=/home/ubuntu/miniconda2/lib/python2.7/site-packages;python processGolayTest.py ',output_fasta,' ',file.path(output.dir,'correctedI1-1.fasta'),' ',file.path(output.dir,'correctedI1-1.done'))
      system(cmd)

  }
  
 
demultiplex_reads <- function(reads, suffix, I1Names, samples, completeMetadata,output.dir) {

RNames <- sapply(strsplit(as.character(ShortRead::id(reads)), " "), "[[", 1)#for some reason we can't dynamically set name/id on ShortRead!
names(RNames) <- NULL

reads <- reads[match(I1Names, RNames)]
reads <- split(reads, samples)
	    
for (i in 1:length(reads)){

writeLog(paste0('Demultiplexing ', suffix, ' read: ', i, '/', length(reads)))

barcode.i <- completeMetadata$bcSeq[ completeMetadata$alias == names(reads)[i] ]
stopifnot(length(barcode.i)==1)
		        
alias_by_barcode <- completeMetadata$alias[ completeMetadata$bcSeq == barcode.i ]
stopifnot(length(alias_by_barcode)>=1)

fqFiles <- file.path(output.dir,paste0("demultiplexedReps/", alias_by_barcode, "_", suffix, ".fastq.gz"))
				        
cat(barcode.i, "\t", paste(fqFiles, collapse=" "), "\n" )

null <- sapply(fqFiles, function(fq) writeFastq(reads[[i]], fq, mode="w") )
}  
}

writeLog <- function(...)
{
   arguments <- list(...)
  
   for(i in arguments){
	   if ( typeof(i) == "character" ){
		write(i,file='intSiteCaller.log',append=T)
         }else{
		 w <- try(write.table(i,file='intSiteCaller.log',append=T,sep="\t",quote=F))
		 if (class(w) == "try-error"){
			 write.table("Could not write requested data item\n",file='intSiteCaller.log',append=T,sep="\t",quote=F)}
	 }
   }
}

 I1 <- readFasta(list.files(output.dir, pattern="correctedI1-.", full.names=T))
  
 completeMetadata <- get(load("./testCases/intSiteValidation/completeMetadata.RData"))
  
 I1 <- I1[as.vector(sread(I1)) %in% completeMetadata$bcSeq]
 
 samples <- completeMetadata[match(as.character(sread(I1)), completeMetadata$bcSeq), "alias"]
    
  #only necessary if using native data - can parse out description w/ python
 I1Names <-  sapply(strsplit(as.character(ShortRead::id(I1)), " "), "[[", 1)#for some reason we can't dynamically set name/id on ShortRead!
      
  unlink(file.path(output.dir,"demultiplexedReps"), recursive=TRUE,  force=TRUE)
  suppressWarnings(dir.create(file.path(output.dir,"demultiplexedReps")))

  writeLog('Starting to demultiplex R1')
 R1 <- readFastq("./testCases/intSiteValidation/Data/Undetermined_S0_L001_R1_001.fastq.gz")
  
  demultiplex_reads(R1, "R1", I1Names, samples, completeMetadata,output.dir)
  
  writeLog('completed demultiplexing R1')
  
  writeLog('Starting to demultiplex R2')
  
 R2 <- readFastq("./testCases/intSiteValidation/Data/Undetermined_S0_L001_R2_001.fastq.gz")
  
  demultiplex_reads(R2, "R2", I1Names, samples, completeMetadata,output.dir)
  
  writeLog('completed demultiplexing R2')
  
  file.create(file.path(output.dir,'demultiplex.done'))

