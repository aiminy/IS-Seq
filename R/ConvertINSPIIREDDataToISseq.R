# Example:

# check which python is used

# echo $(python -c "import site; print(site.getsitepackages()[0])")/home/ubuntu/miniconda2/lib/python2.7/site-packages

# PYTHONPATH='/home/ubuntu/miniconda2/lib/python2.7/site-packages'

# Rscript /home/ubuntu/intsitecaller/ConvertINSPIIREDDataToISseq.R /home/ubuntu/intsitecaller/testCases/intSiteValidation/Data/Undetermined_S0_L001_I1_001.fastq.gz /home/ubuntu/intsitecaller/testCases/intSiteValidation/Data/Undetermined_S0_L001_R1_001.fastq.gz /home/ubuntu/intsitecaller/testCases/intSiteValidation/Data/Undetermined_S0_L001_R2_001.fastq.gz /home/ubuntu/SHARE/ISseqOutput/INSPIIREDData/test R2.fq R1.fq /home/ubuntu/miniconda2/lib/python2.7/site-packages

initial.options <- commandArgs(trailingOnly = FALSE)
print(initial.options)

file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.dirname <- dirname(script.name)
print(script.dirname)

library("ShortRead")

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)>=1) {

  # input
  input.index=args[1]
  input.R1=args[2]
  input.R2=args[3]

  # output
  output.dir=args[4]
  output.R1=args[5]
  output.R2=args[6]
  PYTHONPATH=args[7]

}

if(!dir.exists(output.dir)){dir.create(output.dir,recursive = TRUE)}

output.R1= file.path(output.dir,'R2.fq')
output.R2=file.path(output.dir,'R1.fq')

I1 <- readFastq(input.index)

I1 <- trimTailw(I1, 2, "0", 12)
I1 <- I1[width(I1)==max(width(I1))]
I1 <- split(I1, ceiling(seq_along(I1)/500000))

PYTHONPATH=PYTHONPATH

for(chunk in names(I1)){

  output_fasta=file.path(output.dir,paste0("trimmedI1-", chunk, ".fasta"))
  writeFasta(I1[[chunk]], file=output_fasta)

  cmd <- paste0('export PYTHONPATH=',PYTHONPATH,';','python ',file.path(script.dirname,'processGolayTest.py '),output_fasta,' ',file.path(output.dir,'correctedI1-1.fasta'),' ',file.path(output.dir,'correctedI1-1.done'))

  print(cmd)

  system(cmd)

}

GetIndex <- function(Index.input.fq,Index.input.fa) {

  I1.fa <- readFasta(Index.input.fa)

  I1.fq <- readFastq(Index.input.fq)
  I1.fq <- trimTailw(I1.fq, 2, "0", 12)
  I1.fq <- I1.fq[width(I1.fq)==max(width(I1.fq))]

  I1.fq.filter <- I1.fq[match(id(I1.fa),sapply(strsplit(as.character(ShortRead::id(I1.fq)), " "), "[[", 1))]

  I1.fq.filter@sread <- I1.fa@sread

  I1.fq.filter@id <- I1.fa@id


  I1.fq.matched <- I1.fq.filter
}



Index.input.fa <- file.path(output.dir,'correctedI1-1.fasta')

Index.input.fq <- input.index

I1.fq.matched <- GetIndex(Index.input.fq,Index.input.fa)

AppendIndexToRead <- function(I1.fq.matched,Input.read,Output.read) {

  R1 <- readFastq(Input.read)

  R1.fq.matched <- R1[match(id(I1.fq.matched), sapply(strsplit(as.character(ShortRead::id(R1)), " "), "[[", 1))]

  Index_fq_matched <- tempfile()
  Read_fq_matched <- tempfile()

  writeFastq(I1.fq.matched, file=Index_fq_matched, mode="w", full=FALSE, compress=FALSE)
  writeFastq(R1.fq.matched, file=Read_fq_matched, mode="w", full=FALSE, compress=FALSE)

  cmd <- paste('sh ',file.path(script.dirname,'PasteIndexToRead.sh '),Index_fq_matched,Read_fq_matched,Output.read)
  print(cmd)
  system(cmd)
}

Input.read <- input.R1
Output.read <- output.R1
AppendIndexToRead(I1.fq.matched,Input.read,Output.read)

Input.read <- input.R2
Output.read <- output.R2
AppendIndexToRead(I1.fq.matched,Input.read,Output.read)
