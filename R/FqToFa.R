#!/usr/bin/env Rscript


#Rscript /home/ubuntu/IS-Seq/R/FqToFa.R /home/ubuntu/DEMO/ISseqOutput/DEMO/R1_fastq_trim12nt_qcTrimmed_MatchBlastLtrLcDEMULTIPLEXINGTofq/R1_fastq_trim12nt_qcTrimmed_MatchBlastLtrLc_Barcode_FB-P5-Rd1-LTR.9.fq_trimwithCutAdapt /home/ubuntu/DEMO/ISseqOutput/DEMO/CutAdapt/R1_R2_Barcode_FB-P5-Rd1-LTR.9_FB-P7-Rd2-LC.9_trimmedID demo /home/ubuntu/DEMO/ISseqOutput/DEMO/INSPIIRED /home/ubuntu/DEMO/IS-Seq/utilsRefData/IsSeq/hg38/hg38ChrOnly.fa 

#Rscript /home/ubuntu/IS-Seq/R/FqToFa.R /home/ubuntu/DEMO/ISseqOutput/DEMO/R2_fastq_trim12nt_qcTrimmed_MatchBlastLtrLcDEMULTIPLEXINGTofq/R2_fastq_trim12nt_qcTrimmed_MatchBlastLtrLc_Barcode_FB-P7-Rd2-LC.9.fq_trimwithCutAdapt /home/ubuntu/DEMO/ISseqOutput/DEMO/CutAdapt/R1_R2_Barcode_FB-P5-Rd1-LTR.9_FB-P7-Rd2-LC.9_trimmedID demo /home/ubuntu/DEMO/ISseqOutput/DEMO/INSPIIRED /home/ubuntu/DEMO/IS-Seq/utilsRefData/IsSeq/hg38/hg38ChrOnly.fa

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)>=1) {
  # default output file
  input_LTR=args[1]  	
  input_LTR_LC=args[2]
  sampleName=args[3]
  output.dir=args[4]
  ref.genome=args[5]
}

print(input_LTR)
print(input_LTR_LC)
print(output.dir)

if(!dir.exists(output.dir)){dir.create(output.dir,recursive = TRUE)}

temp <- paste0(file.path(output.dir,basename(input_LTR)),'_',sampleName,'_ReadyToAlign')

cmd <- paste0('seqtk subseq ',input_LTR,' ', input_LTR_LC,' > ',temp)

t.input <- file.info(input_LTR)$mtime
t.output <- file.info(temp)$mtime

if(is.na(t.output)|(t.input>t.output)){
print(cmd)
system(cmd)
}

temp1 <- paste0(temp,'Sort')

cmd <- paste0('export PYTHONPATH="/home/ubuntu/miniconda2/lib/python2.7/site-packages";fastqutils sort ',temp,' > ',temp1)

t.input <- file.info(temp)$mtime
t.output <- file.info(temp1)$mtime

if(is.na(t.output)|(t.input>t.output)){
print(cmd)
system("which python")
#library(reticulate)
#use_condaenv("base") 
system(cmd)
}

temp2 <- paste0(temp1,'.fa')

cmd <- paste0('fastq_to_fasta -r -i ',temp1,' -o ',temp2)

t.input <- file.info(temp1)$mtime
t.output <- file.info(temp2)$mtime

if(is.na(t.output)|(t.input>t.output)){
print(cmd)
system("which python")
system(cmd)
}

temp3 <- paste0(temp2,'.psl')

#ref_g='/home/ubuntu/SHARE/D32_Platform_Development/test/ISAtest/MiSeqTest/utilsRefData/hg38/hg38ChrOnly.fa'
ref_g=ref.genome

cmd <- paste0('blat ',ref_g,' ',temp2,' ',temp3,' ','-tileSize=11 -stepSize=9 -minIdentity=85 -maxIntron=5 -minScore=27 -dots=1000 -out=psl -noHead')

t.input <- file.info(temp2)$mtime
t.output <- file.info(temp3)$mtime

if(is.na(t.output)|(t.input>t.output)){
print(cmd)
system("which python")
system(cmd)
}

