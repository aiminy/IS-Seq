#!/usr/bin/env Rscript

#

# Rscript ~/intsitecaller/MakeKeys.R /home/ubuntu/SHARE/ISseqOutput/Oct2721/IsaByINSPIIRED/fa/R1_fastq_trim12nt_qcTrimmed_MatchBlastLtrLc_Barcode_FB-P5-Rd1-LTR.13.fq_trimwithCutAdapt_HL60_Cl6_ReadyToAlignSort /home/ubuntu/SHARE/ISseqOutput/Oct2721/IsaByINSPIIRED/fa/R2_fastq_trim12nt_qcTrimmed_MatchBlastLtrLc_Barcode_FB-P7-Rd2-LC.20.fq_trimwithCutAdapt_HL60_Cl6_ReadyToAlignSort HL60Cl6 /home/ubuntu/SHARE/ISseqOutput/Oct2721/IsaByINSPIIRED


args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)>=1) {
  # default output file
  input_R1=args[1]
  input_R2=args[2]
  sampleName=args[3]
  output.dir=args[4]
}

#output.dir <- file.path(output.dir,sampleName)

print(input_R1)
print(input_R2)
print(output.dir)


if(!dir.exists(output.dir)){dir.create(output.dir,recursive = TRUE)}


library(data.table)
#aws_root <- '/home/ubuntu/SHARE/ISseqOutput/Oct2721/IsaByINSPIIRED'

#R1.readName <- data.table::fread(paste0("cat < /home/ubuntu/ISseqOutput/Oct2721/R1_fastq_trim12nt_qcTrimmed_MatchBlastLtrLcDEMULTIPLEXINGTofq/R1_fastq_trim12nt_qcTrimmed_MatchBlastLtrLc_Barcode_FB-P5-Rd1-LTR.13.fq_trimwithCutAdapt | grep '@M'"),header=F)

#R2.readName <- data.table::fread(paste0("cat < /home/ubuntu/ISseqOutput/Oct2721/RandomBarcodRemovalOutPut4R2CuReRun/R2_fastq_trim12nt_qcTrimmed_MatchBlastLtrLc_Barcode_FB-P7-Rd2-LC.15.fq_trimwithCutAdapt | grep '@M'"),header=F)

cmd1 <- paste0("cat < ",input_R1," | grep '@M'")

cmd2 <- paste0("cat < ",input_R2," | grep '@M'")


R1.readName <- data.table::fread(cmd1,header=F)

R2.readName <- data.table::fread(cmd2,header=F)


R1Name <- sapply(strsplit(as.character(R1.readName$V1), "@"), "[[", 2)

R1 <- data.frame(R1=seq(1:length(R1Name)),names=R1Name)

R2Name <- sapply(strsplit(as.character(R2.readName$V1), "@"), "[[", 2)

R2 <- data.frame(R2=seq(1:length(R2Name)),names=R2Name)

cat('R1','\n')
print(head(R1))

cat('R2','\n')
print(head(R2))

# 


#keys <- intersect(R1Name,R2Name)

keys <- merge(R2,R1,by='names')

#cat('R1R2','\n')
#print(head(keys))

keys$readPairKey <- paste0(keys$R2,'_',keys$R1)

keys <- keys[,c(2,3,1,4)]

keys$names <- paste0(sampleName,'%',keys$names)

cat('R1R2','\n')
print(head(keys))

saveRDS(keys,file.path(output.dir,'keys.rds'))


# 
# /home/ubuntu/ISseqOutput/Oct2721/R1_fastq_trim12nt_qcTrimmed_MatchBlastLtrLcDEMULTIPLEXINGTofq/R1_fastq_trim12nt_qcTrimmed_MatchBlastLtrLc_Barcode_FB-P5-Rd1-LTR.13.fq_trimwithCutAdapt_ReadyToAlignSort 
# 
# /home/ubuntu/ISseqOutput/Oct2721/RandomBarcodRemovalOutPut4R2CuReRun/R2_fastq_trim12nt_qcTrimmed_MatchBlastLtrLc_Barcode_FB-P7-Rd2-LC.15.fq_trimwithCutAdapt_ReadyToAlignSort
#  
