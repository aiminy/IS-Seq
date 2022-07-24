#!/usr/bin/env Rscript

packages <- c("optparse","stringr","GenomicRanges")

zzz<-lapply(packages, function(xxx) suppressMessages(require(xxx, character.only = TRUE,quietly=TRUE,warn.conflicts = FALSE)))


args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)>=1) {
  # default output file
  input_sam=args[1]  	
  input_repeatMasker=args[2]
  input_MAPQ=args[3]
  output.dir=args[4]
  ISA_run=args[5]
}

# Example:
# Rscript ~/ispipe/R/sam2filterNo.R /home/ubuntu/SHARE/ISseqOutput/Dec282021/CutAdapt/align/R1_R2_Barcode_FB-P5-Rd1-LTR.9_FB-P7-Rd2-LC.9_aligned_mem.sam /home/ubuntu/SHARE/D32_Platform_Development/MANUSCRIPTS_DATA/ISAtest/MiSeqTest/utilsRefData/hg38/repeatMaskerhg38BED 0 /home/ubuntu/SHARE/Aimin/TestSimulation/CL6 POOL-ISA-AVRO-6-Preclin

# grep 49461738 /home/ubuntu/SHARE/Aimin/TestSimulation/CL6/POOL-ISA-AVRO-6-Preclin_FB-P5-Rd1-LTR.9_FB-P7-Rd2-LC.9_final_parse_filterNo.txt | wc -l


# Rscript ~/ispipe/R/sam2filterNo.R /home/ubuntu/SHARE/Aimin/TestSimulation/R1_R2_Barcode_FB-P5-Rd1-LTR.9_FB-P7-Rd2-LC.9_aligned_mem.sam /home/ubuntu/SHARE/D32_Platform_Development/MANUSCRIPTS_DATA/ISAtest/MiSeqTest/utilsRefData/hg38/repeatMaskerhg38BED 0 /home/ubuntu/SHARE/Aimin/TestSimulation POOL-ISA-AVRO-6-Preclin

# grep 49461738  /home/ubuntu/SHARE/Aimin/TestSimulation/POOL-ISA-AVRO-6-Preclin_FB-P5-Rd1-LTR.9_FB-P7-Rd2-LC.9_final_parse_filterNo.txt | wc -l

# Rscript ~/ispipe/R/sam2filterNo.R /home/ubuntu/SHARE/Aimin/TestSimulation/R1_R2_Barcode_FB-P5-Rd1-LTR.9_FB-P7-Rd2-LC.9_aligned_mem.sam /home/ubuntu/SHARE/D32_Platform_Development/MANUSCRIPTS_DATA/ISAtest/MiSeqTest/utilsRefData/hg38/repeatMaskerhg38BED 0 /home/ubuntu/SHARE/Aimin/TestSimu4 POOL-ISA-AVRO-6-Preclin

#Rscript ~/ispipe/R/sam2filterNo.R /home/ubuntu/SHARE/Aimin/TestSimulation/FragR1_R2_Barcode_FB-P5-Rd1-LTR.9_FB-P7-Rd2-LC.9_aligned_mem.sam /home/ubuntu/SHARE/D32_Platform_Development/MANUSCRIPTS_DATA/ISAtest/MiSeqTest/utilsRefData/hg38/repeatMaskerhg38BED 0 /home/ubuntu/SHARE/Aimin/TestSimu5 POOL-ISA-AVRO-6-Preclin

# Rscript ~/ispipe/R/sam2filterNo.R /home/ubuntu/SHARE/Aimin/TestSimulation/FragR1_R2_Barcode_FB-P5-Rd1-LTR.9_FB-P7-Rd2-LC.9_aligned_mem.sam /home/ubuntu/SHARE/D32_Platform_Development/MANUSCRIPTS_DATA/ISAtest/MiSeqTest/utilsRefData/hg38/repeatMaskerhg38BED 0 /home/ubuntu/SHARE/Aimin/TestSimu7 POOL-ISA-AVRO-6-Preclin

# Rscript ~/ispipe/R/sam2filterNo.R /home/ubuntu/SHARE/Aimin/TestSimulation/FragRevR1_R2_Barcode_FB-P5-Rd1-LTR.9_FB-P7-Rd2-LC.9_aligned_mem.sam /home/ubuntu/SHARE/D32_Platform_Development/MANUSCRIPTS_DATA/ISAtest/MiSeqTest/utilsRefData/hg38/repeatMaskerhg38BED 0 /home/ubuntu/SHARE/Aimin/TestSimuRev POOL-ISA-AVRO-6-Preclin

#/home/ubuntu/SHARE/Aimin/TestSimulation/UpFragR1_R2_Barcode_FB-P5-Rd1-LTR.9_FB-P7-Rd2-LC.9_aligned_mem.sam

# Rscript ~/ispipe/R/sam2filterNo.R /home/ubuntu/SHARE/Aimin/TestSimulation/UpFragR1_R2_Barcode_FB-P5-Rd1-LTR.9_FB-P7-Rd2-LC.9_aligned_mem.sam /home/ubuntu/SHARE/D32_Platform_Development/MANUSCRIPTS_DATA/ISAtest/MiSeqTest/utilsRefData/hg38/repeatMaskerhg38BED 0 /home/ubuntu/SHARE/Aimin/TestSimuUpFrag POOL-ISA-AVRO-6-Preclin

# Rscript ~/ispipe/R/sam2filterNo.R /home/ubuntu/SHARE/Aimin/TestSimulation/UpFrag49461738R1_R2_Barcode_FB-P5-Rd1-LTR.9_FB-P7-Rd2-LC.9_aligned_mem.sam /home/ubuntu/SHARE/D32_Platform_Development/MANUSCRIPTS_DATA/ISAtest/MiSeqTest/utilsRefData/hg38/repeatMaskerhg38BED 0 /home/ubuntu/SHARE/Aimin/TestSimuUpFrag Up-Frag-POOL-ISA-AVRO-6-Preclin




print(input_sam)
print(input_repeatMasker)
print(input_MAPQ)
print(output.dir)

if(!dir.exists(output.dir)){dir.create(output.dir,recursive = TRUE)}

output.bam <- file.path(output.dir,'BAM')

if(!dir.exists(output.bam)){dir.create(output.bam,recursive = TRUE)}

#input_sam = '/home/ubuntu/SHARE/ISseqOutput/Dec282021/CutAdapt/align/R1_R2_Barcode_FB-P5-Rd1-LTR.9_FB-P7-Rd2-LC.9_aligned_mem.sam'

sampleParser <- function(input_sam) {
  x <- basename(input_sam)  
  s0 <- unlist(lapply(str_locate_all(x,"aligned"),"[",c(1)))-2
  sampleName <- str_sub(x,1,s0)
  sampleName
}

sampleName <- sampleParser(input_sam)

output.file.0 <- file.path(output.bam,paste0(sampleName,'_aligned_mem.bam'))
  
t.input <- file.info(input_sam)$mtime
t.output.0 <- file.info(output.file.0)$mtime

if(is.na(t.output.0)|(t.input>t.output.0)){
cmd= paste0('samtools view -bS ',input_sam,' > ',output.file.0)  
cat(cmd,'\n')
system(cmd)
}

# samtools view -bS /home/ubuntu/SHARE/ISseqOutput/Dec282021/CutAdapt/align/R1_R2_Barcode_FB-P5-Rd1-LTR.9_FB-P7-Rd2-LC.9_aligned_mem.sam > /home/ubuntu/SHARE/ISseqOutput/Dec282021/CutAdapt/BAM/R1_R2_Barcode_FB-P5-Rd1-LTR.9_FB-P7-Rd2-LC.9_aligned_mem.bam
# 
# 

output.file.1<- file.path(output.bam,paste0(sampleName,'_aligned_mem_mapped.bam'))

t.output.1 <- file.info(output.file.1)$mtime

if(is.na(t.output.1)|(t.output.0>t.output.1)){
  cmd= paste0('samtools view -b -F 4 ',output.file.0,' > ',output.file.1)  
  cat(cmd,'\n')
  system(cmd)
}

# samtools view -b -F 4 /home/ubuntu/SHARE/ISseqOutput/Dec282021/CutAdapt/BAM/R1_R2_Barcode_FB-P5-Rd1-LTR.9_FB-P7-Rd2-LC.9_aligned_mem.bam > /home/ubuntu/SHARE/ISseqOutput/Dec282021/CutAdapt/BAM/R1_R2_Barcode_FB-P5-Rd1-LTR.9_FB-P7-Rd2-LC.9_aligned_mem_mapped.bam

output.file.2<- file.path(output.bam,paste0(sampleName,'_aligned_mem_mapped_primary.bam'))

t.output.2 <- file.info(output.file.2)$mtime

if(is.na(t.output.2)|(t.output.1>t.output.2)){
  cmd= paste0('samtools view -b -F 256 ',output.file.1,' > ',output.file.2)  
  cat(cmd,'\n')
  system(cmd)
}


# samtools view -b -F 256 /home/ubuntu/SHARE/ISseqOutput/Dec282021/CutAdapt/BAM/R1_R2_Barcode_FB-P5-Rd1-LTR.9_FB-P7-Rd2-LC.9_aligned_mem_mapped.bam > /home/ubuntu/SHARE/ISseqOutput/Dec282021/CutAdapt/BAM/R1_R2_Barcode_FB-P5-Rd1-LTR.9_FB-P7-Rd2-LC.9_aligned_mem_mapped_primary.bam
 

output.bam.1 <- file.path(output.dir,'BAMSorted')

if(!dir.exists(output.bam.1)){dir.create(output.bam.1,recursive = TRUE)}

output.file.3<- file.path(output.bam.1,paste0(sampleName,'_aligned_mem_mapped_primary_sort.bam'))

t.output.3 <- file.info(output.file.3)$mtime

if(is.na(t.output.3)|(t.output.2>t.output.3)){
  cmd= paste0('samtools sort ',output.file.2,' -o ',output.file.3)  
  cat(cmd,'\n')
  system(cmd)
}

# samtools sort /home/ubuntu/SHARE/ISseqOutput/Dec282021/CutAdapt/BAM/R1_R2_Barcode_FB-P5-Rd1-LTR.9_FB-P7-Rd2-LC.9_aligned_mem_mapped_primary.bam -o /home/ubuntu/SHARE/ISseqOutput/Dec282021/CutAdapt/BAMSorted/R1_R2_Barcode_FB-P5-Rd1-LTR.9_FB-P7-Rd2-LC.9_aligned_mem_mapped_primary_sort.bam

output.file.4<- file.path(output.bam.1,paste0(sampleName,'_aligned_mem_mapped_primary_sort.bam.index'))

t.output.4 <- file.info(output.file.4)$mtime

if(is.na(t.output.4)|(t.output.3>t.output.4)){
  cmd= paste0('samtools index ',output.file.3)  
  cat(cmd,'\n')
  system(cmd)
}


# 
# samtools index /home/ubuntu/SHARE/ISseqOutput/Dec282021/CutAdapt/BAMSorted/R1_R2_Barcode_FB-P5-Rd1-LTR.9_FB-P7-Rd2-LC.9_aligned_mem_mapped_primary_sort.bam
# 

output.file.5 <- file.path(output.bam.1,paste0(sampleName,'_aligned_mem_sort_inMask.bam'))

t.output.5 <- file.info(output.file.5)$mtime

if(is.na(t.output.5)|(t.output.3>t.output.5)){
  cmd= paste0('bedtools intersect -abam ',output.file.3,' -b ',input_repeatMasker,' > ',output.file.5)  
  cat(cmd,'\n')
  system(cmd)
}


# bedtools intersect -abam /home/ubuntu/SHARE/ISseqOutput/Dec282021/CutAdapt/BAMSorted/R1_R2_Barcode_FB-P5-Rd1-LTR.9_FB-P7-Rd2-LC.9_aligned_mem_mapped_primary_sort.bam -b /home/ubuntu/SHARE/D32_Platform_Development/test/ISAtest/MiSeqTest/utilsRefData/hg38/repeatMaskerhg38BED > /home/ubuntu/SHARE/ISseqOutput/Dec282021/CutAdapt/BAMSorted/R1_R2_Barcode_FB-P5-Rd1-LTR.9_FB-P7-Rd2-LC.9_aligned_mem_sort_inMask.bam

output.file.6 <- file.path(output.bam.1,paste0(sampleName,'_aligned_mem_sort_nonMask.bam'))

t.output.6 <- file.info(output.file.6)$mtime

if(is.na(t.output.6)|(t.output.3>t.output.6)){
  cmd= paste0('bedtools intersect -v -abam ',output.file.3,' -b ',input_repeatMasker,' > ',output.file.6)  
  cat(cmd,'\n')
  system(cmd)
}

# bedtools intersect -v -abam /home/ubuntu/SHARE/ISseqOutput/Dec282021/CutAdapt/BAMSorted/R1_R2_Barcode_FB-P5-Rd1-LTR.9_FB-P7-Rd2-LC.9_aligned_mem_mapped_primary_sort.bam -b /home/ubuntu/SHARE/D32_Platform_Development/test/ISAtest/MiSeqTest/utilsRefData/hg38/repeatMaskerhg38BED > /home/ubuntu/SHARE/ISseqOutput/Dec282021/CutAdapt/BAMSorted/R1_R2_Barcode_FB-P5-Rd1-LTR.9_FB-P7-Rd2-LC.9_aligned_mem_sort_nonMask.bam

 
output.file.7 <- file.path(output.bam.1,paste0(sampleName,'_aligned_mem_sort_inMask_qual.bam'))

t.output.7 <- file.info(output.file.7)$mtime

if(is.na(t.output.7)|(t.output.3>t.output.7)){
  cmd= paste0('samtools view -bq ',input_MAPQ,' ',output.file.5,' > ',output.file.7)  
  cat(cmd,'\n')
  system(cmd)
}




# samtools view -bq 0 /home/ubuntu/SHARE/ISseqOutput/Dec282021/CutAdapt/BAMSorted/R1_R2_Barcode_FB-P5-Rd1-LTR.9_FB-P7-Rd2-LC.9_aligned_mem_sort_inMask.bam > /home/ubuntu/SHARE/ISseqOutput/Dec282021/CutAdapt/BAMSorted/R1_R2_Barcode_FB-P5-Rd1-LTR.9_FB-P7-Rd2-LC.9_aligned_mem_sort_inMask_qual.bam


output.file.8 <- file.path(output.bam.1,paste0(sampleName,'_aligned_mem_allFilter.bam'))

t.output.8 <- file.info(output.file.8)$mtime

if(is.na(t.output.8)|(t.output.7>t.output.8)){
  cmd= paste0('samtools merge ',output.file.8,' ',output.file.6,' ',output.file.7)  
  cat(cmd,'\n')
  system(cmd)
}

# samtools merge /home/ubuntu/SHARE/ISseqOutput/Dec282021/CutAdapt/BAMSorted/R1_R2_Barcode_FB-P5-Rd1-LTR.9_FB-P7-Rd2-LC.9_aligned_mem_allFilter.bam /home/ubuntu/SHARE/ISseqOutput/Dec282021/CutAdapt/BAMSorted/R1_R2_Barcode_FB-P5-Rd1-LTR.9_FB-P7-Rd2-LC.9_aligned_mem_sort_nonMask.bam /home/ubuntu/SHARE/ISseqOutput/Dec282021/CutAdapt/BAMSorted/R1_R2_Barcode_FB-P5-Rd1-LTR.9_FB-P7-Rd2-LC.9_aligned_mem_sort_inMask_qual.bam
# 

output.file.9 <- file.path(output.bam.1,paste0(sampleName,'_aligned_mem_allFilter_rehead.bam'))

t.output.9 <- file.info(output.file.9)$mtime

if(is.na(t.output.9)|(t.output.8>t.output.9)){
  cmd= paste0('PicardCommandLine AddOrReplaceReadGroups I=',output.file.8,' O=',output.file.9,' SORT_ORDER=coordinate RGID=foo RGLB=bar RGPL=illumina RGPU=shoot RGSM=DePristo')
  cat(cmd,'\n')
  system(cmd)
}

# PicardCommandLine AddOrReplaceReadGroups I=/home/ubuntu/SHARE/ISseqOutput/Dec282021/CutAdapt/BAMSorted/R1_R2_Barcode_FB-P5-Rd1-LTR.9_FB-P7-Rd2-LC.9_aligned_mem_allFilter.bam O=/home/ubuntu/SHARE/ISseqOutput/Dec282021/CutAdapt/BAMSorted/R1_R2_Barcode_FB-P5-Rd1-LTR.9_FB-P7-Rd2-LC.9_aligned_mem_allFilter_rehead.bam SORT_ORDER=coordinate RGID=foo RGLB=bar RGPL=illumina RGPU=shoot RGSM=DePristo

output.file.10 <- file.path(output.bam.1,paste0(sampleName,'_aligned_mem_allFilter_rehead.bam.index'))

t.output.10 <- file.info(output.file.10)$mtime

if(is.na(t.output.10)|(t.output.9>t.output.10)){
  cmd= paste0('samtools index ',output.file.9)
  cat(cmd,'\n')
  system(cmd)
}
 
# samtools index /home/ubuntu/SHARE/ISseqOutput/Dec282021/CutAdapt/BAMSorted/R1_R2_Barcode_FB-P5-Rd1-LTR.9_FB-P7-Rd2-LC.9_aligned_mem_allFilter_rehead.bam

output.file.11 <- file.path(output.bam.1,paste0(sampleName,'_aligned_mem_allFilter_rehead_exact3nt.bam'))

t.output.11 <- file.info(output.file.11)$mtime

if(is.na(t.output.11)|(t.output.9>t.output.11)){
  cmd= paste0('export PYTHONPATH=/home/ubuntu/miniconda2/lib/python2.7/site-packages;python2 /home/ubuntu/ispipe/utils/pysam_parse.py ',output.file.9)
  cat(cmd,'\n')
  system(cmd)
}

# python2 /home/ubuntu/ispipe/utils/pysam_parse.py /home/ubuntu/SHARE/ISseqOutput/Dec282021/CutAdapt/BAMSorted/R1_R2_Barcode_FB-P5-Rd1-LTR.9_FB-P7-Rd2-LC.9_aligned_mem_allFilter_rehead.bam

output.file.12 <- file.path(output.bam.1,paste0(sampleName,'_aligned_mem_allFilter_rehead_exact3nt_unmapped.bam'))

t.output.12 <- file.info(output.file.12)$mtime

if(is.na(t.output.12)|(t.output.11>t.output.12)){
  cmd= paste0('samtools view -b -f 4 ',output.file.11,' > ',output.file.12)
  cat(cmd,'\n')
  system(cmd)
}

# samtools view -b -f 4 /home/ubuntu/SHARE/ISseqOutput/Dec282021/CutAdapt/BAMSorted/R1_R2_Barcode_FB-P5-Rd1-LTR.9_FB-P7-Rd2-LC.9_aligned_mem_allFilter_rehead_exact3nt.bam > /home/ubuntu/SHARE/ISseqOutput/Dec282021/CutAdapt/BAMSorted/R1_R2_Barcode_FB-P5-Rd1-LTR.9_FB-P7-Rd2-LC.9_aligned_mem_allFilter_rehead_exact3nt_unmapped.bam
# 

output.file.13 <- file.path(output.bam.1,paste0(sampleName,'_aligned_mem_allFilter_rehead_exact3nt_supplementary.bam'))

t.output.13 <- file.info(output.file.13)$mtime

if(is.na(t.output.13)|(t.output.11>t.output.13)){
  cmd= paste0('samtools view -b -f 2048 ',output.file.11,' > ',output.file.12)
  cat(cmd,'\n')
  system(cmd)
}


# samtools view -b -f 2048 /home/ubuntu/SHARE/ISseqOutput/Dec282021/CutAdapt/BAMSorted/R1_R2_Barcode_FB-P5-Rd1-LTR.9_FB-P7-Rd2-LC.9_aligned_mem_allFilter_rehead_exact3nt.bam > /home/ubuntu/SHARE/ISseqOutput/Dec282021/CutAdapt/BAMSorted/R1_R2_Barcode_FB-P5-Rd1-LTR.9_FB-P7-Rd2-LC.9_aligned_mem_allFilter_rehead_exact3nt_supplementary.bam
# 

output.file.14 <- file.path(output.bam.1,paste0(sampleName,'_aligned_mem_allFilter_rehead_exact3nt_nonSupplementary.bam'))

t.output.14 <- file.info(output.file.14)$mtime

if(is.na(t.output.14)|(t.output.11>t.output.14)){
  cmd= paste0('samtools view -b -F 2048 ',output.file.11,' > ',output.file.14)
  cat(cmd,'\n')
  system(cmd)
}




# samtools view -b -F 2048 /home/ubuntu/SHARE/ISseqOutput/Dec282021/CutAdapt/BAMSorted/R1_R2_Barcode_FB-P5-Rd1-LTR.9_FB-P7-Rd2-LC.9_aligned_mem_allFilter_rehead_exact3nt.bam > /home/ubuntu/SHARE/ISseqOutput/Dec282021/CutAdapt/BAMSorted/R1_R2_Barcode_FB-P5-Rd1-LTR.9_FB-P7-Rd2-LC.9_aligned_mem_allFilter_rehead_exact3nt_nonSupplementary.bam


output.file.15 <- file.path(output.bam.1,paste0(sampleName,'_aligned_mem_allFilter_rehead_exact3nt_nonSupplementary.bam.index'))

t.output.15 <- file.info(output.file.15)$mtime

if(is.na(t.output.15)|(t.output.14>t.output.15)){
  cmd= paste0('samtools index ',output.file.14)
  cat(cmd,'\n')
  system(cmd)
}

# samtools index /home/ubuntu/SHARE/ISseqOutput/Dec282021/CutAdapt/BAMSorted/R1_R2_Barcode_FB-P5-Rd1-LTR.9_FB-P7-Rd2-LC.9_aligned_mem_allFilter_rehead_exact3nt_nonSupplementary.bam
# R1_R2_Barcode_FB-P5-Rd1-LTR.9_FB-P7-Rd2-LC.9_aligned_mem_allFilter_rehead_exact3nt_nonSupplementary.bam

#POOL-ISA-AVRO-6-Preclin_FB-P5-Rd1-LTR.9_FB-P7-Rd2-LC.9_final_parse_filterNo.txt

output.file.16 <- file.path(output.dir,paste0(ISA_run,'_FB-P5-Rd1-LTR.9_FB-P7-Rd2-LC.9_final_parse_filterNo.txt'))

t.output.16 <- file.info(output.file.16)$mtime

if(is.na(t.output.16)|(t.output.14>t.output.16)){
  cmd= paste0('export PYTHONPATH=/home/ubuntu/miniconda2/lib/python2.7/site-packages;python /home/ubuntu/ispipe/utils/try_pysam.py ',output.file.14,' ',ISA_run)
  cat(cmd,'\n')
  system(cmd)
}
 
# python /home/ubuntu/ispipe/utils/try_pysam.py /home/ubuntu/SHARE/ISseqOutput/Dec282021/CutAdapt/BAMSorted/R1_R2_Barcode_FB-P5-Rd1-LTR.9_FB-P7-Rd2-LC.9_aligned_mem_allFilter_rehead_exact3nt_nonSupplementary.bam POOL-ISA-AVRO-6-Preclin
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
