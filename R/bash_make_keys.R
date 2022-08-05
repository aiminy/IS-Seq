#!/usr/bin/env Rscript

initial.options <- commandArgs(trailingOnly = FALSE)
print(initial.options)

file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.dirname <- dirname(script.name)
print(script.dirname)

# Rscript /home/ubuntu/DEMO/IS-Seq/R/bash_make_keys.R /home/ubuntu/SHARE/D32_Platform_Development/MANUSCRIPTS_DATA/ISAtest/MiSeqTest/20210121_AssociationFIle_POOL6_Preclinical.csv /home/ubuntu/SHARE/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa /home/ubuntu/SHARE/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)>=1) {
  # default output file
  input.file=args[1]
  input.dir=args[2]
  output.dir=args[3]
}

print(input.file)
print(input.dir)
print(output.dir)

if(!dir.exists(output.dir)){dir.create(output.dir,recursive = TRUE)}

input <- read.table(input.file,sep=',',header=T)

#print(input)


#print(input)

SampleName=gsub('-','',input$Sample.Type)

input.R1 = file.path(input.dir,SampleName,paste0('R1_fastq_trim12nt_qcTrimmed_MatchBlastLtrLc_Barcode_',input$Fusion.Primer.LTR..ID,'.fq_trimwithCutAdapt','_',SampleName,'_ReadyToAlignSort'))

#print(input.R1)

input.R2 = file.path(input.dir,SampleName,paste0('R2_fastq_trim12nt_qcTrimmed_MatchBlastLtrLc_Barcode_',input$Fusion.Primer.LC..ID,'.fq_trimwithCutAdapt','_',SampleName,'_ReadyToAlignSort'))

#print(input.R2)


input.all <- data.frame(R1=input.R1,R2=input.R2,SampleName=gsub('-','',input$Sample.Type),sampleDir=file.path(output.dir,SampleName))

#print(input.all)


if(TRUE){

cat('Start conversion\n')

#input.all <- data.frame(R1=input.R1,R2=input.R2,SampleName=gsub('-','',input$Sample.Type))


#print(input.all)

#input.all <- input.all[-c(4,9,10),]

null <- lapply(1:dim(input.all)[1], function(u){
 
  
  sample.dir <- input.all[u,]$sampleDir

  #print(sample.dir)

  if(!dir.exists(sample.dir)){dir.create(sample.dir,recursive = TRUE)}

  cmd.R1 <- paste0('Rscript ',file.path(script.dirname,'MakeKeys.R '),input.all[u,]$R1,' ',input.all[u,]$R2,' ',input.all[u,]$SampleName,' ', sample.dir)

  #print(cmd.R1)
  
  #system(cmd.R1)

  cmd.R2 <- paste0('Rscript ',file.path(script.dirname,'PslToIs_one_replicate_change_sequence_similarity.R '),paste0(input.all[u,]$R2,'.fa.psl'),' ',paste0(input.all[u,]$R1,'.fa.psl'),' ',file.path(sample.dir,'keys.rds '),'/home/ubuntu/DEMO/IS-Seq/utilsRefData/INSPIIRED/completeMetadata.RData ',file.path(sample.dir,'IS803 '),'hg38 ',1,' ',80)


  print(cmd.R2)
  system(cmd.R2)

})
}

