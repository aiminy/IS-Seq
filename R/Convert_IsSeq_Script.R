#!/usr/bin/env Rscript

initial.options <- commandArgs(trailingOnly = FALSE)
print(initial.options)

file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.dirname <- dirname(script.name)
print(script.dirname)

# Rscript /home/ubuntu/intsitecaller/Convert_IsSeq_Script.R /home/ubuntu/SHARE/D32_Platform_Development/test/ISAtest/MiSeqTest/20210121_AssociationFIle_POOL6_Preclinical.csv /home/ubuntu/SHARE/ISseqOutput/Dec282021 /home/ubuntu/SHARE/ISseqOutput/Dec282021/IsaByINSPIIRED/fa

# Rscript /home/ubuntu/intsitecaller/Convert_IsSeq_Script.R /home/ubuntu/SHARE/D32_Platform_Development/test/ISAtest/MiSeqTest/Association_pool_ISA_AVRO_TEST1_add_hg38.csv /home/ubuntu/ISseqOutput/Apr13TestHg38 /home/ubuntu/SHARE/ISseqOutput/Apr13TestHg38/IsaByINSPIIRED/fa

# Rscript /home/ubuntu/DEMO/IS-Seq/R/Convert_IsSeq_Script.R /home/ubuntu/SHARE/D32_Platform_Development/MANUSCRIPTS_DATA/ISAtest/MiSeqTest/20210121_AssociationFIle_POOL6_Preclinical.csv /home/ubuntu/SHARE/ISseqOutput/Dec282021 /home/ubuntu/SHARE/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa /home/ubuntu/DEMO/IS-Seq/utilsRefData/IsSeq/hg38/hg38ChrOnly.fa

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)>=1) {
  # default output file
  input.file=args[1]
  input.dir=args[2]
  output.dir=args[3]
  ref.genome=args[4]
}

print(input.file)
print(input.dir)
print(output.dir)
print(ref.genome)

if(!dir.exists(output.dir)){dir.create(output.dir,recursive = TRUE)}

input <- read.table(input.file,sep=',',header=T)

#print(input)

input.R1 = file.path(input.dir,'R1_fastq_trim12nt_qcTrimmed_MatchBlastLtrLcDEMULTIPLEXINGTofq',paste0('R1_fastq_trim12nt_qcTrimmed_MatchBlastLtrLc_Barcode_',input$Fusion.Primer.LTR..ID,'.fq_trimwithCutAdapt'))

#print(input.R1)

input.R1R2 = file.path(input.dir, paste0('CutAdapt/R1_R2_Barcode_',input$Fusion.Primer.LTR..ID,'_',input$Fusion.Primer.LC..ID,'_trimmedID'))
#print(input.R1R2)

input.R2 = file.path(input.dir,'RandomBarcodRemovalOutPut4R2CuReRun',paste0('R2_fastq_trim12nt_qcTrimmed_MatchBlastLtrLc_Barcode_',input$Fusion.Primer.LC..ID,'.fq_trimwithCutAdapt'))

#print(input.R2)

if(TRUE){

cat('Start conversion\n')

input.all <- data.frame(R1=input.R1,R1R2=input.R1R2,R2=input.R2,SampleName=gsub('-','',input$Sample.Type))


#print(input.all)

#input.all <- input.all[-c(4,9,10),]




null <- lapply(1:dim(input.all)[1], function(u){
 
  
  sample.dir <- file.path(output.dir,input.all[u,]$SampleName)

  if(!dir.exists(sample.dir)){dir.create(sample.dir,recursive = TRUE)}

  cmd.R1 <- paste0('Rscript ',file.path(script.dirname,'FqToFa.R '),input.all[u,]$R1,' ',input.all[u,]$R1R2,' ',input.all[u,]$SampleName,' ', sample.dir,' ',ref.genome)

  print(cmd.R1)
  system(cmd.R1)
  
  cmd.R2 <- paste0('Rscript ',file.path(script.dirname,'FqToFa.R '),input.all[u,]$R2,' ',input.all[u,]$R1R2,' ',input.all[u,]$SampleName,' ', sample.dir,' ',ref.genome)

  print(cmd.R2)
  system(cmd.R2)
  
})


}

