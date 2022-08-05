df.5.IS.random <- read.table('/Users/aiminyan/Aimin/aws_root/home/ubuntu/DEMO/ISseqOutput/DEMOPos5ISRandom/CutAdapt/filterNo/db/DEMOPos5ISRandom/FinalOut_DEMOPos5ISRandom/POOL-ISA-AVRO-6-Preclin_HL60POS-CTRL-1CL-6_HL60_DEMOPos5ISRandom_CollisionClean_CollisionTable.txt',header = T)

df.5.IS.random$Relative <- df.5.IS.random$CL.6_POS.CTRL.1_1/sum(df.5.IS.random$CL.6_POS.CTRL.1_1)

df.5.IS.random.1 <- df.5.IS.random[order(df.5.IS.random$V7),]

df.5.IS.random.1$IS <- paste0(df.5.IS.random.1$chr,'_',df.5.IS.random.1$pos,'_',df.5.IS.random.1$strand)

write.table(df.5.IS.random.1[,c(8,5,6,7)],file ='/Volumes/avrobio/Shared/D32_PlatformDevelopment_BIASCO/MANUSCRIPTS_DATA/ISSeq_paper/Simulation4.txt',quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = FALSE,col.names = FALSE)