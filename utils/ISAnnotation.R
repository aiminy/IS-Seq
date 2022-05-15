
args = commandArgs(trailingOnly=TRUE)

if (length(args)<1) {
  stop("Two argument first working directory and suffix name for file(date)", call.=FALSE)
} else{
  # default output file
  wd=args[1]
}
library("reshape")

files=list.files(args[1],pattern="*_CollisionClean$",all.files = TRUE,full.names = TRUE,recursive = TRUE,include.dirs = TRUE) 

output.dir <- file.path(args[1],paste0("FinalOut_",args[2]))

cat(" please find the final output in :",output.dir,"\n")

if(!dir.exists(output.dir)){dir.create(output.dir,recursive = TRUE)}

for(f in 1:length(files)){
  so=read.table(files[f],sep = "\t",header=T)
  an=read.table(paste(files[f],"_BedFormat_closest_knownGenesV19_VectMask_uniq.txt",sep=""),sep = "\t",header = F)
  so$label=apply(so,1,function(x) {paste(as.character(x[1]),as.numeric(x[2]),sep='_',collapse='_')})
  an$label=apply(an,1,function(x) {paste(as.character(x[1]),as.numeric(x[2]),sep='_',collapse='_')})
  unik <- !duplicated(an$label)
  anU=an[unik ,c(10,7)]
  allsoUnique=merge( anU,so, by.x ='label', by.y = 'label', all.x = T, all.y = F)
  tl=dim(allsoUnique)[2]
  allsoUniqueP=allsoUnique[c(3,4,5,6,2,c(7:tl))]
  
  saveRDS(allsoUniqueP,file= file.path(output.dir,basename(paste(files[f],"_CollisionTable.rds",sep=""))))

  write.table(allsoUniqueP,file = file.path(output.dir,basename(paste(files[f],"_CollisionTable.txt",sep=""))),sep = "\t",row.names = F,quote = F,col.names = T)
}    
